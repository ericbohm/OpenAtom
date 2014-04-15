#include "pcBuilder.h"
#include "pcMaps.h"
#include "utility/MapFile.h"
#include "ckPairCalculator.decl.h"
#include "InputDataHandler.h"

#include "ckmulticast.h"

extern Config config;
/*
 * @addtogroup PairCalculator
 *    @{
 */

namespace cp {
    namespace paircalc {

InstanceIDs Builder::build(const startup::PCMapConfig &mapCfg)
{
    traceRegisterUserEvent("calcpairDGEMM", 210);
    traceRegisterUserEvent("calcpairContrib", 220);
    traceRegisterUserEvent("multiplyResultDGEMM1", 230);
    traceRegisterUserEvent("multiplyResultDGEMM2", 240);
    traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

    createMap(mapCfg);
    createPairCalcs();

    pcHandle.mCastMgrGID = CProxy_CkMulticastMgr::ckNew(cfg.inputSpanningTreeFactor);

    #ifdef USE_COMLIB
        // Setup the appropriate multicast strategy
        Strategy *multistrat = new DirectMulticastStrategy();
        if(cfg.isSymmetric)
            mcastInstanceCP=ComlibRegister(multistrat);
        else
            mcastInstanceACP=ComlibRegister(multistrat);
    #endif

    return pcHandle;
}


// Hacky way to allow PC builders of each instance to access the maptable for instance 0 pc
namespace impl { MapType4 *dirtyGlobalMapTable4PCsym, *dirtyGlobalMapTable4PCasym; }

/**
 * Create the map for placing the paircalculator chare array elements. Also perform other housekeeping chores like dumping the maps to files etc.
 */
void Builder::createMap(const startup::PCMapConfig &mapCfg)
{
    bool maptype = cfg.isSymmetric;
    int achunks = config.numChunksAsym;
    if(cfg.isSymmetric && cfg.arePhantomsOn)
    { // evil trickery to use asym map code for phantom sym
        maptype=false;
        achunks=config.numChunksSym;
    }

    MapType4 **inst0MapTable = cfg.isSymmetric ? &impl::dirtyGlobalMapTable4PCsym : &impl::dirtyGlobalMapTable4PCasym;
    // Generate a map name
    std::string mapName = cfg.isSymmetric ? "SymScalcMap" : "AsymScalcMap";
    /// Get an appropriately constructed PeList from the supplied factory functor
    PeList *availGlobG = mapCfg.getPeList();
    availGlobG->reset();

    // Compute num PEs along the states dimension of GSpace
    int pl = cfg.numStates / config.Gstates_per_pe;
    // Compute num PEs along the planes dimension of GSpace
    int pm = config.numPesPerInstance / pl;
    // Compute the num of GSpace planes per PE
    int planes_per_pe = cfg.numPlanes / pm;

    int size[4];
    size[0] = cfg.numPlanes;
    size[1] = cfg.numStates/cfg.grainSize;
    size[2] = cfg.numStates/cfg.grainSize;
    size[3] = achunks;
    //-------------------------------------------------------------
    // Populate maptable for the PC array

    // Start timing the map creation
    double mapCreationTime = CmiWallTimer();

    MapType4 mapTable;
    if (cfg.instanceIndex == 0)
    {
        mapTable.buildMap(cfg.numPlanes, cfg.numStates/cfg.grainSize, cfg.numStates/cfg.grainSize, achunks, cfg.grainSize);

        int success = 0;
        if(config.loadMapFiles)
        {
            MapFile *mf = new MapFile(mapName.c_str(), 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, cfg.grainSize);
            success = mf->loadMap(mapName.c_str(), &mapTable);
            delete mf;
        }

        // If loading the map from a file failed, create a maptable
        if(success == 0)
        {
            SCalcMapTable symTable = SCalcMapTable(&mapTable, availGlobG, cfg.numStates, cfg.numPlanes, cfg.grainSize, maptype, config.scalc_per_plane,
                            planes_per_pe, achunks, config.numChunksSym, mapCfg.gSpaceMap, config.useCuboidMap, config.useCentroidMap, mapCfg.boxSize);
        }

        // Save a globally visible handle to the mapTable that builders of other PC instances can access
        *inst0MapTable = new MapType4(mapTable);
    }
    // else, this is not the first instance. Simply translate the 0th instance
    else
    {
        int x = mapCfg.mapOffset.getx();
        int y = mapCfg.mapOffset.gety();
        int z = mapCfg.mapOffset.getz();
        if((CkNumPes()==1) && !mapCfg.isTorusFake)
            mapTable = *(new MapType4(**inst0MapTable)); ///< Lazy way to avoid having to write an assignment operator for MapType4
        else
            mapTable.translate(*inst0MapTable, x, y, z, mapCfg.isTorusMap);
    }

    /// Create a map group that will read and use this map table
    CProxy_SCalcMap pcMapGrp = CProxy_SCalcMap::ckNew(mapTable);

    mapCreationTime = CmiWallTimer() - mapCreationTime;
    CkPrintf("PairCalculator[%dx%dx%dx%d,%d] map created in %g\n", size[0], size[1], size[2], cfg.numChunks, cfg.isSymmetric, mapCreationTime);

    // If the user wants map dumps
    if(config.dumpMapFiles)
    {
        MapFile *mf = new MapFile(mapName.c_str(), 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, cfg.grainSize);
        mf->dumpMap(&mapTable, cfg.instanceIndex);
        delete mf;
    }

    // If the user wants map coordinate dumps
    if(config.dumpMapCoordFiles)
    {
        MapFile *mf = new MapFile((mapName+"_coord").c_str(), 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, cfg.grainSize);
        mf->dumpMapCoords(&mapTable, cfg.instanceIndex);
        delete mf;
    }

    // Record the group that will provide the procNum mapping function
    pcHandle.mapperGID  = pcMapGrp.ckGetGroupID();
    delete availGlobG;
}



/*
void Builder::createInputHandler(const pcConfig &cfg)
{
    CkArrayOptions handlerOpts;
    /// Create an empty input handler chare array that will accept all incoming messages from GSpace
    handlerOpts.bindTo(pairCalculatorProxy);
    CProxy_InputDataHandler<CollatorType,CollatorType> inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);
}
*/



void Builder::createPairCalcs()
{
    CkArrayOptions paircalcOpts,handlerOpts;
    CProxy_PairCalculator pairCalculatorProxy;
    CProxy_InputDataHandler<CollatorType,CollatorType> inputHandlerProxy;

    // Create an empty array but specify element locations using the map
    paircalcOpts.setMap(pcHandle.mapperGID);
    paircalcOpts.setAnytimeMigration(false);
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(inputHandlerProxy, cfg, paircalcOpts);

    #ifdef DEBUG_CP_PAIRCALC_CREATION
        CkPrintf("Builder: Creating a%s paircalc instance\n", (cfg.isSymmetric?" symmetric":"n asymmetric") );
    #endif
    #ifdef CP_PAIRCALC_USES_COMPLEX_MATH
        CkPrintf("Builder: Creating paircalcs for instance %d that use complex math\n", cfg.instanceIndex);
    #else
        CkPrintf("Builder: Creating paircalcs for instance %d that do not use complex math\n", cfg.instanceIndex);
    #endif

    /// Create an empty input handler chare array that will accept all incoming messages from GSpace
    handlerOpts.bindTo(pairCalculatorProxy);
    inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);

    // Initialize my set of array / group IDs
    pcHandle.pcAID = pairCalculatorProxy.ckGetArrayID();
    pcHandle.handlerAID = inputHandlerProxy.ckGetArrayID();

    // Compute the max value of the state dimension indices of the paircalc array
    int pcMaxStateDimIndex = (cfg.numStates / cfg.grainSize - 1) * cfg.grainSize;
    // Populate the sparse, 4D paircalc array
    for(int numX = 0; numX < cfg.numPlanes; numX ++)
    {
        for (int s1 = 0; s1 <= pcMaxStateDimIndex; s1 += cfg.grainSize)
        {
            // Make the symmetric array a triangular prism of chares only if phantoms are not needed
            int s2start = (cfg.isSymmetric && !cfg.arePhantomsOn) ? s1 : 0;
            for (int s2 = s2start; s2 <= pcMaxStateDimIndex; s2 += cfg.grainSize)
            {
                for (int c = 0; c < cfg.numChunks; c++)
                {
                    #ifdef DEBUG_CP_PAIRCALC_CREATION
                        CkPrintf("\tInserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,cfg.isSymmetric);
                    #endif
                    pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, cfg);
                }
            }
        }
    }

    /// Notify the runtime that we're done inserting all the PC elements
    pairCalculatorProxy.doneInserting();
}




    } // end namespace paircalc
} // end namespace cp
/*@}*/
#include "pcMaps.def.h"

