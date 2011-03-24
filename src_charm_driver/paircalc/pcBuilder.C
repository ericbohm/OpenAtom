#include "pcBuilder.h"
#include "pcMaps.h"
#include "utility/MapFile.h"
#include "ckPairCalculator.decl.h"
#include "InputDataHandler.h"

#include "ckmulticast.h"

extern Config config;

namespace cp {
    namespace paircalc {

InstanceIDs Builder::build(const int boxSize, PeListFactory getPeList, MapType2 *gSpaceMap)
{
    traceRegisterUserEvent("calcpairDGEMM", 210);
    traceRegisterUserEvent("calcpairContrib", 220);
    traceRegisterUserEvent("multiplyResultDGEMM1", 230);
    traceRegisterUserEvent("multiplyResultDGEMM2", 240);
    traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

    createMap(boxSize, getPeList, gSpaceMap);
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




/**
 * Create the map for placing the paircalculator chare array elements. Also perform other housekeeping chores like dumping the maps to files etc.
 */
void Builder::createMap(const int boxSize, PeListFactory getPeList, MapType2 *gSpaceMap)
{
    bool maptype = cfg.isSymmetric;
    int achunks = config.numChunksAsym;
    if(cfg.isSymmetric && cfg.arePhantomsOn)
    { // evil trickery to use asym map code for phantom sym
        maptype=false;
        achunks=config.numChunksSym;
    }

    // Generate a map name
    std::string mapName = cfg.isSymmetric ? "SymScalcMap" : "AsymScalcMap";
    /// Get an appropriately constructed PeList from the supplied factory functor
    PeList *availGlobG = getPeList();
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
                        planes_per_pe, achunks, config.numChunksSym, gSpaceMap, config.useCuboidMap, config.useCentroidMap, boxSize);
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

    int cores_per_node = CkNumPes()/CmiNumPhysicalNodes();

    // Compute the max value of the state dimension indices of the paircalc array
    int pcMaxStateDimIndex = (cfg.numStates / cfg.grainSize - 1) * cfg.grainSize;
    int size = cfg.numPlanes*cfg.numStates*cfg.numStates*cfg.numChunks;

    CProxy_NodeMapPCArray pcCrayMap = CProxy_NodeMapPCArray::ckNew(cfg.numPlanes,cfg.numStates,cfg.numChunks,cfg.grainSize,(cfg.isSymmetric && !cfg.arePhantomsOn),cores_per_node,0,CmiNumPhysicalNodes(),0);

    // Create an empty array but specify element locations using the map
    paircalcOpts.setMap(pcCrayMap);
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(inputHandlerProxy, cfg, paircalcOpts);

    #ifdef DEBUG_CP_PAIRCALC_CREATION
        CkPrintf("Builder: Creating a%s paircalc instance\n", (cfg.isSymmetric?" symmetric":"n asymmetric") );
    #endif

    /// Create an empty input handler chare array that will accept all incoming messages from GSpace
    handlerOpts.bindTo(pairCalculatorProxy);
    inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);

    // Initialize my set of array / group IDs
    pcHandle.pcAID = pairCalculatorProxy.ckGetArrayID();
    pcHandle.handlerAID = inputHandlerProxy.ckGetArrayID();

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

#include "pcMaps.def.h"

