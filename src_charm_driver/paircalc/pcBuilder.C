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
      if (cfg.instanceIndex == 0 || config.simpleTopo)
      {
        /// Get an appropriately constructed PeList from the supplied factory functor
        PeList *availGlobG = mapCfg.getPeList();
        availGlobG->reset();
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
        delete availGlobG;
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
    }

    void Builder::createPairCalcs()
    {
      CkArrayOptions paircalcOpts, handlerOpts;
      CProxy_PairCalculator pairCalculatorProxy;
      CProxy_InputDataHandler<CollatorType,CollatorType> inputHandlerProxy;

#ifdef DEBUG_CP_PAIRCALC_CREATION
      CkPrintf("Builder: Creating a%s paircalc instance.\n",
          (cfg.isSymmetric?" symmetric":"n asymmetric") );
#ifdef CP_PAIRCALC_USES_COMPLEX_MATH
      CkPrintf("Builder: Creating paircalcs for instance %d that use complex math\n", cfg.instanceIndex);
#else
      CkPrintf("Builder: Creating paircalcs for instance %d that do not use complex math\n", cfg.instanceIndex);
#endif
#endif

      // Set up indices for start, end, step, and bounds. This will determine
      // how many chares to create and at which indices to create them when we
      // call ckNew().
      short numPlanes = (short)cfg.numPlanes;
      short numStates = (short)cfg.numStates;
      short grainSize = (short)cfg.grainSize;
      short numChunks = (short)cfg.numChunks;
      // The state indices are set up in such a way that if numStates isn't
      // divisible by grainSize, the remainder is stitched to the last chare.
      short stateEnd  = (numStates / grainSize) * grainSize;
      // Because we always start the 3rd index at zero, we must have phantoms
      // on at all times. The current bulk construction in Charm++ is not
      // flexible enough to handle creating the lower triangle of an array.
      CkIndex4D start = {0,0,0,0};
      CkIndex4D end = {numPlanes,stateEnd,stateEnd,numChunks};
      CkIndex4D step = {1,grainSize,grainSize,1};

#ifdef DEBUG_CP_PAIRCALC_CREATION
      CkPrintf("Builder: numPlanes = %d, numStates = %d, grainSize = %d, numChunks = %d\n",
          numPlanes, numStates, grainSize, numChunks);
      CkPrintf("Builder: Start = (%d, %d, %d, %d).\n",
          start.w, start.x, start.y, start.z);
      CkPrintf("Builder: End = (%d, %d, %d, %d).\n",
          end.w, end.x, end.y, end.z);
      CkPrintf("Builder: Step = (%d, %d, %d, %d).\n",
          step.w, step.x, step.y, step.z);
#endif

      // Set up paircalc to use our custom mapping, and the indices set above.
      paircalcOpts.setMap(pcHandle.mapperGID);
      paircalcOpts.setAnytimeMigration(false);
      paircalcOpts.setStaticInsertion(true);
      paircalcOpts.setStart(CkArrayIndex4D(start)).setEnd(CkArrayIndex4D(end)).setStep(CkArrayIndex4D(step));
      paircalcOpts.setBounds(end.w, end.x, end.y, end.z);

      pairCalculatorProxy = CProxy_PairCalculator::ckNew(cfg, paircalcOpts);

      // The handler will be bound to the paircalc and use the same options.
      handlerOpts.bindTo(pairCalculatorProxy);
      handlerOpts.setAnytimeMigration(false);
      handlerOpts.setStaticInsertion(true);
      handlerOpts.setStart(CkArrayIndex4D(start)).setEnd(CkArrayIndex4D(end)).setStep(CkArrayIndex4D(step));
      handlerOpts.setBounds(end.w, end.x, end.y, end.z);

      inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType>::ckNew(pairCalculatorProxy,handlerOpts);

      // Initialize my set of array / group IDs.
      pcHandle.pcAID = pairCalculatorProxy.ckGetArrayID();
      pcHandle.handlerAID = inputHandlerProxy.ckGetArrayID();
    }

  } // end namespace paircalc
} // end namespace cp
/*@}*/
#include "pcMaps.def.h"

