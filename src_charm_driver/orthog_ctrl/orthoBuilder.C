#include "orthoBuilder.h"
#include "orthoMap.h"
#include "CLA_Matrix.h"

#include "paircalc/pcMaps.h"
#include "main/TimeKeeper.h"
#include "utility/MapFile.h"
#include "ortho.decl.h"

#include "ckmulticast.h"

namespace cp {
    namespace ortho {

// Hacky way to allow Ortho builders of each instance to access the maptable for instance 0 pc
namespace impl { MapType2 *dirtyGlobalMapTable4Ortho, *dirtyGlobalMapTable4OrthoHelper; }

/**
 * Create the map objects and also all the chare arrays needed for an Ortho instance
 */
CkArrayID Builder::build(cp::paircalc::InstanceIDs &asymmHandle, const startup::PCMapConfig mapCfg)
{
    CkPrintf("Building Ortho Chares\n");

    // For using the multicast library :  Set some reduction clients
    /// MulticastMgr that handles CLA_Matrix collectives
    CkGroupID mCastGID;
    /// MulticastMgr that handles Ortho --> PC mcasts
    CkGroupID orthoMcastGID;
    /// MulticastMgr that handles PC --> Ortho redns
    CkGroupID orthoRedGID;
    /// The step 2 helper chare array proxy
    CProxy_OrthoHelper orthoHelperProxy;
    /// The step 2 helper chare AID to pass to ortho
    CkArrayID helperAID;
    /// The step 2 helper map logic object
    MapType2 helperMapTable;

    PeList *availGlobR = mapCfg.getPeList();

    //-------------------------------------------------------------------------
    // Create maps for placing the Ortho chare array elements

    PeList *excludePes= new PeList(1);
    excludePes->TheList[0]=config.numPes;
    int nOrtho= (cfg.numStates/cfg.grainSize);
    nOrtho *= nOrtho;
    double Timer=CmiWallTimer();

    availGlobR->reset();
    PeList *avail= new PeList();
    MapType2 orthoMapTable;
    if (cfg.instanceIndex == 0)
    {
        orthoMapTable.buildMap(cfg.numStates/cfg.grainSize, cfg.numStates/cfg.grainSize);

        int success = 0;
        if(config.loadMapFiles)
        {
            int size[2];
            size[0] = size[1] = cfg.numStates/cfg.grainSize;
            MapFile *mf = new MapFile("OrthoMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
            success = mf->loadMap("OrthoMap", &orthoMapTable);
            delete mf;
        }
        if(success == 0)
        {
            SCalcMap *asymmMap = CProxy_SCalcMap(asymmHandle.mapperGID).ckLocalBranch();
            MapType4 *maptable = asymmMap->getMapTable();
            OrthoMapTable Otable = OrthoMapTable(&orthoMapTable, avail, cfg.numStates, cfg.grainSize, maptable, config.nchareG, config.numChunks, config.sGrainSize, excludePes);
        }

        // Save a globally visible handle to the mapTable that builders of other ortho instances can access
        impl::dirtyGlobalMapTable4Ortho = new MapType2(orthoMapTable);
    }
    // else, simply translate the instance 0 ortho map
    else
    {
        int x = mapCfg.mapOffset.getx();
        int y = mapCfg.mapOffset.gety();
        int z = mapCfg.mapOffset.getz();
        orthoMapTable.translate(impl::dirtyGlobalMapTable4Ortho, x, y, z, mapCfg.isTorusMap);
    }

    double newtime=CmiWallTimer();
    CkPrintf("OrthoMap created in %g\n\n", newtime-Timer);

    // If map files need to be dumped
    if(config.dumpMapFiles)
    {
        int size[2];
        size[0] = size[1] = cfg.numStates/cfg.grainSize;
        MapFile *mf = new MapFile("OrthoMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        mf->dumpMap(&orthoMapTable, cfg.instanceIndex);
        delete mf;
    }

    // If the map coordinates need to be dumped
    if(config.dumpMapCoordFiles)
    {
        int size[2];
        size[0] = size[1] = cfg.numStates/cfg.grainSize;
        MapFile *mf = new MapFile("OrthoMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        mf->dumpMapCoords(&orthoMapTable, cfg.instanceIndex);
        delete mf;
    }

    // Create the ortho map group
    CProxy_OrthoMap orthoMap = CProxy_OrthoMap::ckNew(orthoMapTable);
    CkArrayOptions orthoOpts;
    orthoOpts.setMap(orthoMap);
    CProxy_Ortho orthoProxy = CProxy_Ortho::ckNew(orthoOpts);



    // Create maps for the Ortho helper chares
    if(config.useOrthoHelpers)
    {
        double Timer=CmiWallTimer();

        if (cfg.instanceIndex == 0)
        {
            helperMapTable.buildMap(cfg.numStates/cfg.grainSize, cfg.numStates/cfg.grainSize);
            int success = 0;
            if(config.loadMapFiles)
            {
                int size[2];
                size[0] = size[1] = cfg.numStates/cfg.grainSize;
                MapFile *mf = new MapFile("OrthoHelperMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
                success = mf->loadMap("OrthoHelperMap", &helperMapTable);
                delete mf;
            }

            if(success == 0)
            {
                OrthoHelperMapTable OHtable = OrthoHelperMapTable(&helperMapTable, cfg.numStates, cfg.grainSize, &orthoMapTable, avail, excludePes);
            }
            // Save a globally visible handle to the mapTable that builders of other ortho instances can access
            impl::dirtyGlobalMapTable4OrthoHelper = new MapType2(orthoMapTable);
        }
        else
        {
            int x = mapCfg.mapOffset.getx();
            int y = mapCfg.mapOffset.gety();
            int z = mapCfg.mapOffset.getz();
            helperMapTable.translate(impl::dirtyGlobalMapTable4OrthoHelper, x, y, z, mapCfg.isTorusMap);
        }

        double newtime=CmiWallTimer();
        CkPrintf("OrthoHelperMap created in %g\n", newtime-Timer);
        CProxy_OrthoHelperMap orthoHMap = CProxy_OrthoHelperMap::ckNew(helperMapTable);
        CkArrayOptions orthoHOpts;
        orthoHOpts.setMap(orthoHMap);
        orthoHelperProxy = CProxy_OrthoHelper::ckNew(orthoHOpts);
        helperAID = orthoHelperProxy.ckGetArrayID();

        if(config.dumpMapFiles)
        {
            int size[2];
            size[0] = size[1] = cfg.numStates/cfg.grainSize;
            MapFile *mf = new MapFile("OrthoHelperMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
            mf->dumpMap(&helperMapTable, cfg.instanceIndex);
            delete mf;
        }

        if(config.dumpMapCoordFiles)
        {
            int size[2];
            size[0] = size[1] = cfg.numStates/cfg.grainSize;
            MapFile *mf = new MapFile("OrthoHelperMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
            mf->dumpMapCoords(&helperMapTable, cfg.instanceIndex);
            delete mf;
        }
    }

    //-------------------------------------------------------------------------
    // Delegate collectives within the ortho array
    #ifdef USE_COMLIB
        Strategy *multistrat = new DirectMulticastStrategy();
        orthoInstance=ComlibRegister(multistrat);
        //  ComlibAssociateProxy(orthoInstance, orthoProxy);
    #endif

    // Set the root of array reductions within Ortho
    CkCallback ocb= CkCallback(CkIndex_Ortho::collect_error(NULL), orthoProxy(0, 0));
    orthoProxy.ckSetReductionClient(&ocb);

    // Each multiplier calls ortho back to notify that its ready
    CkCallback ortho_ready_cb = CkCallback(CkIndex_Ortho::all_ready(), orthoProxy(0, 0));
    // The multicast group that will handle CLA_Matrix collectives 
    mCastGID = CProxy_CkMulticastMgr::ckNew(config.numMulticastMsgs);

    /// Create multicast manager groups for Ortho to use
    orthoMcastGID = (CProxy_CkMulticastMgr::ckNew(config.OrthoMcastSpanFactor));
    orthoRedGID   = (CProxy_CkMulticastMgr::ckNew(config.OrthoRedSpanFactor));

    //-------------------------------------------------------------------------
    // Create matrix multiplication objects

    // extra triangle ortho elements are really a waste of our time
    // and resources, but we don't have a triangular solver for
    // inv_square, so we'll just make do.
    // They need to exist solely so that the inv_sq method can work.
    // So we need to copy their mirror elements data into them.
    // then when complete they need to know not to call finishpaircalc.
    // Because their redundant data has nowhere to go.
    // We've made use of them anyway to handle: lambda reduction, the
    // gamma multiply, and communication balancing for phantoms, so they
    // aren't completely horrible.

    CLA_Matrix_interface matA1, matB1, matC1;
    CLA_Matrix_interface matA2, matB2, matC2;
    CLA_Matrix_interface matA3, matB3, matC3;

    make_multiplier(&matA1, &matB1, &matC1,
                    orthoProxy, orthoProxy, orthoProxy,
                    cfg.numStates, cfg.numStates, cfg.numStates,
                    cfg.grainSize, cfg.grainSize, cfg.grainSize,
                    1, 1, 1,
                    ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                    mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                   );

    if(config.useOrthoHelpers)
    {
        make_multiplier(&matA2, &matB2, &matC2,
                        orthoHelperProxy, orthoHelperProxy, orthoHelperProxy,
                        cfg.numStates, cfg.numStates, cfg.numStates,
                        cfg.grainSize, cfg.grainSize, cfg.grainSize,
                        1, 1, 1,
                        ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                        mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                       );
    }
    else
    {
        make_multiplier(&matA2, &matB2, &matC2,
                        orthoProxy, orthoProxy, orthoProxy,
                        cfg.numStates, cfg.numStates, cfg.numStates,
                        cfg.grainSize, cfg.grainSize, cfg.grainSize,
                        1, 1, 1,
                        ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                        mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                       );
    }

    make_multiplier(&matA3, &matB3, &matC3,
                    orthoProxy, orthoProxy, orthoProxy,
                    cfg.numStates, cfg.numStates, cfg.numStates,
                    cfg.grainSize, cfg.grainSize, cfg.grainSize,
                    1, 1, 1,
                    ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                    mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                   );

    // Register with the time keeper
    int timekeep=keeperRegister("Ortho S to T");
    int maxorthoindex=(cfg.numStates/cfg.grainSize-1);
    int maxorthostateindex=(cfg.numStates/cfg.grainSize-1) * cfg.grainSize;
    // Insert each element of the Ortho array
    for (int s1 = 0; s1 <= maxorthostateindex; s1 += cfg.grainSize)
        for (int s2 = 0; s2 <= maxorthostateindex; s2 += cfg.grainSize)
        {
            int indX = s1 / cfg.grainSize;
            int indY = s2 / cfg.grainSize;
            indX = (indX>maxorthoindex) ? maxorthoindex : indX;
            indY = (indY>maxorthoindex) ? maxorthoindex : indY;

            orthoProxy(indX, indY).insert(
                                                          cfg.grainSize, cfg.grainSize,
                                                          matA1, matB1, matC1,
                                                          matA2, matB2, matC2,
                                                          matA3, matB3, matC3,
                                                          cfg,
                                                          helperAID,
                                                          timekeep,
                                                          orthoMcastGID, orthoRedGID);

            if(config.useOrthoHelpers)
            {
                CkCallback endOfStep2CB(CkIndex_Ortho::recvStep2(NULL), CkArrayIndex2D(indX,indY), orthoProxy);
                orthoHelperProxy(indX, indY).insert(
                                                          cfg.grainSize, cfg.grainSize,
                                                          matA2, matB2, matC2,
                                                          endOfStep2CB);
            }
        }
    // Notify that you're done inserting the elements
    orthoProxy.doneInserting();
    if(config.useOrthoHelpers)
        orthoHelperProxy.doneInserting();

    delete avail;
    delete excludePes;

    return orthoProxy.ckGetArrayID();
}

    } // end namespace ortho
} // end namespace cp

