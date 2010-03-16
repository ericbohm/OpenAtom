#include "arrayBuilder.h"
#include "orthoMap.h"

#include "main/TimeKeeper.h"
#include "main/CLA_Matrix.h"
#include "utility/MapFile.h"
#include "ortho.decl.h"

#include "ckmulticast.h"

namespace cp {
    namespace ortho {

/**
 * Create the map objects and also all the chare arrays needed for an Ortho instance
 */
CkArrayID ArrayBuilder::build(int nstates, PeListFactory getPeList, UberCollection thisInstance)
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

    PeList *availGlobR = getPeList();

    //-------------------------------------------------------------------------
    // Create maps for placing the Ortho chare array elements

    PeList *excludePes= new PeList(1);
    excludePes->TheList[0]=config.numPes;
    int nOrtho= (nstates/config.orthoGrainSize);
    nOrtho *= nOrtho;
    double Timer=CmiWallTimer();

    availGlobR->reset();
    #ifdef USE_INT_MAP
        OrthoImaptable[thisInstance.getPO()].buildMap(nstates/config.orthoGrainSize, nstates/config.orthoGrainSize);
    #endif

    int success = 0;
    if(config.loadMapFiles)
    {
        int size[2];
        size[0] = size[1] = nstates/config.orthoGrainSize;
        MapFile *mf = new MapFile("OrthoMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        #ifdef USE_INT_MAP
            success = mf->loadMap("OrthoMap", &OrthoImaptable[thisInstance.getPO()]);
        #endif
        delete mf;
    }
    PeList *avail= new PeList();
    if(success == 0)
    {
        #ifdef USE_INT_MAP
            OrthoMapTable Otable = OrthoMapTable(&OrthoImaptable[thisInstance.getPO()], avail, nstates, config.orthoGrainSize, &AsymScalcImaptable[thisInstance.getPO()], config.nchareG, config.numChunks, config.sGrainSize, excludePes);
        #endif
    }

    double newtime=CmiWallTimer();
    CkPrintf("OrthoMap created in %g\n\n", newtime-Timer);

    // If map files need to be dumped
    if(config.dumpMapFiles)
    {
        int size[2];
        size[0] = size[1] = nstates/config.orthoGrainSize;
        MapFile *mf = new MapFile("OrthoMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        #ifdef USE_INT_MAP
            mf->dumpMap(&OrthoImaptable[thisInstance.getPO()], thisInstance.getPO());
        #endif
        delete mf;
    }
    
    // If the map coordinates need to be dumped
    if(config.dumpMapCoordFiles)
    {
        int size[2];
        size[0] = size[1] = nstates/config.orthoGrainSize;
        MapFile *mf = new MapFile("OrthoMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
        #ifdef USE_INT_MAP
            mf->dumpMapCoords(&OrthoImaptable[thisInstance.getPO()], thisInstance.getPO());
        #endif
        delete mf;
    }

    // Create the ortho map group
    CProxy_OrthoMap orthoMap = CProxy_OrthoMap::ckNew(thisInstance);
    CkArrayOptions orthoOpts;
    orthoOpts.setMap(orthoMap);
    CProxy_Ortho orthoProxy = CProxy_Ortho::ckNew(orthoOpts);

    // Create maps for the Ortho helper chares
    if(config.useOrthoHelpers)
    {
        #ifdef USE_INT_MAP
            OrthoHelperImaptable[thisInstance.getPO()].buildMap(nstates/config.orthoGrainSize, nstates/config.orthoGrainSize);
        #endif
        double Timer=CmiWallTimer();
        success = 0;
        if(config.loadMapFiles)
        {
            int size[2];
            size[0] = size[1] = nstates/config.orthoGrainSize;
            MapFile *mf = new MapFile("OrthoHelperMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
            #ifdef USE_INT_MAP
                success = mf->loadMap("OrthoHelperMap", &OrthoHelperImaptable[thisInstance.getPO()]);
            #endif
            delete mf;
        }

        if(success == 0)
        {
            #ifdef USE_INT_MAP
                OrthoHelperMapTable OHtable = OrthoHelperMapTable(&OrthoHelperImaptable[thisInstance.getPO()], nstates, config.orthoGrainSize, &OrthoImaptable[thisInstance.getPO()], avail, excludePes);
            #endif
        }
        double newtime=CmiWallTimer();
        CkPrintf("OrthoHelperMap created in %g\n", newtime-Timer);
        CProxy_OrthoHelperMap orthoHMap = CProxy_OrthoHelperMap::ckNew(thisInstance);
        CkArrayOptions orthoHOpts;
        orthoHOpts.setMap(orthoHMap);
        orthoHelperProxy = CProxy_OrthoHelper::ckNew(orthoHOpts);
        helperAID = orthoHelperProxy.ckGetArrayID();

        if(config.dumpMapFiles)
        {
            int size[2];
            size[0] = size[1] = nstates/config.orthoGrainSize;
            MapFile *mf = new MapFile("OrthoHelperMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
            #ifdef USE_INT_MAP
                mf->dumpMap(&OrthoHelperImaptable[thisInstance.getPO()], thisInstance.getPO());
            #endif
            delete mf;
        }

        if(config.dumpMapCoordFiles)
        {
            int size[2];
            size[0] = size[1] = nstates/config.orthoGrainSize;
            MapFile *mf = new MapFile("OrthoHelperMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
            #ifdef USE_INT_MAP
                mf->dumpMapCoords(&OrthoHelperImaptable[thisInstance.getPO()], thisInstance.getPO());
            #endif
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
                    nstates, nstates, nstates,
                    config.orthoGrainSize, config.orthoGrainSize, config.orthoGrainSize,
                    1, 1, 1,
                    ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                    mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                   );

    if(config.useOrthoHelpers)
    {
        make_multiplier(&matA2, &matB2, &matC2,
                        orthoHelperProxy, orthoHelperProxy, orthoHelperProxy,
                        nstates, nstates, nstates,
                        config.orthoGrainSize, config.orthoGrainSize, config.orthoGrainSize,
                        1, 1, 1,
                        ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                        mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                       );
    }
    else
    {
        make_multiplier(&matA2, &matB2, &matC2,
                        orthoProxy, orthoProxy, orthoProxy,
                        nstates, nstates, nstates,
                        config.orthoGrainSize, config.orthoGrainSize, config.orthoGrainSize,
                        1, 1, 1,
                        ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                        mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                       );
    }

    make_multiplier(&matA3, &matB3, &matC3,
                    orthoProxy, orthoProxy, orthoProxy,
                    nstates, nstates, nstates,
                    config.orthoGrainSize, config.orthoGrainSize, config.orthoGrainSize,
                    1, 1, 1,
                    ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
                    mCastGID, MM_ALG_2D, config.gemmSplitOrtho
                   );

    // Register with the time keeper
    int timekeep=keeperRegister("Ortho S to T");
    int maxorthoindex=(nstates/config.orthoGrainSize-1);
    int maxorthostateindex=(nstates/config.orthoGrainSize-1) * config.orthoGrainSize;
    // Insert each element of the Ortho array
    for (int s1 = 0; s1 <= maxorthostateindex; s1 += config.orthoGrainSize)
        for (int s2 = 0; s2 <= maxorthostateindex; s2 += config.orthoGrainSize)
        {
            int indX = s1 / config.orthoGrainSize;
            int indY = s2 / config.orthoGrainSize;
            indX = (indX>maxorthoindex) ? maxorthoindex : indX;
            indY = (indY>maxorthoindex) ? maxorthoindex : indY;

            orthoProxy(indX, indY).insert(
                                                          config.orthoGrainSize, config.orthoGrainSize,
                                                          matA1, matB1, matC1,
                                                          matA2, matB2, matC2,
                                                          matA3, matB3, matC3,
                                                          helperAID,
                                                          timekeep, thisInstance,
                                                          orthoMcastGID, orthoRedGID);

            if(config.useOrthoHelpers)
            {
                CkCallback endOfStep2CB(CkIndex_Ortho::recvStep2(NULL), CkArrayIndex2D(indX,indY), orthoProxy);
                orthoHelperProxy(indX, indY).insert(
                                                          config.orthoGrainSize, config.orthoGrainSize,
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

