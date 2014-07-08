#ifdef _PAIRCALC_STARTUP_H
#define _PAIRCALC_STARTUP_H

     /* choose whether ortho should use local callback */
     Ortho_use_local_cb = true;


    // Stuff it with the actual configurations
    cfgSymmPC.isDynamics         = (sim->cp_min_opt==1)? false: true;
    cfgSymmPC.useComplexMath     = false;

    cfgSymmPC.numPlanes          = config.nchareG;
    cfgSymmPC.numStates          = nstates;
    cfgSymmPC.grainSize          = config.sGrainSize;
    cfgSymmPC.orthoGrainSize     = config.orthoGrainSize;

    cfgSymmPC.conserveMemory     = config.conserveMemory;
    cfgSymmPC.isLBon             = config.lbpaircalc;

    cfgSymmPC.areBWTilesCollected= config.PCCollectTiles;
    cfgSymmPC.isBWstreaming      = config.PCstreamBWout;
    cfgSymmPC.isBWbarriered      = config.useBWBarrier;
    cfgSymmPC.shouldDelayBWsend  = config.PCdelayBWSend;
    cfgSymmPC.isInputMulticast   = !config.usePairDirectSend;
    cfgSymmPC.isOutputReduced    = !config.gSpaceSum;
    cfgSymmPC.inputSpanningTreeFactor = config.PCSpanFactor;

    cfgSymmPC.gemmSplitFWk       = config.gemmSplitFWk;
    cfgSymmPC.gemmSplitFWm       = config.gemmSplitFWm;
    cfgSymmPC.gemmSplitBW        = config.gemmSplitBW;



    // Configurations specific to the symmetric PC instance
    cfgSymmPC.isSymmetric        = true;
    cfgSymmPC.arePhantomsOn      = config.phantomSym;
    cfgSymmPC.numChunks          = config.numChunksSym;
    cfgSymmPC.isDoublePackOn     = doublePack;
    cfgSymmPC.inputMsgPriority   = config.psipriority;
    cfgSymmPC.resultMsgPriority  = config.gsfftpriority;

    // Configurations specific to the asymmetric PC instance
    cfgAsymmPC.isSymmetric        = false;
    cfgAsymmPC.arePhantomsOn      = false;
    cfgAsymmPC.numChunks          = config.numChunksAsym;
    cfgAsymmPC.isDoublePackOn     = 0;
    cfgAsymmPC.inputMsgPriority   = config.lambdapriority;
    cfgAsymmPC.resultMsgPriority  = config.lambdapriority+2;

    // Configure the GSpace entry methods that the PCs will callback
    if(cfgSymmPC.isOutputReduced)
    {
        cfgSymmPC.gSpaceEP        = CkIndex_CP_State_GSpacePlane::acceptNewPsi ((CkReductionMsg*)NULL);
        cfgSymmPC.PsiVEP          = CkIndex_CP_State_GSpacePlane::acceptNewPsiV((CkReductionMsg*)NULL);
    }
    else
    {
        cfgSymmPC.gSpaceEP        = CkIndex_CP_State_GSpacePlane::acceptNewPsi ((partialResultMsg*)NULL);
        cfgSymmPC.PsiVEP          = CkIndex_CP_State_GSpacePlane::acceptNewPsiV((partialResultMsg*)NULL);
    }

    if(cfgAsymmPC.isOutputReduced)
    {
        cfgAsymmPC.gSpaceEP       = CkIndex_CP_State_GSpacePlane::acceptLambda ((CkReductionMsg*)NULL);
        cfgAsymmPC.PsiVEP         = 0;
    }
    else
    {
        cfgAsymmPC.gSpaceEP       = CkIndex_CP_State_GSpacePlane::acceptLambda ((partialResultMsg*)NULL);
        cfgAsymmPC.PsiVEP         = 0;
    }

#ifdef _CP_SUBSTEP_TIMING_
    //symmetric AKA Psi
    cfgSymmPC.forwardTimerID      = keeperRegister("Sym Forward");
    cfgSymmPC.backwardTimerID     = keeperRegister("Sym Backward");
    cfgSymmPC.beginTimerCB        = CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
    cfgSymmPC.endTimerCB          = CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
    //asymmetric AKA Lambda AKA Gamma
    cfgAsymmPC.forwardTimerID     = keeperRegister("Asym Forward");
    cfgAsymmPC.backwardTimerID    = keeperRegister("Asym Backward");
    cfgAsymmPC.beginTimerCB       = CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
    cfgAsymmPC.endTimerCB         = CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
#endif

#endif
