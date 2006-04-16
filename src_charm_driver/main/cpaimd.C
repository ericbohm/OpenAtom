//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file cpaimd.C
 *                         cpaimd-charm-driver
 *    Software developed by the Parallel Programing Laboratory, UIUC.
 *    in collaboration with IBM and NYU.
 *   
 *    This file contains cpaimd-charm-driver main. It creates and 
 *    initializes all the arrays and libraries. 
 */      
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \mainpage
<TABLE width="100%" border="0" cellspacing="0" cellpadding="5">
<TR bgcolor="#1C097D"><TD><font color="#FFFFFF" size="4"><table width=100%><tr><td><font size=6 color="#FFFFFF">Ab Initio Molecular Dynamics
</font></td><td align=right><img src="../images/image.jpg"></td></tr></table></FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><br>Many important problems in material science, chemistry, solid-state physics, and biophysics require a modeling
 approach based on fundamental quantum mechanical principles. A particular approach that has proved to be relatively
 efficient and useful is Car-Parrinello ab initio molecular dynamics (CPAIMD). Parallelization of this approach beyond a
 few hundred processors is challenging, due to the complex dependencies among various subcomputations, which lead to
 complex communication optimization and load balancing problems. We are parallelizing CPAIMD using Charm++. The computation is modeled using 
 a large number of virtual processors, which are mapped flexibly to
 available processors with assistance from the Charm++ runtime system.
 <br>
 <br>
 This project is a large NSF funded collaboration involving us (PPL : Laxmikant Kale)
 and Drs. Roberto Car, Michael Klein, Glenn Martyna, Mark Tuckerman, Nick Nystrom and Josep Torrellas.
 <br>
 
 <center>
 <img src="../images/Slide1.JPG">
 </center>

 
 
<br>&nbsp;</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">Software</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD>Currently we have a Charm++ implementation of the core of the CP method. You
 can check out the latest build using CVS. The code is available under the
 module name "new_leanCP". Charm++ and <a href="http://www.fftw.org">FFTW</a>
 are required to run the code.   
</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">People</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><UL>
<LI><A href="mailto:ebohm AT uiuc.edu
">Eric Bohm</A></LI>
<LI><A href="mailto:yanshi AT uiuc.edu
">Yan Shi</A></LI>

<LI><A href="mailto:kale AT cs.uiuc.edu
">L. V. Kale</A></LI>
<LI><A href="mailto:kunzman2 AT uiuc.edu
">David Kunzman</A></LI>
<LI><A href="mailto:skumar2 AT uiuc.edu
">Sameer Kumar</A></LI>

</UL>
</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">Papers</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><UL>
<LI><A href="/papers/CpaimdTR03.shtml">03-06
</A>&nbsp;&nbsp;&nbsp;Ramkumar Vadali, L. V. Kale, Glenn Martyna, Mark Tuckerman,&nbsp;&nbsp;<B>Scalable Parallelization of Ab Initio Molecular Dynamics</B>,
&nbsp;&nbsp;<I>Technical Report, communicated to SC 2003</I></LI>

</UL>
</TD></TR>
<TR bgcolor="#1C097D"><TD><FONT color="#FFFFFF" size="4">Related Links</FONT></TD></TR>
<TR bgcolor="#FFFFFF"><TD><UL>
<LI><A href="http://homepages.nyu.edu/~mt33/PINY_MD/PINY.html">PINY MD</A></LI>
<LI><a href="http://charm.cs.uiuc.edu/presentations/cpaimd03/cpaimd.ppt">Presentation on latest state of the code</a></LI>

</UL>
</TD></TR>
<tr bgcolor=#1C097D><td><table width=100%><tr bgcolor=#1C097D>
<td><font color=white>This page maintained by  <A href="mailto:ebohm AT uiuc.edu
">Eric Bohm</A>.</font></td>

<td align=right><font color=white>Back to the <A HREF="/research" onclick="document.research.submit();return false;">PPL Research Page</A></font></td>
</tr></table></td></tr>
</TABLE>
</td></tr></table>
 *
 * end of mainpage
 */


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================ 
#include <math.h>
#include "charm++.h"
#include "ckarray.h"
#include "util.h"
//============================================================================
#include "cpaimd.h"
#include "groups.h"
#include "ortho.h"
#include "sim_subroutines.h"
#include "StructFactorCache.h"
#include "StructureFactor.h"
#include "CP_State_Plane.h"
#include "MeshStreamingStrategy.h"
#include "MultiRingMulticast.h"
#ifdef USE_TOPOMAP
#include "bgltorus.h"
#include "FindProcessor.h"
#endif
//============================================================================
#include "../include/CPcharmParaInfo.h"
#include "../../src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "../../src_piny_physics_v1.0/include/charm_defs/Interface_ctrl.decl.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h"
//============================================================================



//============================================================================
/** 
 * \defgroup  piny_vars Defining all Charm++ readonly variables for PINY physics 
 *
 */
//============================================================================
/* @{ */
extern MDINTEGRATE  readonly_mdintegrate;
extern MDATOMS      readonly_mdatoms;
extern MDINTER      readonly_mdinter;
extern MDINTRA      readonly_mdintra;
extern GENERAL_DATA readonly_general_data;
extern CP           readonly_cp; 
/* @} */
//============================================================================


//============================================================================
/**
 * \defgroup proxy_comlib_vars
 * Defining all the Charm++ readonly variables, which include proxies
 * to access the arrays and groups and the Communication Library
 * handles.
 */
//============================================================================

Config config;
PairCalcID pairCalcID1;
PairCalcID pairCalcID2;
CProxy_main mainProxy;
CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
CProxy_CP_State_ParticlePlane particlePlaneProxy;
CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
CProxy_CP_Rho_GSpacePlane rhoGProxy;
CProxy_CP_Rho_GHartExt rhoGHartExtProxy;
CProxy_Ortho orthoProxy;
CProxy_CPcharmParaInfoGrp scProxy;
CProxy_AtomsGrp atomsGrpProxy;
CProxy_EnergyGroup egroupProxy;
CProxy_FFTcache fftCacheProxy;
CProxy_StructFactCache sfCacheProxy;
CProxy_StructureFactor sfCompProxy;

//============================================================================

//============================================================================
int nstates;              // readonly globals
int sizeX;
int nchareG;
int Ortho_UE_step2;
int Ortho_UE_step3;
int Ortho_UE_error;
bool Ortho_use_local_cb;
int done_init=0;
int planes_per_pe;
#ifdef USE_TOPOMAP
CkHashtableT<intdual, int> *maptable;
#if CMK_VERSION_BLUEGENE
BGLTorusManager *bgltm;
#endif
#endif
//============================================================================


//============================================================================
// For using the multicast library :  Set some reduction clients

CkGroupID mCastGrpId; 
// For using the communication library
ComlibInstanceHandle commGHartInstance;
ComlibInstanceHandle commGInstance0;
ComlibInstanceHandle commGInstance1;
ComlibInstanceHandle commGInstance2;
ComlibInstanceHandle commGInstance3;
ComlibInstanceHandle commGByrdInstance;
ComlibInstanceHandle commRealInstance;
ComlibInstanceHandle commRealIGXInstance;
ComlibInstanceHandle commRealIGYInstance;
ComlibInstanceHandle commRealIGZInstance;
ComlibInstanceHandle mcastInstance;
ComlibInstanceHandle ssInstance;
ComlibInstanceHandle mssInstance;
ComlibInstanceHandle gssInstance;
ComlibInstanceHandle mcastInstancePP;
CkReduction::reducerType complexVectorAdderType;
#include "ReductionClients.h"
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** The Main of CPAIMD. It calls all the init functions.
 *
 */
//============================================================================

main::main(CkArgMsg *msg) {
//============================================================================
/** Check arguments : Tell people what we are doing */

    done_init=0;
    if (msg->argc < 3) {
      CkAbort("Usage: pgm config_file_charm config_file_piny");
    }//endif

    CkPrintf("\n================================================\n");
    CkPrintf("Starting Cpaimd-Charm-Driver Setup Phase \n");
    CkPrintf("---------------------------------------------------\n");
    CkPrintf("  Cpaimd-Charm-Driver running on %d processors. \n", CkNumPes());
    CkPrintf("  Reading Physics input from %s\n",msg->argv[2]);
    CkPrintf("  Reading Driver  input from %s\n",msg->argv[1]);
    CkPrintf("---------------------------------------------------\n\n");

//============================================================================    
/* Invoke PINY input class */

    CkCallback piny_callback (CkCallback::ignore);
    Interface_ctrl piny_interface (msg->argv[2],piny_callback);

    CPcharmParaInfo *sim  = new CPcharmParaInfo();
    PhysicsParamTransfer::ParaInfoInit(sim);
    int ibinary_opt = sim->ibinary_opt;
    int natm_nl     = sim->natm_nl;

//============================================================================    
/* Invoke Ramkumar input class */

    CkPrintf("\n======================================================\n");
    CkPrintf("Cpaimd-Charm-Driver input started \n");
    CkPrintf("---------------------------------------------------------\n\n");

    Config::readConfig(msg->argv[1],config,sim->nstates,
                       sim->sizeX,sim->sizeY,sim->sizeZ,
                       sim->ntime,ibinary_opt,natm_nl);

    int numSfGrps    = config.numSfGrps;  // local copies are nice
    int doublePack   = config.doublePack;
    size2d sizeYZ    = size2d(sim->sizeY,sim->sizeZ);
    nchareG          = config.nchareG;
    sim->nchareG     = nchareG; 


    nstates          = config.nstates;    // globals on all procs
    sizeX            = sim->sizeX;

    config.print(msg->argv[1]);

    CkPrintf("\n------------------------------------------------\n");
    CkPrintf("Cpaimd-Charm-Driver input completed \n");
    CkPrintf("================================================\n\n");

//============================================================================
// Set user trace events for projections optimizations

     traceRegisterUserEvent("doRealFwFFT", doRealFwFFT_);
     traceRegisterUserEvent("doRealBwFFT", doRealBwFFT_);
     traceRegisterUserEvent("GspaceFwFFT", GspaceFwFFT_);
     traceRegisterUserEvent("GspaceBwFFT", GspaceBwFFT_);
     traceRegisterUserEvent("fwFFTGtoR0", fwFFTGtoR0_);
     traceRegisterUserEvent("fwFFTGtoRnot0", fwFFTGtoRnot0_);
     traceRegisterUserEvent("GradCorrGGA", GradCorrGGA_);
     traceRegisterUserEvent("WhiteByrdFFTX", WhiteByrdFFTX_);
     traceRegisterUserEvent("WhiteByrdFFTY", WhiteByrdFFTY_);
     traceRegisterUserEvent("WhiteByrdFFTZ", WhiteByrdFFTZ_);
     traceRegisterUserEvent("PostByrdfwFFTGtoR", PostByrdfwFFTGtoR_);
     traceRegisterUserEvent("RhoRtoGFFT", RhoRtoGFFT_);
     traceRegisterUserEvent("BwFFTRtoG", BwFFTRtoG_);
     traceRegisterUserEvent("ByrdanddoFwFFTGtoR",ByrdanddoFwFFTGtoR_);

     traceRegisterUserEvent("HartExcVksG",HartExcVksG_);
     traceRegisterUserEvent("divRhoVksGspace",divRhoVksGspace_);

     traceRegisterUserEvent("DoFFTContribute", DoFFTContribute_);
     traceRegisterUserEvent("IntegrateModForces", IntegrateModForces_);
     traceRegisterUserEvent("Scalcmap", Scalcmap_);
     traceRegisterUserEvent("AcceptStructFact", AcceptStructFact_);
     Ortho_UE_step2 = traceRegisterUserEvent("Ortho step 2");
     Ortho_UE_step3 = traceRegisterUserEvent("Ortho step 3");
     Ortho_UE_error = traceRegisterUserEvent("Ortho error");

     /* choose whether ortho should use local callback */
     Ortho_use_local_cb = true;

//============================================================================    
// Compute structure factor grp parameters and static map for chare arrays

    sim->numSfGrps   = numSfGrps;
    int natm_nl_grp_max;
    PhysicsParamTransfer::get_Sfgrp_max(natm_nl,config.numSfGrps, 
                                        &natm_nl_grp_max);
    sim->natm_nl_grp_max = natm_nl_grp_max;

    create_line_decomp_descriptor(sim);

    PhysicsParamTransfer::control_new_mapping_function(sim,doublePack);

    make_rho_runs(sim);

    scProxy  = CProxy_CPcharmParaInfoGrp::ckNew(*sim);
//============================================================================    
// Create the multicast/reduction manager for array sections
// Create the parainfo group from sim
// Initialize chare arrays for real and g-space of states 

    mCastGrpId = CProxy_CkMulticastMgr::ckNew(config.numMulticastMsgs);

    init_state_chares(sizeYZ,natm_nl,natm_nl_grp_max,numSfGrps,doublePack,sim);

//============================================================================    
// Transfer parameters from physics to driver
//    read in atoms : create atoms group 
    
    control_physics_to_driver();

//============================================================================ 
// Initialize the density chare arrays

    init_rho_chares(sizeYZ,sim);

//============================================================================ 
// Create mapping classes for Paircalcular

    mainProxy=thishandle;
  //-------------------------------------------------------------
    int indexSize = nchareG;

    int* indexZ = new int[indexSize];
    for(int i=0, count=0; i<nchareG; i++){
        indexZ[count] = i;
        count++;
    }//endif

//============================================================================ 
// Initialize paircalculators for Psi and Lambda

    init_pair_calculators( nstates,indexSize,indexZ,doublePack,sim);

//============================================================================ 
// initialize Ortho

    init_ortho_chares(nstates, indexSize, indexZ);

//============================================================================ 
// Initialize commlib strategies for later association and delegation

    init_commlib_strategies(sim->nchareRhoG, sizeYZ[1]);

//============================================================================
// clean up

    delete msg;
    delete sim;
    delete [] indexZ;

//============================================================================

    CkPrintf("\n-----------------------------------------------------\n");
    CkPrintf("Cpaimd-Charm-Driver setup phase complete\n");
    CkPrintf("======================================================\n\n");

    CkPrintf("======================================================\n");
    CkPrintf("Launching chare arrays and obtaining data sets  \n");
    CkPrintf("------------------------------------------------------\n\n");

//============================================================================
   }// end Main
//============================================================================


//============================================================================    
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================    
/**
 * Initialize paircalc1 Psi (sym) and paircalc2 Lambda (asym)
 */
//============================================================================    
void init_pair_calculators(int nstates, int indexSize, int *indexZ ,
                           int doublePack, CPcharmParaInfo *sim)
//============================================================================    
   {
//============================================================================    

  PRINT_LINE_STAR;
  PRINTF("Building Psi and Lambda Pair Calculators\n");
  PRINT_LINE_DASH;printf("\n");
  //-------------------------------------------------------------
  // Create mapping classes for Paircalcular

    CProxy_SCalcMap scMap_sym = CProxy_SCalcMap::ckNew(config.nstates,
                   config.nchareG, config.sGrainSize, CmiTrue, sim->nchareG, 
                   sim->lines_per_chareG, sim->pts_per_chareG, config.scalc_per_plane, planes_per_pe, config.numChunks) ;
    
    CProxy_SCalcMap scMap_asym = CProxy_SCalcMap::ckNew(config.nstates,
  	           config.nchareG, config.sGrainSize, CmiFalse, sim->nchareG, 
                   sim->lines_per_chareG, sim->pts_per_chareG, config.scalc_per_plane, planes_per_pe, config.numChunks);
    
    CkGroupID scalc_sym_id  = scMap_sym.ckGetGroupID();
    CkGroupID scalc_asym_id = scMap_asym.ckGetGroupID();


  //-------------------------------------------------------------
  // Register the PCs

    int gsp_ep =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_CkReductionMsg;
    int gsp_ep_tol =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsiV_CkReductionMsg;

   //symmetric AKA Psi
    createPairCalculator(true, nstates, config.sGrainSize, indexSize, indexZ,  CkCallback(CkIndex_Ortho::start_calc(NULL), orthoProxy), &pairCalcID1, gsp_ep, gsp_ep_tol, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_sym_id, doublePack, config.conserveMemory,config.lbpaircalc, config.psipriority, mCastGrpId, config.numChunks );

    CkArrayIndex2D myindex(0, 0);

    gsp_ep = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_CkReductionMsg;
    int myPack=0;

    //asymmetric AKA Lambda AKA Gamma
    createPairCalculator(false, nstates,  config.sGrainSize, indexSize, indexZ,CkCallback(CkIndex_CP_State_GSpacePlane::acceptAllLambda(NULL), myindex, gSpacePlaneProxy.ckGetArrayID()), &pairCalcID2, gsp_ep, 0, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_asym_id, myPack, config.conserveMemory,config.lbpaircalc, config.lambdapriority, mCastGrpId, config.numChunks);

//============================================================================ 
   }//end routine
//============================================================================ 

//============================================================================ 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================ 
/**
 * Initialize Commlib communication strategies  
 */ 
//                                          
//============================================================================
void init_commlib_strategies(int numRhoG, int numReal){
//============================================================================

  PRINT_LINE_STAR;
  PRINTF("Building Commlib strategies\n");
  PRINT_LINE_DASH;printf("\n");
  int i = 0;
  int rhoGhelpers = config.rhoGHelpers;
  int numRhoGHart = rhoGhelpers*numRhoG;

//============================================================================
    if (config.useCommlib) {        
       CkPrintf("Making real_strategy with src numReal %d dest numRhoG %d and numHartG %d\n",
                  numReal,numRhoG,numRhoGHart);
     //--------------------------------------------------------------
     //  For rho(r) to rho(g)
        CkArrayIndexMax *rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
	    CkArrayIndex2D idx2d(i,0);
            rhoGElements[i] = idx2d;
        }//endfor

        CkArrayIndexMax *rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
	    CkArrayIndex2D idx2d(i,0);
            rhoRealElements[i] = idx2d; 
        }//endfor

        CharmStrategy *real_strat = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoRealProxy.ckGetArrayID(), rhoGProxy.ckGetArrayID(),
             numReal, rhoRealElements, numRhoG, rhoGElements);

        commRealInstance= ComlibRegister(real_strat);
     //--------------------------------------------------------------
     //  For drho(r)/dx to igx*rho(g)
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }//endfor
        
        CharmStrategy *real_strat_igx = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoRealProxy.ckGetArrayID(), rhoGProxy.ckGetArrayID(),
             numReal, rhoRealElements, numRhoG,rhoGElements);
        commRealIGXInstance= ComlibRegister(real_strat_igx);
     //--------------------------------------------------------------
     //  For drho(r)/dy to igy*rho(g)
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }//endfor
        
        CharmStrategy *real_strat_igy = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoRealProxy.ckGetArrayID(), rhoGProxy.ckGetArrayID(),
             numReal, rhoRealElements, numRhoG,rhoGElements);
        commRealIGYInstance= ComlibRegister(real_strat_igy);
     //--------------------------------------------------------------
     //  For drho(r)/dz to igz*rho(g)
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }//endfor
        
        CharmStrategy *real_strat_igz = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoRealProxy.ckGetArrayID(), rhoGProxy.ckGetArrayID(),
             numReal, rhoRealElements, numRhoG,rhoGElements);
        commRealIGZInstance= ComlibRegister(real_strat_igz);
     //--------------------------------------------------------------
     // For hartree-Ext(g) to vks(r)
        rhoGElements = new CkArrayIndexMax[numRhoGHart];
        for (i = 0; i < numRhoGHart; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        CharmStrategy *gstrathart = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoGHartExtProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             numRhoGHart, rhoGElements, numReal, rhoRealElements);
        commGHartInstance = ComlibRegister(gstrathart);
     //--------------------------------------------------------------
     // vks(g), igxrho igyrho igzrho and white byrd to g-space
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }//endfor

        CharmStrategy *gstrat0 = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             numRhoG, rhoGElements, numReal, rhoRealElements);
        commGInstance0 = ComlibRegister(gstrat0);
     //--------------------------------------------------------------
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }

        CharmStrategy *gstrat1 = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             numRhoG, rhoGElements, numReal, rhoRealElements);
        commGInstance1 = ComlibRegister(gstrat1);
     //--------------------------------------------------------------
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }

        CharmStrategy *gstrat2 = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             numRhoG, rhoGElements, numReal, rhoRealElements);
        commGInstance2 = ComlibRegister(gstrat2);
     //--------------------------------------------------------------
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }

        CharmStrategy *gstrat3 = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             numRhoG, rhoGElements, numReal, rhoRealElements);
        commGInstance3 = ComlibRegister(gstrat3);
     //--------------------------------------------------------------
     // For white-byrd
        rhoGElements = new CkArrayIndexMax[numRhoG];
        for (i = 0; i < numRhoG; i++) {
            rhoGElements[i] = CkArrayIndex2D(i,0);
        }

        rhoRealElements = new  CkArrayIndexMax[numReal];
        for(i = 0; i < numReal; i++) {
            rhoRealElements[i] = CkArrayIndex2D(i,0);
        }

        CharmStrategy *gstratByrd = new EachToManyMulticastStrategy
            (USE_DIRECT, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             numRhoG, rhoGElements, numReal, rhoRealElements);
        commGByrdInstance = ComlibRegister(gstratByrd);

    }//endif : use commlib

//============================================================================
// Real state space to gspace state and particle plane comm.

    if (config.useCommlibMulticast) {
        DirectMulticastStrategy *dstrat = new DirectMulticastStrategy
	    (realSpacePlaneProxy.ckGetArrayID(),1);
        
        RingMulticastStrategy *rstrat = new RingMulticastStrategy
            (realSpacePlaneProxy.ckGetArrayID(),1);
        
        MultiRingMulticast *mrstrat = new MultiRingMulticast
            (realSpacePlaneProxy.ckGetArrayID(),1);
        
        DirectMulticastStrategy *d1strat = new DirectMulticastStrategy
            (particlePlaneProxy.ckGetArrayID(),1);

        RingMulticastStrategy *r1strat = new RingMulticastStrategy
            (particlePlaneProxy.ckGetArrayID(),1);

        MultiRingMulticast *mr1strat = new MultiRingMulticast
            (particlePlaneProxy.ckGetArrayID(),1);

         //multiring should be good on large runs, but not on BG/L
	if(CkNumNodes()>64){
	      mcastInstance=ComlibRegister(dstrat);
	      mcastInstancePP=ComlibRegister(mr1strat);
        }else{
	      mcastInstance=ComlibRegister(rstrat);
	      mcastInstancePP=ComlibRegister(r1strat);
        }//endif
	
    }// end Sameer's new communication strategies 

//============================================================================
   }//end routine
//============================================================================


//============================================================================
/**
 ** Create stuff for ortho which PC invokes by section reduction
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_ortho_chares(int nstates, int indexSize, int *indexZ){
//============================================================================

  PRINT_LINE_STAR;
  PRINTF("Building Ortho Chares\n");
  PRINT_LINE_DASH;printf("\n");

    int chunks = (nstates + config.sGrainSize - 1) / config.sGrainSize;
    CProxy_OrthoMap orthoMap = CProxy_OrthoMap::ckNew(chunks);
    CkArrayOptions orthoOpts;
    orthoOpts.setMap(orthoMap);
    orthoProxy = CProxy_Ortho::ckNew(orthoOpts);

    CkCallback *orthoReduction = new CkCallback(CkIndex_Ortho::collect_error(NULL), orthoProxy(0, 0));
    orthoProxy.ckSetReductionClient(orthoReduction);
    
    // extra triangle ortho elements are really a waste of our time
    // and resources, but we don't have a triangular solver for
    // inv_square, so we'll just make do.

    // They need to exist solely so that the inv_sq method can work.
    // So we need to copy their mirror elements data into them.
    // then when complete they need to know not to call finishpaircalc.
    // Because their redundant data has nowhere to go.

    /* create matrix multiplication objects */
    CLA_Matrix_interface matA1, matB1, matC1;
    CLA_Matrix_interface matA2, matB2, matC2;
    CLA_Matrix_interface matA3, matB3, matC3;

    CkCallback ortho_ready_cb = CkCallback(CkIndex_Ortho::all_ready(),
     orthoProxy(0, 0));
    make_multiplier(&matA1, &matB1, &matC1, orthoProxy, orthoProxy, orthoProxy,
     nstates, nstates, nstates, config.sGrainSize, config.sGrainSize,
     config.sGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
     mCastGrpId, MM_ALG_2D);
    make_multiplier(&matA2, &matB2, &matC2, orthoProxy, orthoProxy, orthoProxy,
     nstates, nstates, nstates, config.sGrainSize, config.sGrainSize,
     config.sGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
     mCastGrpId, MM_ALG_2D);
    make_multiplier(&matA3, &matB3, &matC3, orthoProxy, orthoProxy, orthoProxy,
     nstates, nstates, nstates, config.sGrainSize, config.sGrainSize,
     config.sGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
     mCastGrpId, MM_ALG_2D);


    for (int s1 = 0; s1 < nstates; s1 += config.sGrainSize)
      for (int s2 = 0; s2 < nstates; s2 += config.sGrainSize) {
	int indX = s1 / config.sGrainSize;
	int indY = s2 / config.sGrainSize;
	orthoProxy(indX, indY).insert(config.sGrainSize, config.sGrainSize,
         matA1, matB1, matC1, matA2, matB2, matC2, matA3, matB3, matC3);
      }
    orthoProxy.doneInserting();
    orthoProxy.makeSections(indexSize, indexZ);

//============================================================================
   }//end routine 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void main::doneInit(CkReductionMsg *msg){
//============================================================================
    delete msg;

    if(done_init<3){
      CkPrintf("Completed chare instantiation phase %d\n",done_init+1);
    }else{
      CkPrintf("Completed chare data acquisition phase %d\n",done_init+1);
      CkPrintf("\n-----------------------------------------------------\n");
      CkPrintf("Chare array launch and initialization complete       \n");
      CkPrintf("======================================================\n\n");
    }//endif

    if (done_init == 2){
	// kick off file reading in gspace
	CkPrintf("Initiating import of states\n");
	for(int s=0;s<nstates;s++){
	    gSpacePlaneProxy(s,0).readFile();
	}//endfor
    }//endif
    if (done_init >= 3) {
      if (done_init == 3){ 
 	  CkPrintf("\n======================================================\n");
          if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==1){
            CkPrintf("Running Open Atom CP Minimization: \n");
	  }else{
            CkPrintf("Running Open Atom CP Dynamics: \n");
	  }//endif
  	  CkPrintf("======================================================\n\n");
  	  CkPrintf("\n======================================================\n");
	  gSpacePlaneProxy.run();
      }//endif
    }
    done_init++;
}
//============================================================================


//============================================================================
//Create the array elements for the GSpace, Particle and Real Space planes
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_state_chares(size2d sizeYZ, int natm_nl,int natm_nl_grp_max,int numSfGrps,
                 int doublePack, CPcharmParaInfo *sim)

//============================================================================
   { //begin routine 
//============================================================================
    /*
     * Set up the map variables used to control the location of the
     * 2d-array elements over processors   
     */

  PRINT_LINE_STAR;
  PRINTF("Building G-space and R-space Chares state %d sizeYZ %d %d\n",nstates,sizeYZ[0],sizeYZ[1]);
  PRINT_LINE_DASH;printf("\n");


    CProxy_RSMap rsMap= CProxy_RSMap::ckNew(config.nstates, sim->sizeY, config.Rstates_per_pe);

    CkArrayOptions realSpaceOpts;
    realSpaceOpts.setMap(rsMap);
    size2d sizeRealPlane(sizeYZ[1], sizeX);

    realSpacePlaneProxy = CProxy_CP_State_RealSpacePlane::ckNew(sizeRealPlane,1,1,
                                                                realSpaceOpts);

								
    fftCacheProxy = CProxy_FFTcache::ckNew(sizeRealPlane,1);
    sfCacheProxy = CProxy_StructFactCache::ckNew(numSfGrps,natm_nl,natm_nl_grp_max);
    sfCompProxy = CProxy_StructureFactor::ckNew();
    
    CProxy_GSMap gsMap = CProxy_GSMap::ckNew(sim->nchareG, sim->lines_per_chareG, sim->pts_per_chareG, config.nstates, config.Gstates_per_pe);

    CkArrayOptions gSpaceOpts;
    gSpaceOpts.setMap(gsMap);

    gSpacePlaneProxy = CProxy_CP_State_GSpacePlane::ckNew(sizeX, sizeYZ,1, 
                                                          1,config.sGrainSize, config.numChunks, gSpaceOpts);

    // We bind the particlePlane array to the gSpacePlane array migrate together
    CkArrayOptions particleOpts;
    particleOpts.setMap(gsMap); // the maps for both the arrays are the same
    particleOpts.bindTo(gSpacePlaneProxy);
    particlePlaneProxy = CProxy_CP_State_ParticlePlane::ckNew(nchareG, sizeYZ[0], sizeYZ[1],   
		      1,numSfGrps,natm_nl,natm_nl_grp_max,particleOpts);

    /*
     * Insert the planes in the particle plane array, gSpacePlane array
     */

    int s,x;
    CkPrintf("Making (nstates=%d)x(nchareG=%d) gspace objects\n\n",nstates,nchareG);
    for (s = 0; s < nstates; s++){
      for (x = 0; x <nchareG; x++){
             gSpacePlaneProxy(s, x).insert(sizeX, sizeYZ, 1, 
                                    1,config.sGrainSize, config.numChunks);
             particlePlaneProxy(s, x).insert(sizeX, sizeYZ[0], sizeYZ[1],   
				      1,numSfGrps,natm_nl,natm_nl_grp_max);
      }
    }

    gSpacePlaneProxy.doneInserting();
    particlePlaneProxy.doneInserting();
    particlePlaneProxy.setReductionClient(doneCreatingPP, (void *) NULL);
    
    /*
     * Insert the planes in the real space plane array
     */
    int y;
    for (s = 0;  s < nstates; s++)
        {for (y = 0; y < sizeYZ[0]; y += 1)
            {realSpacePlaneProxy(s, y).insert(sizeRealPlane,1,1);}}

    realSpacePlaneProxy.doneInserting();
    // 
    //    gSpacePlaneProxy.setReductionClient(doneInit, (void *) NULL);
    //    realSpacePlaneProxy.setReductionClient(doneInit, (void *) NULL);

    int *nsend   = new int[nchareG];
    int **listpe = new int * [nchareG];
    int numproc  = CkNumPes();
    int *gspace_proc = new int [numproc];
    
   
    for(int i =0;i<numproc;i++){gspace_proc[i]=0;}
    for(int j=0;j<nchareG;j++){   
      listpe[j]= new int[nstates];
      nsend[j]=0;
      for(int i=0;i<nstates;i++){
	listpe[j][i]=cheesyhackgsprocNum(sim, i,j);
	gspace_proc[listpe[j][i]]+=1;
      }//endfor
      lst_sort_clean(nstates, &nsend[j], listpe[j]);
    }//endfor
    
#ifdef USE_TOPOMAP
    if(maptable!=NULL) //clean up cheesyhack
	delete maptable;
#endif
    
    FILE *fp = fopen("gspplane_proc_distrib.out","w");
    for(int i=0;i<numproc;i++){
      fprintf(fp,"%d %d\n",i,gspace_proc[i]);
    }//endfor
    fclose(fp);
    delete [] gspace_proc;


    int minsend=nstates;
    int maxsend=0;
    double avgsend=0.0;
    int chareG_use=0;
    CkPrintf("============================\n");
    CkPrintf("Structure factor chareG dests\n");
    CkPrintf("Number of g-space chares : %d\n",nchareG);    
    CkPrintf("---------------------------\n");
    for(int lsi=0;lsi<nchareG;lsi++){
        chareG_use++;
	CkPrintf("chareG [%d] nsend %d\n",lsi,nsend[lsi]);
	if(nsend[lsi]>maxsend)
	    maxsend=nsend[lsi];
	if(nsend[lsi]<minsend)
	    minsend=nsend[lsi];
	avgsend+=(double) nsend[lsi];
#ifdef SF_SEND_LIST_OUT
        for(int lsj=0;lsj<nsend[lsi];lsj++){
	  CkPrintf("[%d %d] pe %d\n",lsi,lsj,listpe[lsi][lsj]);
	}//endfor
#endif
    }//endfor
    CkPrintf("---------------------------\n");
    CkPrintf("SFSends min %d max %d avg %g\n",minsend,maxsend,avgsend/(double)chareG_use);
    CkPrintf("============================\n");
    // Insert the objects into the StructureFactor array

    int dupmax=maxsend;  // there is no point in ever having more than that

//    if (dupmax>=natm_nl) // stability issues if over 
//      dupmax=natm_nl-1;  
    
    if (config.numSfDups<dupmax)
	dupmax=config.numSfDups;
    config.numSfDups=dupmax;
    CkPrintf("real numSfdups is %d based on maxsend of %d\n",config.numSfDups, maxsend);
    for (int dup=0; dup<dupmax; dup++)
      for (x = 0; x < nchareG; x += 1)
      {
	  int num_dup, istart, iend;
	  get_grp_params( nsend[x],  config.numSfDups,  dup, x ,&num_dup,  &istart, &iend);
	  int pe_ind=istart;
	  if(x%2==0)
	      pe_ind=iend;
	  for (int AtmGrp=0; AtmGrp<numSfGrps; AtmGrp++)
	  {
	      sfCompProxy(AtmGrp, x, dup).insert(numSfGrps, config.numSfDups, natm_nl_grp_max,  num_dup, &(listpe[x][istart]),atmGrpMap(istart, num_dup, nsend[x], listpe[x], AtmGrp, dup,x));	      
	      pe_ind++;
	      if(pe_ind>nsend[x])
		  pe_ind=0;
	  }
      }
    sfCompProxy.doneInserting();

    for(int j=0;j<nchareG;j++)
      delete [] listpe[j];
    delete [] listpe;
    delete [] nsend;
    if(config.useCommlib) {
        // Set some com strategy of Sameer
        ssInstance = CkGetComlibInstance();
        StreamingStrategy *strat = new StreamingStrategy(0.2,5);
        //strat->enableShortArrayMessagePacking();
        ssInstance.setStrategy(strat);
        
        // Set some com strategy of Sameer
        StreamingStrategy *mstrat = new StreamingStrategy(0.2,5);
        //mstrat->enableShortArrayMessagePacking();
        mssInstance= ComlibRegister(mstrat);    
        StreamingStrategy *gmstrat = new StreamingStrategy(0.2,5);
        //mstrat->enableShortArrayMessagePacking();
        gssInstance= ComlibRegister(gmstrat);    

    }

  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed g-space chare build\n");
  PRINT_LINE_STAR;printf("\n");

//----------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
// Creating arrays CP_Rho_GSpacePlane, CP_Rho_GSpacePlaneHelper 
// and CP_Rho_RealSpacePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_rho_chares(size2d sizeYZ, CPcharmParaInfo *sim)

//============================================================================
{//begin routine
//============================================================================
/*
 * create the array for real-space densities (two-dimensional chare array)
 */    
//============================================================================
//  Chare array sizes and offsets 

    int nchareRhoR     = sizeYZ[1];
    int nchareRhoG     = sim->nchareRhoG;
    int rhoGHelpers    = config.rhoGHelpers;
    int nchareRhoGHart = rhoGHelpers*nchareRhoG;
    
    int rhoRstride     = CkNumPes()/(nchareRhoR);
    int rhoGstride     = CkNumPes()/(nchareRhoG);
    int rhoGHartstride = CkNumPes()/(nchareRhoGHart);
    if(rhoRstride<1)    {rhoRstride=1;}
    if(rhoGstride<1)    {rhoGstride=1;}
    if(rhoGHartstride<1){rhoGHartstride=1;}
    if(rhoRstride==rhoGstride && rhoGstride!=1){rhoRstride--;}

    int rhoGHartoff   =0;
    int offsetFromZero=1;  // proc0 has reduction root issues. Avoid it if possible
    if(rhoGstride==1)   {offsetFromZero=0;}
    if(offsetFromZero>0){rhoGHartoff=2;}

//============================================================================
// Maps and options

    CProxy_RhoRSMap rhorsMap = CProxy_RhoRSMap::ckNew(rhoRstride,offsetFromZero, nchareRhoR);
    CkArrayOptions rhorsOpts;
    rhorsOpts.setMap(rhorsMap);

    CProxy_RhoGSMap rhogsMap = CProxy_RhoGSMap::ckNew(rhoGstride,offsetFromZero, rhoRstride, offsetFromZero, nchareRhoG);
    CkArrayOptions rhogsOpts;
    rhogsOpts.setMap(rhogsMap);

    CProxy_RhoGSMap rhogHartMap = CProxy_RhoGHartMap::ckNew(rhoGHartstride, rhoGHartoff, rhoRstride, offsetFromZero, nchareRhoGHart);
    CkArrayOptions rhoghartOpts;
    rhoghartOpts.setMap(rhogHartMap);

//============================================================================
// Instantiate the chares

    bool dummy = true;

//---------------------------------------------------------------------------
// rhoreal
    rhoRealProxy = 
      CProxy_CP_Rho_RealSpacePlane::ckNew(sizeX, sizeYZ, 1, 1, dummy, rhorsOpts);
    for (int i = 0; i < nchareRhoR; i++) {
	rhoRealProxy(i,0).insert(sizeX, sizeYZ, 1,1,dummy);
    }//endfor
    rhoRealProxy.doneInserting();
    rhoRealProxy.setReductionClient(printEnergyEexc, 0);
//---------------------------------------------------------------------------
// rhog
    rhoGProxy = CProxy_CP_Rho_GSpacePlane::ckNew(sizeX, sizeYZ, 1, 
						 1, dummy, 
						 rhogsOpts);
    for (int i = 0; i < nchareRhoG; i++){
	rhoGProxy(i,0).insert(sizeX, sizeYZ,1,1,dummy );
    }//endfor
    rhoGProxy.doneInserting();
//---------------------------------------------------------------------------
// rhoghart
    rhoGHartExtProxy = CProxy_CP_Rho_GHartExt::ckNew(sizeYZ, rhoghartOpts);
    for (int i = 0; i < nchareRhoGHart; i++){
      rhoGHartExtProxy(i,0).insert(sizeYZ);
    }//endfor
    rhoGHartExtProxy.setReductionClient(printEnergyHart, NULL);
    rhoGHartExtProxy.doneInserting();

//===========================================================================
  }//end routine
//============================================================================



//============================================================================
// Get the atoms and the parainfo
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void control_physics_to_driver(){

//============================================================================
// make a group : create a proxy for the atom class and also a reduction client

    PhysicsAtomPosInit *PhysicsAtom  = new PhysicsAtomPosInit();
    int natm          = PhysicsAtom->natm_tot;
    int natm_nl       = PhysicsAtom->natm_nl;
    int len_nhc       = PhysicsAtom->len_nhc;
    int iextended_on  = PhysicsAtom->iextended_on;
    int cp_min_opt    = PhysicsAtom->cp_min_opt;
    int cp_wave_opt   = PhysicsAtom->cp_wave_opt;
    int isokin_opt    = PhysicsAtom->isokin_opt;
    double kT         = PhysicsAtom->kT;
  
    Atom *atoms       = new Atom[natm];
    AtomNHC *atomsNHC = new AtomNHC[natm];

    PhysicsAtom->DriverAtomInit(natm,atoms,atomsNHC);
    atomsGrpProxy = CProxy_AtomsGrp::ckNew(natm,natm_nl,len_nhc,iextended_on,
                                           cp_min_opt,cp_wave_opt,isokin_opt,
                                           kT,atoms,atomsNHC);

    delete [] atoms;
    delete [] atomsNHC;
    delete PhysicsAtom;

//=====================================================================
// Make a group for the energies

    egroupProxy = CProxy_EnergyGroup::ckNew(); 

//----------------------------------------------------------------------------
    }
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
static CkReductionMsg *complexSum(int nMsg, CkReductionMsg **msgs)
{
    // get the size of the messages from the first message
    int matrixSize = msgs[0]->getSize() / sizeof(complex);
    complex *matrix = new complex[matrixSize];
    int i, m;
    for (m = 0; m < nMsg; m++) {
        // sanity check
        CkAssert(msgs[m]->getSize() == matrixSize * sizeof(complex));
        complex *data = (complex *) msgs[m]->getData();
        for (i = 0; i < matrixSize; i++)
			matrix[i] += data[i];
    }
    CkReductionMsg *ret = CkReductionMsg::buildNew(matrixSize * 
                                                   sizeof(complex), matrix);
    delete [] matrix;
    return ret;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void get_grp_params(int natm_nl, int numSfDups, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret)
{

   int n     = (natm_nl/numSfDups);
   int m     = (natm_nl % numSfDups);

   int istrt = n*indexSfGrp;
   if(indexSfGrp>=m){istrt += m;}
   if(indexSfGrp<m) {istrt += indexSfGrp;}
   if(indexSfGrp<m) {n++;}
   int iend  = n+istrt;
   if(numSfDups>natm_nl)
     if(m>=indexSfGrp)
       {
	 n=0;
	 istrt=natm_nl+1;
	 iend=natm_nl;
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
	 CkPrintf("Redundant DupSF for chare-G %d\n",planeIndex);
         CkPrintf("At present this hangs, so out you go\n");
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
       }
  
   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;
}
//============================================================================



//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex)
{

  int numSfDups=config.numSfDups;
  if(listsize <= numSfDups)
    { //dup list not unique
      return listpe[(dup % listsize)];
    }
  if(nsend <= config.numSfGrps)
    {  //atom group list not unique
      return listpe[(istart + (AtmGrp % nsend))];
    }
  // nsend must be > numSfGrps and its list is unique
  return listpe[(istart + ((AtmGrp + config.numSfGrps * planeIndex)%nsend))];

}
//============================================================================


void makemap()
{
#ifdef USE_TOPOMAP
	FindProcessor fp;
	int x = CkNumPes();
	int y = 1;
	int z = 1;
	int vn = 0;
        int c = 0;
	
#if CMK_VERSION_BLUEGENE
	bgltm = BGLTorusManager::getObject();
	x = bgltm->getXSize();
	y = bgltm->getYSize();
	z = bgltm->getZSize();
#endif
	fp.count=1;
	fp.nopX=x;
	fp.nopY=y;
	fp.nopZ=z;

	int assign[3]={0, 0, 0};
	int w = 0;
	for(int i=0;i<3;i++)
		fp.start[i]=fp.next[i]=0;
	int gsobjs_per_pe;
	
	if((config.nstates*config.nchareG) % CkNumPes() == 0)
	    gsobjs_per_pe = (config.nstates*config.nchareG)/CkNumPes();
	else
	    gsobjs_per_pe = (config.nstates*config.nchareG)/CkNumPes()+1;
	int l=config.Gstates_per_pe;
	int m, pl, pm, rem;
        
        pl = config.nstates / l;
        pm = CkNumPes() / pl;
       
	if(pm==0)
          CkAbort("Choose a larger Rstates_per_pe\n");\
 
        m = config.nchareG / pm;
        rem = config.nchareG % pm;

        planes_per_pe=m;
        
        //CkPrintf("nstates %d nchareG %d Pes %d, gsobjs_per_pe %d\n", config.nstates, config.nchareG, CkNumPes(), gsobjs_per_pe);	
	//CkPrintf("l %d, m %d pl %d pm %d rem %d\n", l, m, pl, pm, rem);
	
        for(int ychunk=0; ychunk<config.nchareG; ychunk=ychunk+m)
        {
                if(ychunk==(pm-rem)*m)
                  m=m+1;
                for(int xchunk=0; xchunk<config.nstates; xchunk=xchunk+l)
		{
			if(xchunk==0 && ychunk==0) {}
			else
			{
				for(int i=0;i<3;i++)
					fp.start[i]=fp.next[i];
				if(fp.start[2]>x/2)
					assign[2]=fp.start[2]-x;
				else
					assign[2]=fp.start[2];
				if(fp.start[1]>y/2)
					assign[1]=fp.start[1]-y;
				else
					assign[1]=fp.start[1];
				if(fp.start[0]>z/2)
					assign[0]=fp.start[0]-z;
				else
					assign[0]=fp.start[0];
				if(vn==0)
					fp.findNextInTorus(assign);
                                else
                                {
                                	fp.findNextInTorusV(w, assign);
                                        w = fp.w;
                                }	
			}
                        c=0;
			for(int state=xchunk; state<xchunk+l && state<config.nstates; state++)
			{
				for(int plane=ychunk; plane<ychunk+m && plane<config.nchareG; plane++)
				{
					if(xchunk==0 && ychunk==0)
					{
                                                c++;
						maptable->put(intdual(state, plane))=0;
						//CkPrintf("%d %d on 0\n", state, plane);
					}
					else
					{
                                                c++;
						if(vn==0)
						  maptable->put(intdual(state, plane))=fp.next[0]*x*y+fp.next[1]*x+fp.next[2];
                                                else
                                                  maptable->put(intdual(state, plane))=(fp.next[0]*x*y+fp.next[1]*x+fp.next[2])*2+fp.w;
						//CkPrintf("%d %d on %d\n", state, plane, fp.next[0]*x*y+fp.next[1]*x+fp.next[2]);
						
					}
				}
			}
                }
        }
#ifdef MAP_DEBUG
	CkPrintf("Local gsmap created on processor %d\n", CkMyPe());
#endif
#endif
}

#ifdef USE_TOPOMAP
int cheesyhackgsprocNum(CPcharmParaInfo *sim,int state, int plane)
{
	if(maptable==NULL)
	{
		maptable= new CkHashtableT<intdual, int> (config.nstates*config.nchareG);
		makemap();
	}
	return maptable->get(intdual(state, plane));
}

#else

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int cheesyhackgsprocNum(CPcharmParaInfo *sim,int state, int plane) {
//============================================================================

  CkArrayIndex2D idx2d(state,plane);

  int pe        = 0;
  int numChareG = 0;
  
  if(config.doublePack) {
     numChareG = nchareG;
  }else{
     CkAbort("not doublepack is broken\n");
  }//endif
  
  double state_load=0.0;
  if(state_load <= 0.0 ){
    for(int x = 0; x  < numChareG; x++) {
      double curload = 0.0;
      double sload = 0.0;
      hackGSpacePlaneLoad(sim, x, &curload, &sload);
      state_load += curload;
    }//endfor
  }//endif

  int pes_per_state = config.GpesPerState;
  int np            = CkNumPes()/pes_per_state;
  
  if(np < 1){np = 1;}
  
  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0){partition_nstates ++;}
  
  int start_pe    = (idx2d.index[0]/partition_nstates) * pes_per_state;
  int start_state = (idx2d.index[0]/partition_nstates) * partition_nstates;
  
  double cum_load     = 0.0;
  double average_load = state_load * config.nstates / CkNumPes();

  for(int x = 0; x < numChareG; x++) {
    for(int s = 0; s < partition_nstates; s++) {
      double curload = 0.0;
      double sload = 0.0;
      
      hackGSpacePlaneLoad(sim,x, &curload, &sload);

      cum_load += curload;

      if((idx2d.index[0] == s + start_state) && idx2d.index[1] == x) {
      
	double dpe = 0.0;
	dpe = cum_load / average_load;

	pe = (int)dpe;
	pe = pe % pes_per_state;
	pe += start_pe;

	return pe % CkNumPes();
      }//endif
    }//endfor : s
  }//endfor : x

  //  CkPrintf("Warning pe not found for index [%d, %d]\n",idx2d.index[0],idx2d.index[1]);
  return (idx2d.index[0]*1037+idx2d.index[1])%CkNumPes();

//============================================================================
  }//end routine
//============================================================================

#endif

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int idx, double *line_load, 
                         double *pt_load){
//============================================================================

  int nchareG             = sim->nchareG;
  double *lines_per_chareG = sim->lines_per_chareG;
  double *pts_per_chareG   = sim->pts_per_chareG;

  if( (idx < 0) || (idx >= nchareG) ){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nchareG);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  line_load[0] = lines_per_chareG[idx];
  pt_load[0]   =   pts_per_chareG[idx];

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void mapOutput()
//============================================================================
    {//begin routine
//============================================================================

	int i;
        i = 1;

//============================================================================
  }//end routine
//============================================================================


//============================================================================
#include "cpaimd.def.h"
//============================================================================
