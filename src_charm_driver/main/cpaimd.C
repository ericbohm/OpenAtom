/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

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
 *  OpenAtom <A HREF="http://ccharm.cs.uiuc.edu/OpenAtom">Webpage</A>.
 *
 */

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================ 
#include <math.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <string>
#include "charm++.h"
#include "ckarray.h"
#include "util.h"
//============================================================================
#include "cpaimd.h"
#include "InstanceController.h"
#include "ckPairCalculator.h"
#include "groups.h"
#include "orthoHelper.h"
#include "ortho.h"
#include "fftCacheSlab.h"
#include "eesCache.h"
#include "StructFactorCache.h"
#include "StructureFactor.h"
#include "cp_state_ctrl/CP_State_ParticlePlane.h"
#include "CP_State_Plane.h"
#include "MeshStreamingStrategy.h"
#include "MultiRingMulticast.h"
#include "PeList.h"
#include "MapFile.h"
#include "TopoManager.h"
#include "TimeKeeper.h"

//============================================================================
#include "../include/debug_flags.h"
#include "../../src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "../../src_piny_physics_v1.0/include/charm_defs/Interface_ctrl.decl.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h"
//============================================================================

int TimeKeeperID=0;
vector <string> TimeKeeperNames;
UberCollection thisInstance;
//============================================================================
/** 
 * \defgroup torus_vars Defining the size of the torus, handy when debugging 
 * torus map logic on non torus architectures.
 */
int numPes;
bool fakeTorus;

//============================================================================
//extern "C" {int get_memory_allocated_user_total();}
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

bool firstInstance = true;  

// INT_MAPs are the ones actually used so
// these are being changed to CkVec's for beads
CkVec <MapType2> GSImaptable;
CkVec <MapType2> RSImaptable;
CkVec <MapType2> RPPImaptable;
CkVec <MapType2> RhoGSImaptable;
CkVec <MapType2> RhoRSImaptable;
CkVec <MapType2> RhoGHartImaptable;
CkVec <MapType3> RhoRHartImaptable;
//CkVec <MapType1> LSPRhoGSImaptable;
CkVec <MapType2> LSPRhoRSImaptable;
CkVec <MapType2> OrthoImaptable;
CkVec <MapType2> OrthoHelperImaptable;
CkVec <MapType4> AsymScalcImaptable;
CkVec <MapType4> SymScalcImaptable;

CkVec <CkGroupID> mCastGrpIds;
#ifndef USE_INT_MAP
CkHashtableT<intdual, int> GSmaptable(10000,0.25);
CkHashtableT<intdual, int> RSmaptable(10000,0.25);
CkHashtableT<intdual, int> RPPmaptable(10000,0.25);
CkHashtableT<intdual, int> AsymScalcmaptable(10000,0.25);
CkHashtableT<intdual, int> SymScalcmaptable(10000,0.25);
#else
CkHashtableT<intdual, int> GSmaptable;
CkHashtableT<intdual, int> RSmaptable;
CkHashtableT<intdual, int> RPPmaptable;
CkHashtableT<intdual, int> AsymScalcmaptable;
CkHashtableT<intdual, int> SymScalcmaptable;
#endif
CkHashtableT<intdual, int> RhoGSmaptable;
CkHashtableT<intdual, int> RhoRSmaptable;
CkHashtableT<intdual, int> RhoGHartmaptable;
CkHashtableT<inttriple, int> RhoRHartmaptable;
CkHashtableT<intdual, int> Orthomaptable;
CkHashtableT<intdual, int> OrthoHelpermaptable;

CkVec <PairCalcID> UpairCalcID1;
CkVec <PairCalcID> UpairCalcID2;


CProxy_main                       mainProxy;
CProxy_CPcharmParaInfoGrp         scProxy;
Config                            config;
CProxy_TimeKeeper                 TimeKeeperProxy;
CProxy_OrthoMap                   orthoMap;
CProxy_InstanceController         instControllerProxy;

//============================================================================
/** Uber proxies for all the things which change per step 
 *  Indexed by PathIntegral Bead.  Each Bead has its own set of
 *  proxies.  Charm driver startup will construct a different set of
 *  arrays for each bead.  
 *
 *  There is a small flexibility vs performance tradeoff.  If we
 *  assume beads are never co-mapped, then the old readonly can be
 *  overwritten locally on each processor
 *  (e.g. gSpacePlaneProxy=UgSpacePlaneProxy[mybead];)
 *  Otherwise we force a slight indirection penalty to lookup 
 *  U*Proxy[mybead] for every send.  In practice this is probably
 *  noise compared to the real expense of sending a message, but it
 *  does seem a little silly for the default case where there is only
 *  one bead.
 */

CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;
CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
CkVec <CProxy_CP_State_RealSpacePlane>    UrealSpacePlaneProxy;
CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
CkVec <CProxy_CP_Rho_GSpacePlane>         UrhoGProxy;
CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
CkVec <CProxy_CP_Rho_GHartExt>            UrhoGHartExtProxy;
CkVec <CProxy_Ortho>                      UorthoProxy;
CkVec <CProxy_AtomsGrp>                   UatomsGrpProxy;
CkVec <CProxy_EnergyGroup>                UegroupProxy;
CkVec <CProxy_FFTcache>                   UfftCacheProxy;
CkVec <CProxy_StructFactCache>            UsfCacheProxy;
CkVec <CProxy_StructureFactor>            UsfCompProxy;
CkVec <CProxy_eesCache>                   UeesCacheProxy;
CkVec <CProxy_OrthoHelper>                UorthoHelperProxy;
CkVec <CProxy_CP_LargeSP_RhoGSpacePlane>      UlsRhoGProxy;
CkVec <CProxy_CP_LargeSP_RhoRealSpacePlane>      UlsRhoRealProxy;

CkVec <UberCollection>			  UberAlles;

//============================================================================


//============================================================================
// readonly globals

double Timer;
int nstates; 
int sizeX;
int nchareG;
int Ortho_UE_step2;
int Ortho_UE_step3;
int Ortho_UE_error;
bool Ortho_use_local_cb;
int done_init=0;
int planes_per_pe;



CkVec < CkVec <int> > UpeUsedBySF;
CkVec < CkVec <int> > UpeUsedByNLZ;
CkVec < CkVec <int> > UplaneUsedByNLZ;
PeList *availGlobR=NULL;
PeList *availGlobG=NULL;
PeList *excludePes=NULL;
int boxSize;
TopoManager *topoMgr=NULL;

//============================================================================


//============================================================================
// For using the multicast library :  Set some reduction clients

CkGroupID            mCastGrpId; 
CkGroupID            orthomCastGrpId; 
CkGroupID            orthoRedGrpId; 
ComlibInstanceHandle orthoInstance;
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

ComlibInstanceHandle gAsymInstance;
ComlibInstanceHandle gSymInstance;

ComlibInstanceHandle mcastInstance;
ComlibInstanceHandle mcastInstancePP;
ComlibInstanceHandle mcastInstanceRPP;
ComlibInstanceHandle mcastInstancemRPP;

ComlibInstanceHandle mssInstance;
ComlibInstanceHandle gssInstance;

ComlibInstanceHandle gssPInstance;
ComlibInstanceHandle mssPInstance;

ComlibInstanceHandle commRHartGHartIns;
ComlibInstanceHandle commGHartRHartIns0;
ComlibInstanceHandle commGHartRHartIns1;

CkReduction::reducerType complexVectorAdderType;
//============================================================================



//============================================================================
/** The Main of CPAIMD. It calls all the init functions.
 *
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
main::main(CkArgMsg *msg) {
  topoMgr = NULL;
//============================================================================
/** Check arguments : Tell people what we are doing */

    done_init=0;
    if (msg->argc < 3) {
      CkAbort("Usage: cpaimd.x cpaimd_config pinysystem.input");
    }//endif
    CkPrintf("Executing OpenAtom: BINARY - %s\n", msg->argv[0]);
    if(msg->argc >3 && msg->argv[3][0] == 't')
      {

	//get system name
	char tidyphysfname[1024];
	bzero(tidyphysfname,1024);
	char simfname[1024];
	bzero(simfname,1024);
	char *lastslash=strrchr(msg->argv[1],'/');
	if(lastslash==NULL)
	  {
	    lastslash=msg->argv[2];
	    strncat(tidyphysfname,"./tidy ",8);
	  }
	else
	  {
	    strncat(tidyphysfname,msg->argv[2],lastslash-msg->argv[2]);
	    strncpy(simfname,tidyphysfname,1024);
	    strncat(tidyphysfname,"/tidy ",8);
	    strncat(simfname,"/",2);
	    lastslash++;
	  }
	strncat(tidyphysfname,lastslash,strchr(msg->argv[2],'.') - lastslash);
	strncat(simfname,lastslash,strchr(msg->argv[2],'.') - lastslash);
	CkPrintf("  Tidy mode, running %s\n",tidyphysfname);
	unlink(simfname);
	system(tidyphysfname);
	strncpy(tidyphysfname,simfname,1024);
	strncat(tidyphysfname,".coords_out",20);
	CkPrintf("  Tidy mode, unlinking %s\n",tidyphysfname);
	unlink(tidyphysfname);
	strncpy(tidyphysfname,simfname,1024);
	strncat(tidyphysfname,".coords.out",20);
	CkPrintf("  Tidy mode, unlinking %s\n",tidyphysfname);
	unlink(tidyphysfname);
	strncpy(tidyphysfname,msg->argv[2],1024);
	strncat(tidyphysfname,".out",20);
	CkPrintf("  Tidy mode, unlinking %s\n",tidyphysfname);
	unlink(tidyphysfname);
	strncpy(tidyphysfname,simfname,1024);
	strncat(tidyphysfname,".confp",20);
	CkPrintf("  Tidy mode, unlinking %s\n",tidyphysfname);
	unlink(tidyphysfname);
	strncpy(simfname,msg->argv[1],1024);
	strncat(simfname,".out",4);
	CkPrintf("  Tidy mode, unlinking %s\n",simfname);
	unlink(simfname);
	sleep(1);
      }
    CkPrintf("\n");
    PRINT_LINE_STAR;
    CkPrintf("Starting Cpaimd-Charm-Driver Setup Phase\n");
    PRINT_LINE_DASH;
    CkPrintf("  Cpaimd-Charm-Driver running on %d processors. \n", CkNumPes());
    CkPrintf("  Reading Physics input from %s\n",msg->argv[2]);
    CkPrintf("  Reading Driver  input from %s\n",msg->argv[1]);

    PRINT_LINE_DASH; CkPrintf("\n");

//============================================================================    
// check the debug flags for consistency

#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
#ifdef _CP_DEBUG_SFNL_OFF_
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't test non-local by itself and turn it off\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
#endif
#endif

#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
#ifndef _CP_DEBUG_VKS_OFF_
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't test non-local by itself with vks on\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
#endif
#endif

#ifdef _CP_DEBUG_VKS_ONLY_
#ifdef _CP_DEBUG_VKS_OFF_
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't test vks by itself and turn it off\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
#endif
#endif

#ifdef _CP_DEBUG_VKS_ONLY_
#ifndef _CP_DEBUG_SFNL_OFF_
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't test vks by itself with SNFL on\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
#endif
#endif

#ifdef _CP_DEBUG_NON_LOCAL_VKS_ONLY_
#ifdef _CP_DEBUG_SFNL_OFF_
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't test vks and SNFL only without sfnl \n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
#endif
#ifdef _CP_DEBUG_VKS_OFF_
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't test vks and SNFL only without vks \n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
#endif
#endif

//============================================================================    
/* Invoke PINY input class */

    CkCallback piny_callback (CkCallback::ignore);
    Interface_ctrl piny_interface (msg->argv[2],piny_callback);

    CPcharmParaInfo *sim  = new CPcharmParaInfo();
    PhysicsParamTransfer::ParaInfoInit(sim);

    int ibinary_opt    = sim->ibinary_opt;
    int natm_nl        = sim->natm_nl;
    int ees_eext_opt   = sim->ees_eext_on;
    int nchareRhoRHart = sim->ngrid_eext_c;
    int fftopt         = sim->fftopt;
    int natm_typ       = sim->natm_typ;

//============================================================================    
/* Invoke parallel driver input class */

    PRINT_LINE_STAR;
    CkPrintf("Cpaimd-Charm-Driver input started \n");
    PRINT_LINE_DASH; CkPrintf("\n");
    Timer=CmiWallTimer();
    double phase1start=Timer;
    //numPes = 2048; 
    numPes=CkNumPes();
    config.readConfig(msg->argv[1],sim->nstates,sim->sizeX,sim->sizeY,sim->sizeZ,
                      sim->ntime,ibinary_opt,natm_nl,fftopt,numPes,natm_typ,
                      ees_eext_opt,sim->gen_wave,sim->ncoef, sim->cp_min_opt, sim->ngrid_eext_c);
    fakeTorus        = config.fakeTorus>0;
    CkPrintf("for numInstances %d numPes %d numPesPerInstance is %d \n",config.numInstances, config.numPes, config.numPesPerInstance);
    if(fakeTorus)
      numPes=config.torusDimNX * config.torusDimNY * config.torusDimNZ * config.torusDimNT;
    int numSfGrps    = config.numSfGrps;  // local copies are nice
    int doublePack   = config.doublePack;
    nchareG          = config.nchareG;
    sim->nchareG     = nchareG; 

    nstates          = config.nstates;    // globals : avail on all procs
    sizeX            = sim->sizeX;

    double newtime= CmiWallTimer();

    PRINT_LINE_DASH;
    CkPrintf("Cpaimd-Charm-Driver input completed in %g\n",newtime-Timer);
    PRINT_LINE_STAR; CkPrintf("\n");
    Timer=newtime;

//============================================================================
// Set user trace events for projections optimizations

     traceRegisterUserEvent("doRealFwFFT", doRealFwFFT_);

     traceRegisterUserEvent("GspaceFwFFT", GspaceFwFFT_);

     traceRegisterUserEvent("fwFFTGtoR0", fwFFTGtoR0_);
     traceRegisterUserEvent("fwFFTGtoRnot0", fwFFTGtoRnot0_);
     traceRegisterUserEvent("OrthoDGEMM1", OrthoDGEMM1_);
     traceRegisterUserEvent("GradCorrGGA", GradCorrGGA_);
     traceRegisterUserEvent("WhiteByrdFFTX", WhiteByrdFFTX_);
     traceRegisterUserEvent("doRealBwFFT", doRealBwFFT_);
     traceRegisterUserEvent("WhiteByrdFFTY", WhiteByrdFFTY_);
     traceRegisterUserEvent("WhiteByrdFFTZ", WhiteByrdFFTZ_);
     traceRegisterUserEvent("PostByrdfwFFTGtoR", PostByrdfwFFTGtoR_);
     traceRegisterUserEvent("RhoRtoGFFT", RhoRtoGFFT_);
     traceRegisterUserEvent("BwFFTRtoG", BwFFTRtoG_);
     traceRegisterUserEvent("OrthoDGEMM2", OrthoDGEMM2_);
     traceRegisterUserEvent("ByrdanddoFwFFTGtoR",ByrdanddoFwFFTGtoR_);
     traceRegisterUserEvent("eesHartExcG",eesHartExcG_);
     traceRegisterUserEvent("eesEwaldG",eesEwaldG_);
     traceRegisterUserEvent("eesAtmForcR",eesAtmForcR_);
     traceRegisterUserEvent("eesAtmBspline",eesAtmBspline_);
     traceRegisterUserEvent("eesZmatR",eesZmatR_);
     traceRegisterUserEvent("eesEnergyAtmForcR",eesEnergyAtmForcR_);
     traceRegisterUserEvent("eesProjG",eesProjG_);
     traceRegisterUserEvent("doNlFFTGtoR",doNlFFTGtoR_);
     traceRegisterUserEvent("doNlFFTRtoG",doNlFFTRtoG_);
     traceRegisterUserEvent("eesPsiForcGspace",eesPsiForcGspace_);
     traceRegisterUserEvent("enlMatrixCalc",enlMatrixCalc_);
     traceRegisterUserEvent("enlAtmForcCalc",enlAtmForcCalc_);
     traceRegisterUserEvent("enlForcCalc",enlForcCalc_);
     traceRegisterUserEvent("doEextFFTRtoG",doEextFFTRtoG_);
     traceRegisterUserEvent("doEextFFTGtoR",doEextFFTGtoR_);

     traceRegisterUserEvent("HartExcVksG",HartExcVksG_);
     traceRegisterUserEvent("divRhoVksGspace",divRhoVksGspace_);
     traceRegisterUserEvent("GspaceBwFFT", GspaceBwFFT_);
     traceRegisterUserEvent("DoFFTContribute", DoFFTContribute_);
     traceRegisterUserEvent("IntegrateModForces", IntegrateModForces_);
     traceRegisterUserEvent("Scalcmap", Scalcmap_);
     traceRegisterUserEvent("AcceptStructFact", AcceptStructFact_);
     traceRegisterUserEvent("doEextFFTGxtoRx", doEextFFTGxtoRx_);
     traceRegisterUserEvent("doEextFFTRytoGy", doEextFFTRytoGy_);
     traceRegisterUserEvent("doRhoFFTRytoGy", doRhoFFTRytoGy_);
     traceRegisterUserEvent("doRhoFFTGxtoRx", doRhoFFTGxtoRx_);

     traceRegisterUserEvent("GSProcnum", 10000);
     traceRegisterUserEvent("RSProcnum", 20000);
     traceRegisterUserEvent("SCProcnum", 30000);
     traceRegisterUserEvent("GHartAtmForcCopy",GHartAtmForcCopy_);
     traceRegisterUserEvent("GHartAtmForcSend",GHartAtmForcSend_);
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

    // bump all the INT_MAPs to the right size
    GSImaptable.resize(config.numInstances);
    RSImaptable.resize(config.numInstances);
    RPPImaptable.resize(config.numInstances);
    RhoGSImaptable.resize(config.numInstances);
    RhoRSImaptable.resize(config.numInstances);
    RhoGHartImaptable.resize(config.numInstances);
    RhoRHartImaptable.resize(config.numInstances);
    OrthoImaptable.resize(config.numInstances);
    OrthoHelperImaptable.resize(config.numInstances);
    AsymScalcImaptable.resize(config.numInstances);
    SymScalcImaptable.resize(config.numInstances);

    // bump all our proxy vecs to the right size
    UpairCalcID2.resize(config.numInstances);
    UpairCalcID1.resize(config.numInstances);
    UgSpacePlaneProxy.reserve(config.numInstances);
    UgSpaceDriverProxy.reserve(config.numInstances);
    UparticlePlaneProxy.reserve(config.numInstances);
    UrealParticlePlaneProxy.reserve(config.numInstances);
    UrealSpacePlaneProxy.reserve(config.numInstances);
    UrhoRealProxy.reserve(config.numInstances);
    UrhoGProxy.reserve(config.numInstances);
    UrhoRHartExtProxy.reserve(config.numInstances);
    UrhoGHartExtProxy.reserve(config.numInstances);
    UorthoProxy.reserve(config.numInstances);
    UatomsGrpProxy.reserve(config.numInstances);
    UegroupProxy.reserve(config.numInstances);
    UfftCacheProxy.reserve(config.numInstances);
    UsfCacheProxy.reserve(config.numInstances);
    UsfCompProxy.reserve(config.numInstances);
    UeesCacheProxy.reserve(config.numInstances);
    UorthoHelperProxy.reserve(config.numInstances);
    UplaneUsedByNLZ.reserve(config.numInstances);
    UlsRhoRealProxy.reserve(config.numInstances);
    UlsRhoGProxy.reserve(config.numInstances);

    excludePes=NULL;
    mainProxy=thishandle;
    // make one controller chare per instance
    instControllerProxy= CProxy_InstanceController::ckNew(config.numInstances);
    instControllerProxy.doneInserting();
//============================================================================    
// Create the multicast/reduction manager for array sections
// Create the parainfo group from sim
// Initialize chare arrays for real and g-space of states 


    int l=config.Gstates_per_pe;
    int m, pl, pm;
    pl = nstates / l;
    pm = config.numPesPerInstance / pl;
    if(pm == 0) {
      CkPrintf("Choose a larger Gstates_per_pe than %d such that { (no. of processors [%d] / no. of Instances [%d]) / (no. of states [%d] / Gstates_per_pe [%d]) } is > 0 \n", 
	l, config.numPes, config.numInstances, nstates, l);
      CkAssert( pm > 0);
    }
    m = config.nchareG / pm;

    planes_per_pe=m;
    if(planes_per_pe <= 0) {
      CkPrintf("Choose a Gstates_per_pe than %d such that { (no. of processors [%d] / no. of Instances [%d]) / (no. of states [%d] / Gstates_per_pe [%d]) } is > 0 and config.nchareG [%d] / pm [%d] {where pm = config.numInstances [%d]/ pl [%d] }  > 0\n", 
	l, config.numPes, config.numInstances, nstates, l, config.nchareG, pm, config.numPesPerInstance, pl);
      CkAssert( m > 0);
    }

    mCastGrpId = CProxy_CkMulticastMgr::ckNew(config.numMulticastMsgs);
    //    if(pm==0){CkAbort("Choose a larger Gstates_per_pe\n");}
    //    for(int i=0; i<nstates;i++){
    //      peUsedByNLZ.push_back(((i % config.Gstates_per_pe)*planes_per_pe)%nchareG);
    //    }//endfor

    if(config.torusMap==1) {
      PRINT_LINE_STAR; CkPrintf("\n");
      CkPrintf("         Topology Sensitive Mapping being done for RSMap, GSMap, ....\n");
      CkPrintf("            ......., PairCalc, RhoR, RhoG and RhoGHart .........\n\n");
      PRINT_LINE_STAR; CkPrintf("\n");
      CkPrintf("Initializing TopoManager\n");
      if(config.fakeTorus) {
	  topoMgr = new TopoManager(config.torusDimNX, config.torusDimNY, 
				    config.torusDimNZ, config.torusDimNT);
      }
      else {
	  topoMgr = new TopoManager();
      }
      CkPrintf("         Torus %d x %d x %d nodes %d x %d x %d VN %d DimNT %d .........\n", 
             topoMgr->getDimX(), topoMgr->getDimY(), topoMgr->getDimZ(),
             topoMgr->getDimNX(), topoMgr->getDimNY(), topoMgr->getDimNZ(),
             topoMgr->hasMultipleProcsPerNode(), topoMgr->getDimNT());
    }
    CkPrintf("Initializing PeList\n");
    
    PeList *gfoo=NULL;
    PeList *rfoo=NULL;

    if(!config.loadMapFiles && config.useCuboidMap)
      {
	if( config.numPesPerInstance % config.nchareG != 0)
	  {
	    CkPrintf("To use CuboidMap nchareG %d should be chosen as a factor of numprocs %d / numInstances %d = %d\n",config.nchareG, config.numPes, config.numInstances, config.numPesPerInstance);
	    CkExit();
	  }
	int procsPerPlane = config.numPesPerInstance / nchareG;
	int bx, by, bz;
        if(config.torusMap == 1) {
	  boxSize = procsPerPlane;
	  int order;

	  // correction to accomodate multiple instances
	  int dimNX, dimNY, dimNZ, dimNT;
	  int longDim = -1, maxD;

	  dimNX = topoMgr->getDimNX();
	  dimNY = topoMgr->getDimNY();
	  dimNZ = topoMgr->getDimNZ();
	  dimNT = topoMgr->getDimNT();
	  int x, y, z, x1 = dimNX, y1 = dimNY, z1 = dimNZ;

	  maxD = dimNX;
	  longDim = 1;

	  if (dimNY > maxD) { maxD = dimNY; longDim = 2; }
	  if (dimNZ > maxD) { maxD = dimNZ; longDim = 3; }

	  for(int i=0; i<config.numInstances; i++) {
	    switch(longDim) {
	      case 1:
		x = i*(maxD/config.numInstances); y = 0; z = 0;
		x1 = dimNX / config.numInstances;
		break;
	      case 2:
		x = 0; y = i*(maxD/config.numInstances); z = 0;
		y1 = dimNY / config.numInstances;
		break;
	      case 3:
		x = 0; y = 0; z = i*(maxD/config.numInstances);
		z1 = dimNZ / config.numInstances;
		break;
	    }
	    // printf("Origin %d %d %d Size %d %d %d\n", x, y, z, x1, y1, z1);
	  }
       
	  CkPrintf("         Box per Instance: nodes %d x %d x %d .........\n", x1, y1, z1);
	  if(findCuboid(bx, by, bz, order, x1, y1, z1, topoMgr->getDimNT(), boxSize, topoMgr->hasMultipleProcsPerNode()))
	  {
	    CkPrintf("Using %d, %d, %d dimensions for box %d mapping order %d\n", bx, by, bz, boxSize, order);
	    gfoo = new PeList(bx, by, bz, order, x1, y1, z1, dimNT);	// heap it
	  }
	  else
	  {
	    CkPrintf("no box for %d\n", boxSize);
	    config.useCuboidMap = 0;
	    gfoo = new PeList;				// heap it
	  }
	}
	else
	  gfoo = new PeList;				// heap it
      }
    else
      gfoo = new PeList;				// heap it
    if(!config.loadMapFiles && config.useCuboidMapRS)
      rfoo = new PeList(*gfoo);
    else
      rfoo = new PeList;				// heap it
    
    availGlobG = rfoo;
    availGlobR = gfoo;
    newtime = CmiWallTimer();
    CkPrintf("Pelist initialized in %g\n", newtime-Timer);
    // availGlobG->dump();
    Timer = newtime;

    // TODO timekeeper registerees will need to distinguish by instance
    // timekeeper itself doesn't care.
    TimeKeeperProxy = CProxy_TimeKeeper::ckNew();    

    /*
===============================================================================
Per Instance startup BEGIN
===============================================================================
    */
    // used as a signal during startup, to handle things, such as map
    // creation, which will only be done for the first instance.

    // maps will have a transform function to compute the placement
    // for instances after the first.
    
    for(int integral=0; integral< config.UberImax; integral++)
      {
      for(int kpoint=0; kpoint< config.UberJmax; kpoint++)
	{
	for(int temper=0; temper< config.UberKmax; temper++)
	  {
	    // for each new instance we need a new Uber Index
	    CkVec  <int>  peUsedBySF;
	    CkVec  <int>  peUsedByNLZ;
	    CkVec  <int>  planeUsedByNLZ;

	    UberIndex thisInstanceIndex(integral,kpoint,temper);
	    thisInstance=UberCollection(thisInstanceIndex);
	    UberAlles.push_back(thisInstance);
	    //============================================================================    
	    // Transfer parameters from physics to driver
	    //    read in atoms : create atoms group 
	    
	    // we will need a different one of these per instance
	    control_physics_to_driver(thisInstance);
	    
	    //============================================================================ 

	    // and then we make the usual set of chares to which we pass
	    // the Uber Index.
	    init_state_chares(natm_nl,natm_nl_grp_max,numSfGrps,doublePack,sim, thisInstance);

	    CmiNetworkProgressAfter(1);

	    int *usedProc= new int[config.numPesPerInstance];
	    memset(usedProc, 0, sizeof(int)*config.numPesPerInstance);
	    int charperpe = nstates/(config.numPesPerInstance);
	    if(nstates % config.numPesPerInstance != 0)  charperpe++;
	    if(charperpe<1) charperpe=1;
	    for(int state=0; state<nstates; state++) {
	      int plane = nchareG-1;
	      while(plane >= 0) {
		bool used = false;
		int thisstateplaneproc = GSImaptable[thisInstance.getPO()].get(state,plane) % config.numPesPerInstance;
		if(usedProc[thisstateplaneproc]>charperpe) {
		  used=true;
		}
		if(!used || plane==0) {
		  peUsedByNLZ.push_back(thisstateplaneproc);
		  planeUsedByNLZ.push_back(plane);
		  usedProc[thisstateplaneproc]++;
		  plane=-1;
		}
		plane--;
	      }
	    }
	    peUsedByNLZ.quickSort();
	    delete [] usedProc;
	    UpeUsedByNLZ.push_back(peUsedByNLZ);	   
	    UplaneUsedByNLZ.push_back(planeUsedByNLZ);
	    // CkPrintf("UplaneUsedByNLZ length now %d\n",UplaneUsedByNLZ.length());
	    // Create mapping classes for Paircalcular

	    //-------------------------------------------------------------
	    int indexSize = nchareG;

	    int* indexZ = new int[indexSize];
	    for(int i=0, count=0; i<nchareG; i++){
	      indexZ[count] = i;
	      count++;
	    }//endif

	    //============================================================================ 
	    // Initialize paircalculators for Psi and Lambda and ortho

	    init_pair_calculators( nstates,indexSize,indexZ,doublePack,sim, boxSize, thisInstance);
	    CmiNetworkProgressAfter(1);

	    //============================================================================ 
	    // Initialize the density chare arrays
	    init_rho_chares(sim, thisInstance);
	    CmiNetworkProgressAfter(1);
	    //============================================================================ 
	    // Initialize commlib strategies for later association and delegation
	    if(sim->ees_nloc_on)
	      init_eesNL_chares( natm_nl, natm_nl_grp_max, doublePack, excludePes, sim, thisInstance);
	    CmiNetworkProgressAfter(1);
	    delete [] indexZ;

	  } 
	}
      } // end of per instance init
    //============================================================================ 
    // Initialize commlib strategies for later association and delegation
    if(config.numInstances>0)
	CkPrintf("WARNING!!! Commlib does not work for multiple instances\n");

    init_commlib_strategies(sim->nchareRhoG, sim->sizeZ,nchareRhoRHart, thisInstance);

    TimeKeeperProxy.init();

//============================================================================
// clean up

    delete msg;
    delete sim;
    delete rfoo;
    delete gfoo;
    delete excludePes;
    //    delete availGlobR;
    //    delete availGlobG;

//============================================================================

    newtime=CmiWallTimer();
    PRINT_LINE_DASH;
    CkPrintf("Cpaimd-Charm-Driver setup phase completed in %g \n",newtime-phase1start);
    PRINT_LINE_STAR; CkPrintf("\n");
    PRINT_LINE_STAR; 
    PRINT_LINE_DASH;CkPrintf("\n");
    CkPrintf("user mem %d\n",CmiMemoryUsage());
    Timer=newtime;

//============================================================================
   }// end Main
//============================================================================

//============================================================================    
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================    
/**
 * Cleanup stuff in the hopes of getting clean valgrind
 */
main::~main(){
    if (config.useCommlib) {        
	if(config.usePairEtoM){
	}//endif
    }//endif
}//end routine
//============================================================================    


//============================================================================    
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================    
/**
 * Initialize paircalc1 Psi (sym) and paircalc2 Lambda (asym)
 */
//============================================================================    
void init_pair_calculators(int nstates, int indexSize, int *indexZ ,
                           int doublePack, CPcharmParaInfo *sim, int boxSize,
			   UberCollection thisInstance
			   )
//============================================================================    
  {// begin routine
//============================================================================    

  PRINT_LINE_STAR;
  PRINTF("Building Psi and Lambda Pair Calculators\n");
  PRINT_LINE_DASH;printf("\n");

//============================================================================    
  //-------------------------------------------------------------
  // Populate maptable for Symmetric Paircalculators
  Timer =CmiWallTimer();
  availGlobG->reset();
  bool maptype=true;
  int achunks=config.numChunksAsym;
  bool cp_need_orthoT= (sim->cp_min_opt==1) ? false: true;
  if(config.phantomSym)
    { // evil trickery to use asym map code for phantom sym
      maptype=false;
      achunks=config.numChunksSym; 
      
    }//endif
#ifdef USE_INT_MAP
  SymScalcImaptable[thisInstance.getPO()].buildMap(config.nchareG, config.nstates/config.sGrainSize, config.nstates/config.sGrainSize, achunks, config.sGrainSize);
#endif

  int success = 0;
  if(config.loadMapFiles) {
    int size[4];
    size[0] = config.nchareG; size[1] = config.nstates/config.sGrainSize;
    size[2] = config.nstates/config.sGrainSize; size[3] = achunks;
    MapFile *mf = new MapFile("SymScalcMap", 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, config.sGrainSize);
#ifdef USE_INT_MAP
    success = mf->loadMap("SymScalcMap", &SymScalcImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("SymScalcMap", &SymScalcmaptable);
#endif
    delete mf;
  }

  if(success == 0) {
#ifdef USE_INT_MAP
    SCalcMapTable symTable = SCalcMapTable(&SymScalcImaptable[thisInstance.getPO()], 
					 availGlobG, config.nstates,
                   config.nchareG, config.sGrainSize, maptype, 
		 config.scalc_per_plane, planes_per_pe, achunks, config.numChunksSym, &GSImaptable[thisInstance.getPO()], config.useCuboidMap, config.useCentroidMap, boxSize);
#else
    SCalcMapTable symTable = SCalcMapTable(&SymScalcmaptable, 
					 availGlobG, config.nstates,
                   config.nchareG, config.sGrainSize, maptype,
		 config.scalc_per_plane, planes_per_pe, achunks, config.numChunksSym, &GSmaptable, config.useCuboidMap, config.useCentroidMap, boxSize);
#endif
  }

  CProxy_SCalcMap scMap_sym = CProxy_SCalcMap::ckNew(CmiTrue, thisInstance);

  double newtime=CmiWallTimer();
  CkPrintf("SymScalcMap %d x %d x %d x %d created in %g\n",config.nchareG, config.nstates/config.sGrainSize, config.nstates/config.sGrainSize, config.numChunksSym, newtime-Timer);
  
  if(config.dumpMapFiles) {
    int size[4];
    size[0] = config.nchareG; size[1] = config.nstates/config.sGrainSize;
    size[2] = config.nstates/config.sGrainSize; size[3] = achunks;
    MapFile *mf = new MapFile("SymScalcMap", 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, config.sGrainSize);
#ifdef USE_INT_MAP
    mf->dumpMap(&SymScalcImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&SymScalcmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[4];
    size[0] = config.nchareG; size[1] = config.nstates/config.sGrainSize;
    size[2] = config.nstates/config.sGrainSize; size[3] = achunks;
    MapFile *mf = new MapFile("SymScalcMap_coord", 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, config.sGrainSize);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&SymScalcImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&SymScalcmaptable);
#endif
    delete mf;
  }

  //-------------------------------------------------------------
  // Populate maptable for Asymmetric Paircalculators
  Timer=CmiWallTimer();
  availGlobG->reset();
#ifdef USE_INT_MAP
  AsymScalcImaptable[thisInstance.getPO()].buildMap(config.nchareG, config.nstates/config.sGrainSize,config.nstates/config.sGrainSize, config.numChunksAsym, config.sGrainSize);
#endif

  success = 0;
  if(config.loadMapFiles) {
    int size[4];
    size[0] = config.nchareG; size[1] = config.nstates/config.sGrainSize;
    size[2] = config.nstates/config.sGrainSize; size[3] = config.numChunksAsym;
    MapFile *mf = new MapFile("AsymScalcMap", 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, config.sGrainSize);
#ifdef USE_INT_MAP
    success = mf->loadMap("AsymScalcMap", &AsymScalcImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("AsymScalcMap", &AsymScalcmaptable);
#endif
    delete mf;
  }

  if(success == 0) {
#ifdef USE_INT_MAP
    SCalcMapTable asymTable = SCalcMapTable(&AsymScalcImaptable[thisInstance.getPO()], availGlobG,config.nstates,
  	           config.nchareG, config.sGrainSize, CmiFalse, config.scalc_per_plane, 
                   planes_per_pe, config.numChunksAsym, config.numChunksSym, &GSImaptable[thisInstance.getPO()], config.useCuboidMap, config.useCentroidMap, boxSize);
#else
    SCalcMapTable asymTable = SCalcMapTable(&AsymScalcmaptable, availGlobG,config.nstates,
  	           config.nchareG, config.sGrainSize, CmiFalse, config.scalc_per_plane, 
                   planes_per_pe, config.numChunksAsym, config.numChunksSym, &GSmaptable, config.useCuboidMap, config.useCentroidMap,boxSize);
#endif
  }

  CProxy_SCalcMap scMap_asym = CProxy_SCalcMap::ckNew(CmiFalse, thisInstance);
  newtime=CmiWallTimer();
  CkPrintf("AsymScalcMap %d x %d x %d x %d created in %g\n\n",config.nchareG, config.nstates/config.sGrainSize, config.nstates/config.sGrainSize, config.numChunksAsym, newtime-Timer);
  Timer=newtime;

  if(config.dumpMapFiles) {
    int size[4];
    size[0] = config.nchareG; size[1] = config.nstates/config.sGrainSize;
    size[2] = config.nstates/config.sGrainSize; size[3] = config.numChunksAsym;
    MapFile *mf = new MapFile("AsymScalcMap", 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, config.sGrainSize);
#ifdef USE_INT_MAP
    mf->dumpMap(&AsymScalcImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&AsymScalcmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[4];
    size[0] = config.nchareG; size[1] = config.nstates/config.sGrainSize;
    size[2] = config.nstates/config.sGrainSize; size[3] = config.numChunksAsym;
    MapFile *mf = new MapFile("AsymScalcMap_coord", 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, config.sGrainSize);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&AsymScalcImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&AsymScalcmaptable);
#endif
    delete mf;
  }


  //============================================================================ 
  // initialize Ortho  now that we have the PC maps
  init_ortho_chares(nstates, indexSize, indexZ, thisInstance);
  CmiNetworkProgressAfter(1);

  CkGroupID scalc_sym_id  = scMap_sym.ckGetGroupID();
  CkGroupID scalc_asym_id = scMap_asym.ckGetGroupID();
  //-------------------------------------------------------------
  // Register the PCs
  int gsp_ep;
  int gsp_ep_tol;
  if(config.gSpaceSum)
    {
      gsp_ep =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_partialResultMsg;
      gsp_ep_tol =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsiV_partialResultMsg;
    }
  else
    {
      gsp_ep =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_CkReductionMsg;
      gsp_ep_tol =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsiV_CkReductionMsg;
    }
    //    CkGroupID symMcast = CProxy_CkMulticastMgr::ckNew(config.PCSpanFactor);

    for(int i=0; i< nchareG ;i++)
      mCastGrpIds.push_back(CProxy_CkMulticastMgr::ckNew(config.PCSpanFactor));
    //mCastGrpIds.push_back(symMcast);
   //symmetric AKA Psi
    orthomCastGrpId=(CProxy_CkMulticastMgr::ckNew(config.OrthoMcastSpanFactor));
    orthoRedGrpId=(CProxy_CkMulticastMgr::ckNew(config.OrthoRedSpanFactor));

#ifdef _CP_SUBSTEP_TIMING_
    UpairCalcID1[thisInstance.proxyOffset].forwardTimerID=keeperRegister("Sym Forward");
    UpairCalcID1[thisInstance.proxyOffset].backwardTimerID=keeperRegister("Sym Backward");
    UpairCalcID1[thisInstance.proxyOffset].beginTimerCB=  CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
    UpairCalcID1[thisInstance.proxyOffset].endTimerCB=  CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
#endif
    // CkPrintf("creating PC instance %d\n",thisInstance.proxyOffset);
    createPairCalculator(true, nstates, config.sGrainSize, indexSize, indexZ,  CkCallback(CkIndex_Ortho::start_calc(NULL), UorthoProxy[thisInstance.proxyOffset]), &(UpairCalcID1[thisInstance.proxyOffset]), gsp_ep, gsp_ep_tol, UgSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(), 1, &scalc_sym_id, doublePack, config.conserveMemory,config.lbpaircalc, config.psipriority, mCastGrpIds, orthomCastGrpId, orthoRedGrpId, config.numChunksSym, config.orthoGrainSize,  config.PCCollectTiles, config.PCstreamBWout, config.PCdelayBWSend, config.PCstreamFWblock, config.usePairDirectSend, config.gSpaceSum, config.gsfftpriority, config.phantomSym, config.useBWBarrier, config.gemmSplitFWk, config.gemmSplitFWm, config.gemmSplitBW,false, thisInstance.proxyOffset);

    CkArrayIndex2D myindex(0, 0);
    if(config.gSpaceSum)
      gsp_ep = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_partialResultMsg;
    else
      gsp_ep = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_CkReductionMsg;
    int myPack=0;
    CkVec <CkGroupID> mCastGrpIdsA;
    for(int i=0; i< nchareG ;i++)
      mCastGrpIdsA.push_back(CProxy_CkMulticastMgr::ckNew(config.PCSpanFactor));
      //asymmetric AKA Lambda AKA Gamma
#ifdef _CP_SUBSTEP_TIMING_
    UpairCalcID2[thisInstance.proxyOffset].forwardTimerID=keeperRegister("Asym Forward");
    UpairCalcID2[thisInstance.proxyOffset].backwardTimerID=keeperRegister("Asym Backward");
    UpairCalcID2[thisInstance.proxyOffset].beginTimerCB= CkCallback(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
    UpairCalcID2[thisInstance.proxyOffset].endTimerCB=  CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
#endif

    createPairCalculator(false, nstates,  config.sGrainSize, indexSize, indexZ,CkCallback(CkIndex_CP_State_GSpacePlane::acceptAllLambda(NULL), myindex, UgSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID()), &(UpairCalcID2[thisInstance.proxyOffset]), gsp_ep, 0, UgSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(), 1, &scalc_asym_id, myPack, config.conserveMemory,config.lbpaircalc, config.lambdapriority, mCastGrpIdsA, orthomCastGrpId, orthoRedGrpId,config.numChunksAsym, config.lambdaGrainSize,  config.PCCollectTiles, config.PCstreamBWout, config.PCdelayBWSend, config.PCstreamFWblock, config.usePairDirectSend, config.gSpaceSum, config.lambdapriority+2, false, config.useBWBarrier, config.gemmSplitFWk, config.gemmSplitFWm, config.gemmSplitBW, cp_need_orthoT,thisInstance.proxyOffset);
    

//============================================================================ 
   }//end routine
//============================================================================ 



//============================================================================ 
/**
 * Initialize Commlib communication strategies  
 */ 
//============================================================================ 
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_commlib_strategies(int numRhoG, int numReal, int numRhoRhart, UberCollection thisInstance){
//============================================================================

  PRINT_LINE_STAR;
  PRINTF("Building Commlib strategies\n");
  PRINT_LINE_DASH;printf("\n");
  int i = 0;
  int j = 0;
  int rhoGhelpers = config.rhoGHelpers;
  int numRhoGHart = rhoGhelpers*numRhoG;
  int nchareHartAtmT = config.nchareHartAtmT;
//============================================================================
  if (config.useCommlib) {        

    /*
    //StreamingStrategy *cmstrat = new StreamingStrategy(0.1,10);
    MeshStreamingStrategy *cmstrat = new MeshStreamingStrategy(1,5);
    gAsymInstance= ComlibRegister(cmstrat);    

    //StreamingStrategy *csymstrat = new StreamingStrategy(0.1,10);
    MeshStreamingStrategy *csymstrat = new MeshStreamingStrategy(1,5);
    gSymInstance= ComlibRegister(csymstrat);    
    */
    //================================================================
    CkArrayIndexMax *rhoGElements=NULL;
    CkArrayIndexMax *rhoRealElements = NULL;
    if(config.useRInsRhoGP)
      {
	CkPrintf("Making real_strategy with :");
	CkPrintf("src numReal %d dest numRhoG %d and numHartG %d numHartR %d\n",
		 numReal,numRhoG,numRhoGHart,numRhoRhart);
	//--------------------------------------------------------------
	//  For rho(r) to rho(g)
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  CkArrayIndex2D idx2d(i,0);
	  rhoGElements[i] = idx2d;
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  CkArrayIndex2D idx2d(i,0);
	  rhoRealElements[i] = idx2d; 
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *real_strat;
#else
	Strategy *real_strat;
#endif      
	real_strat= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG, rhoGElements);
	commRealInstance= ComlibRegister(real_strat);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    //  For drho(r)/dx to igx*rho(g)
    if(config.useRInsIGXRhoGP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *real_strat_igx;
#else
	Strategy *real_strat_igx;
#endif
	real_strat_igx= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG,rhoGElements);
	commRealIGXInstance= ComlibRegister(real_strat_igx);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    if(config.useRInsIGYRhoGP)
      {
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
#ifdef OLD_COMMLIB
	CharmStrategy *real_strat_igy;
#else
	Strategy *real_strat_igy;
#endif
	real_strat_igy= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG,rhoGElements);
	commRealIGYInstance= ComlibRegister(real_strat_igy);

	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    //  For drho(r)/dz to igz*rho(g)
    if(config.useRInsIGZRhoGP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}//endfor
#ifdef OLD_COMMLIB        
	CharmStrategy *real_strat_igz;
#else
	Strategy *real_strat_igz;
#endif
	real_strat_igz= new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(),
	   numReal, rhoRealElements, numRhoG,rhoGElements);
	commRealIGZInstance= ComlibRegister(real_strat_igz);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // For hartree-Ext(g) to vks(r)
    if(config.useGHartInsRhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoGHart*nchareHartAtmT];
	for (j= 0; j < nchareHartAtmT; j++) {
	  for (i = 0; i < numRhoGHart; i++) {
	    rhoGElements[i+j*numRhoGHart] = CkArrayIndex2D(i,j);
	  }//endfor
	}
	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *gstrathart = new EachToManyMulticastStrategy
#else
	Strategy *gstrathart = new EachToManyMulticastStrategy
#endif	

	  (USE_DIRECT, UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoGHart, rhoGElements, numReal, rhoRealElements);
	commGHartInstance = ComlibRegister(gstrathart);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // vks(g), igxrho igyrho igzrho and white byrd to g-space
    if(config.useGIns0RhoRP)
      {

	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}//endfor

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}//endfor
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat0;
#else
	Strategy *gstrat0;
#endif
	gstrat0 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance0 = ComlibRegister(gstrat0);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // rhog sends to rhor : div_x rho
    if(config.useGIns1RhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat1;
#else
	Strategy *gstrat1;
#endif

	gstrat1 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance1 = ComlibRegister(gstrat1);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // rhog sends to rhor : div_y rho
    if(config.useGIns2RhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat2;
#else
	Strategy *gstrat2;
#endif

	gstrat2 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance2 = ComlibRegister(gstrat2);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // rhog sends to rhor : div_z rho
    if(config.useGIns3RhoRP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstrat3;
#else
	Strategy *gstrat3;
#endif

	gstrat3 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGInstance3 = ComlibRegister(gstrat3);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // For white-byrd :  rhog sends to rhor : white byrd
    if(config.useGByrdInsRhoRBP)
      {
	rhoGElements = new CkArrayIndexMax[numRhoG];
	for (i = 0; i < numRhoG; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}

	rhoRealElements = new  CkArrayIndexMax[numReal];
	for(i = 0; i < numReal; i++) {
	  rhoRealElements[i] = CkArrayIndex2D(i,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstratByrd;
#else
	Strategy *gstratByrd;
#endif

	gstratByrd = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoG, rhoGElements, numReal, rhoRealElements);
	commGByrdInstance = ComlibRegister(gstratByrd);
	delete [] rhoGElements;
	delete [] rhoRealElements;
      }
    //--------------------------------------------------------------
    // For rhogHart send to rhorHart SF atmtyp
    if(config.useGHartInsRHart)
      {
	rhoGElements = new CkArrayIndexMax[numRhoGHart];
	for (i = 0; i < numRhoGHart; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}

	rhoRealElements = new  CkArrayIndexMax[numRhoRhart];
	for(i = 0; i < numRhoRhart; i++) {
	  rhoRealElements[i] = CkArrayIndex3D(i,0,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstratEext0;
#else
	Strategy *gstratEext0;
#endif

	gstratEext0 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoGHart, rhoGElements, numRhoRhart, rhoRealElements);

	commGHartRHartIns0 = ComlibRegister(gstratEext0);
      }
    //--------------------------------------------------------------
    // For rhogHart send to rhoRHart SF tot
    if(config.useGHartInsRHart)
      {
	rhoGElements = new CkArrayIndexMax[numRhoGHart];
	for (i = 0; i < numRhoGHart; i++) {
	  rhoGElements[i] = CkArrayIndex2D(i,0);
	}

	rhoRealElements = new  CkArrayIndexMax[numRhoRhart];
	for(i = 0; i < numRhoRhart; i++) {
	  rhoRealElements[i] = CkArrayIndex3D(i,0,0);
	}
#ifdef OLD_COMMLIB
	CharmStrategy *gstratEext1;
#else
	Strategy *gstratEext1;
#endif

	gstratEext1 = new EachToManyMulticastStrategy
	  (USE_DIRECT, UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), 
	   numRhoGHart, rhoGElements, numRhoRhart, rhoRealElements);

	commGHartRHartIns1 = ComlibRegister(gstratEext1);
      }
    //--------------------------------------------------------------
    // For rhoRHart send to rhoGHart SF
    if(config.useRHartInsGHart){
      rhoGElements = new CkArrayIndexMax[numRhoGHart];
      for (i = 0; i < numRhoGHart; i++) {
	rhoGElements[i] = CkArrayIndex2D(i,0);
      }//endfor

      rhoRealElements = new  CkArrayIndexMax[numRhoRhart];
      for(i = 0; i < numRhoRhart; i++) {
	rhoRealElements[i] = CkArrayIndex3D(i,0,0);
      }//endfor
#ifdef OLD_COMMLIB
      CharmStrategy *real_strat_eext;
#else
	Strategy *real_strat_eext;
#endif
        
      real_strat_eext = new EachToManyMulticastStrategy
	(USE_DIRECT, UrhoRHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(), UrhoGHartExtProxy[thisInstance.proxyOffset].ckGetArrayID(),
	 numRhoRhart, rhoRealElements, numRhoGHart,rhoGElements);
      commRHartGHartIns  = ComlibRegister(real_strat_eext);
    }

  }//endif : use commlib

    //============================================================================
    // Real state space to gspace state and particle plane comm.

    if (config.useCommlibMulticast) {
      DirectMulticastStrategy *dstrat = new DirectMulticastStrategy
	(UrealSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
        
      RingMulticastStrategy *rstrat = new RingMulticastStrategy
	(UrealSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
        
      RingMulticastStrategy *r1strat = new RingMulticastStrategy
	(UparticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);

#ifdef OLD_COMMLIB
      MultiRingMulticast *mr1strat = new MultiRingMulticast
	(UparticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
#else
      MultiRingMulticastStrategy *mr1strat = new MultiRingMulticastStrategy
	();
#endif 
      DirectMulticastStrategy *ppdstrat = new DirectMulticastStrategy
	(UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
 
      RingMulticastStrategy *pprstrat = new RingMulticastStrategy
	(UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
#ifdef OLD_COMMLIB 
      MultiRingMulticast *ppmr1strat = new MultiRingMulticast
	(UrealParticlePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(),1);
#else
      MultiRingMulticastStrategy *ppmr1strat = new MultiRingMulticastStrategy
	(1);
#endif

      //multiring should be good on large runs, but not on BG/L
      if(CkNumNodes()>64){
	mcastInstance=ComlibRegister(dstrat);
	mcastInstancePP=ComlibRegister(mr1strat);
	mcastInstanceRPP=ComlibRegister(ppdstrat);
	mcastInstancemRPP=ComlibRegister(ppmr1strat);
      }else{
	mcastInstance=ComlibRegister(rstrat);
	mcastInstancePP=ComlibRegister(r1strat);
	mcastInstanceRPP=ComlibRegister(pprstrat);
	mcastInstancemRPP=ComlibRegister(ppmr1strat);
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
void init_ortho_chares(int nstates, int indexSize, int *indexZ, UberCollection thisInstance) {
//============================================================================

  PRINT_LINE_STAR;
  PRINTF("Building Ortho Chares\n");
  PRINT_LINE_DASH;printf("\n");
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
  if(config.loadMapFiles) {
    int size[2];
    size[0] = size[1] = nstates/config.orthoGrainSize;
    MapFile *mf = new MapFile("OrthoMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    success = mf->loadMap("OrthoMap", &OrthoImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("OrthoMap", &Orthomaptable);
#endif
    delete mf;
  }
  PeList *avail= new PeList();
  if(success == 0) {
#ifdef USE_INT_MAP
    OrthoMapTable Otable = OrthoMapTable(&OrthoImaptable[thisInstance.getPO()], avail, nstates, config.orthoGrainSize, &AsymScalcImaptable[thisInstance.getPO()], config.nchareG, config.numChunks, config.sGrainSize, excludePes);
#else
    OrthoMapTable Otable = OrthoMapTable(&Orthomaptable, avail, nstates, config.orthoGrainSize, &AsymScalcmaptable, config.nchareG, config.numChunks, config.sGrainSize, excludePes);
#endif
  }
  double newtime=CmiWallTimer();
  CkPrintf("OrthoMap created in %g\n\n", newtime-Timer);

  //CProxy_OrthoMap orthoMap = CProxy_OrthoMap::ckNew(chunks, nOrtho, stride);
  orthoMap = CProxy_OrthoMap::ckNew(thisInstance);
  CkArrayOptions orthoOpts;
  orthoOpts.setMap(orthoMap);

  UorthoProxy.push_back( CProxy_Ortho::ckNew(orthoOpts));
#ifdef OLD_COMMLIB
  CharmStrategy *multistrat = new DirectMulticastStrategy(UorthoProxy[thisInstance.proxyOffset].ckGetArrayID());
#else
  Strategy *multistrat = new DirectMulticastStrategy(UorthoProxy[thisInstance.proxyOffset].ckGetArrayID());
#endif
  orthoInstance=ComlibRegister(multistrat);

  CkCallback ocb= CkCallback(CkIndex_Ortho::collect_error(NULL), UorthoProxy[thisInstance.proxyOffset](0, 0));
  UorthoProxy[thisInstance.proxyOffset].ckSetReductionClient(&ocb);
    
  if(config.dumpMapFiles) {
    int size[2];
    size[0] = size[1] = nstates/config.orthoGrainSize;
    MapFile *mf = new MapFile("OrthoMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&OrthoImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&Orthomaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = size[1] = nstates/config.orthoGrainSize;
    MapFile *mf = new MapFile("OrthoMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&OrthoImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&Orthomaptable);
#endif
    delete mf;
  }

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

  /* create matrix multiplication objects */
  CLA_Matrix_interface matA1, matB1, matC1;
  CLA_Matrix_interface matA2, matB2, matC2;
  CLA_Matrix_interface matA3, matB3, matC3;

  CkCallback ortho_ready_cb = CkCallback(CkIndex_Ortho::all_ready(),
   UorthoProxy[thisInstance.proxyOffset](0, 0));

  make_multiplier(&matA1, &matB1, &matC1, UorthoProxy[thisInstance.proxyOffset], UorthoProxy[thisInstance.proxyOffset], UorthoProxy[thisInstance.proxyOffset],
   nstates, nstates, nstates, config.orthoGrainSize, config.orthoGrainSize,
   config.orthoGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
   mCastGrpId, MM_ALG_2D, config.gemmSplitOrtho);
  if(config.useOrthoHelpers)
    {
#ifdef USE_INT_MAP
      OrthoHelperImaptable[thisInstance.getPO()].buildMap(nstates/config.orthoGrainSize, nstates/config.orthoGrainSize);
#endif
      double Timer=CmiWallTimer();

      success = 0;
      if(config.loadMapFiles) {
	int size[2];
	size[0] = size[1] = nstates/config.orthoGrainSize;
	MapFile *mf = new MapFile("OrthoHelperMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
	success = mf->loadMap("OrthoHelperMap", &OrthoHelperImaptable[thisInstance.getPO()]);
#else
	success = mf->loadMap("OrthoHelperMap", &OrthoHelpermaptable);
#endif
	delete mf;
      }

      if(success == 0) {
#ifdef USE_INT_MAP
	OrthoHelperMapTable OHtable = OrthoHelperMapTable(&OrthoHelperImaptable[thisInstance.getPO()], nstates, config.orthoGrainSize, &OrthoImaptable[thisInstance.getPO()], avail, excludePes);
#else
	OrthoHelperMapTable OHtable = OrthoHelperMapTable(&OrthoHelperImaptable, nstates, config.orthoGrainSize, &Orthomaptable, avail, excludePes);
#endif
      }
      double newtime=CmiWallTimer();
      CkPrintf("OrthoHelperMap created in %g\n", newtime-Timer);

      CProxy_OrthoHelperMap orthoHMap = CProxy_OrthoHelperMap::ckNew(thisInstance);
      CkArrayOptions orthoHOpts;
      orthoHOpts.setMap(orthoHMap);

      if(config.dumpMapFiles) {
	int size[2];
	size[0] = size[1] = nstates/config.orthoGrainSize;
	MapFile *mf = new MapFile("OrthoHelperMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
	mf->dumpMap(&OrthoHelperImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
	mf->dumpMap(&OrthoHelpermaptable);
#endif
	delete mf;
      }
      if(config.dumpMapCoordFiles) {
	int size[2];
	size[0] = size[1] = nstates/config.orthoGrainSize;
	MapFile *mf = new MapFile("OrthoHelperMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
	mf->dumpMapCoords(&OrthoHelperImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
	mf->dumpMapCoords(&OrthoHelpermaptable);
#endif
	delete mf;
      }
      UorthoHelperProxy.push_back( CProxy_OrthoHelper::ckNew(orthoHOpts));
      make_multiplier(&matA2, &matB2, &matC2, UorthoHelperProxy[thisInstance.proxyOffset], UorthoHelperProxy[thisInstance.proxyOffset], UorthoHelperProxy[thisInstance.proxyOffset],
		      nstates, nstates, nstates, config.orthoGrainSize, config.orthoGrainSize,
		      config.orthoGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, 
		      ortho_ready_cb,	mCastGrpId, MM_ALG_2D, config.gemmSplitOrtho);
    }
  else  //no helpers
  {
    make_multiplier(&matA2, &matB2, &matC2, UorthoProxy[thisInstance.proxyOffset], UorthoProxy[thisInstance.proxyOffset], UorthoProxy[thisInstance.proxyOffset],
	nstates, nstates, nstates, config.orthoGrainSize, config.orthoGrainSize,
	config.orthoGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
	mCastGrpId, MM_ALG_2D, config.gemmSplitOrtho);
  }

  make_multiplier(&matA3, &matB3, &matC3, UorthoProxy[thisInstance.proxyOffset], UorthoProxy[thisInstance.proxyOffset], UorthoProxy[thisInstance.proxyOffset],
   nstates, nstates, nstates, config.orthoGrainSize, config.orthoGrainSize,
   config.orthoGrainSize, 1, 1, 1, ortho_ready_cb, ortho_ready_cb, ortho_ready_cb,
   mCastGrpId, MM_ALG_2D, config.gemmSplitOrtho);

  int timekeep=keeperRegister("Ortho S to T");
  int maxorthoindex=(nstates/config.orthoGrainSize-1);
  int maxorthostateindex=(nstates/config.orthoGrainSize-1) * config.orthoGrainSize;
  for (int s1 = 0; s1 <= maxorthostateindex; s1 += config.orthoGrainSize)
    for (int s2 = 0; s2 <= maxorthostateindex; s2 += config.orthoGrainSize) {
      int indX = s1 / config.orthoGrainSize;
      int indY = s2 / config.orthoGrainSize;
      indX = (indX>maxorthoindex) ? maxorthoindex : indX;
      indY = (indY>maxorthoindex) ? maxorthoindex : indY;

      UorthoProxy[thisInstance.proxyOffset](indX, indY).insert(config.orthoGrainSize, config.orthoGrainSize,
      matA1, matB1, matC1, matA2, matB2, matC2, matA3, matB3, matC3,timekeep, thisInstance);
      if(config.useOrthoHelpers)
      {
	UorthoHelperProxy[thisInstance.proxyOffset](indX, indY).insert(config.orthoGrainSize, config.orthoGrainSize,
		   matA2, matB2, matC2, thisInstance);
      }
    }
  UorthoProxy[thisInstance.proxyOffset].doneInserting();
  if(config.useOrthoHelpers)
    UorthoHelperProxy[thisInstance.proxyOffset].doneInserting();
    
  UorthoProxy[thisInstance.proxyOffset].makeSections(indexSize, indexZ);
  delete avail;
  delete excludePes;
//============================================================================
  }//end routine 
//============================================================================




//============================================================================
//Create the array elements for the GSpace, Particle and Real Space planes
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_state_chares(int natm_nl,int natm_nl_grp_max,int numSfGrps,
                       int doublePack, CPcharmParaInfo *sim, UberCollection thisInstance)
//============================================================================
   { //begin routine 
//============================================================================
/*
 * Set up the map variables used to control the location of the
 * 2d-array elements over processors   
*/
//============================================================================
// Useful Local variables


  int ngrida          = sim->sizeX; 
  int ngridb          = sim->sizeY; 
  int ngridc          = sim->sizeZ;

  int ngridaNl        = sim->ngrid_nloc_a;
  int ngridbNl        = sim->ngrid_nloc_b;
  int ngridcNl        = sim->ngrid_nloc_c;
  int ees_nonlocal_on = sim->ees_nloc_on;

  int ngridaEext      = sim->ngrid_eext_a;
  int ngridbEext      = sim->ngrid_eext_b;
  int ngridcEext      = sim->ngrid_eext_c;
  int ees_eext_on     = sim->ees_eext_on;

  int nchareG         = sim->nchareG;
  int nchareR         = sim->sizeZ;
  int nchareRPP       = ngridcNl;

  int numIterNL       = sim->nlIters;

  int nchareRhoG      = sim->nchareRhoG;
  int rhoGHelpers     = config.rhoGHelpers;
  int nchareGHart     = rhoGHelpers*nchareRhoG;
  int nchareRHart     = ngridcEext;

  int Rstates_per_pe  = config.Rstates_per_pe;
  int Gstates_per_pe  = config.Gstates_per_pe;
  int sGrainSize      = config.sGrainSize; 
  int numChunks       = config.numChunks;

  //Need our maps and groups to exist before anyone tries to use them

  //============================================================================

  //============================================================================
  //need some groups to exist before we kick off the state which use them
  //--------------------------------------------------------------------------------
  // Groups : no placement required 

  UsfCacheProxy.push_back( CProxy_StructFactCache::ckNew(numSfGrps,natm_nl,natm_nl_grp_max, thisInstance));
  if(firstInstance) CkPrintf("created sfcache proxy\n");
  UsfCompProxy.push_back(CProxy_StructureFactor::ckNew());
  if(firstInstance) CkPrintf("created sfcomp proxy\n");
  UeesCacheProxy.push_back(CProxy_eesCache::ckNew(nchareRPP,nchareG,nchareRHart,nchareGHart,
					     nstates,nchareRhoG, thisInstance));
  if(firstInstance) CkPrintf("created eescache proxy\n");

  int nchareRRhoTot  = nchareR*(config.rhoRsubplanes);
  int nchareRHartTot = nchareRHart*(config.rhoRsubplanes);

  int *numGState     = sim->nlines_per_chareG;
  int *numGNL        = sim->nlines_per_chareG;
  int *numGRho       = sim->nlines_per_chareRhoG;
  int *numGEext      = sim->nlines_per_chareRhoGEext;

  if(firstInstance)
    {
      int *numRXState    = new int [nchareR];
      int *numRYState    = new int [nchareR];
      int *numRXNL       = new int [nchareR];
      int *numRYNL       = new int [nchareR];
      for(int i=0;i<nchareR;i++){
	numRXState[i] = sim->sizeY;
	numRYState[i] = sim->nplane_x;
	numRXNL[i]    = ngridbNl;
	numRYNL[i]    = sim->nplane_x;
      }//endfor
      int *numRXRho      = new int [nchareRRhoTot];
      int *numRYRho      = new int [nchareRRhoTot];
      int *numRXEext     = new int [nchareRHartTot];
      int *numRYEext     = new int [nchareRHartTot];
      int *numSubGx      = sim->numSubGx;
      create_Rho_fft_numbers(nchareR,nchareRHart,config.rhoRsubplanes,
			     sim->nplane_rho_x,sim->sizeY,ngridbEext,
			     numRXRho,numRYRho,numRXEext,numRYEext,numSubGx);

      UfftCacheProxy.push_back(CProxy_FFTcache::ckNew(
					     sim->sizeX,sim->sizeY,sim->sizeZ,
					     ngridaEext,ngridbEext,ngridcEext,ees_eext_on,
					     ngridaNl,  ngridbNl,  ngridcNl,  ees_nonlocal_on, 
					     sim->nlines_max, sim->nlines_max_rho,
					     config.nchareG,nchareR,
					     config.nchareG,nchareRPP, 
					     nchareRhoG,    nchareR,    nchareRRhoTot,
					     nchareGHart,   nchareRHart,nchareRHartTot,
					     numGState,     numRXState, numRYState,
					     numGNL,        numRXNL,    numRYNL,
					     numGRho,       numRXRho,   numRYRho,
					     numGEext,      numRXEext,  numRYEext,
					     config.fftopt,config.fftprogresssplitReal,config.fftprogresssplit,
					     config.rhoRsubplanes, thisInstance));
      CkPrintf("created fftcache proxy\n");
      delete [] numRXState;
      delete [] numRYState;
      delete [] numRXNL;
      delete [] numRYNL;
      delete [] numRXRho;
      delete [] numRYRho;
      delete [] numRXEext;
      delete [] numRYEext;
    }      

  //============================================================================
  // Instantiate the Chares with placement determined by the maps

  //---------------------------------------------------------------------------
  // state g-space

  PRINT_LINE_STAR;
  PRINTF("Building G-space (%d %d) and R-space (%d %d/%d) state Chares\n",
          nstates, nchareG, nstates, nchareR, nchareRPP);
  PRINT_LINE_DASH;printf("\n");
  availGlobG->reset();
  double newtime=CmiWallTimer();

  int numInst = thisInstance.getPO();
    int x, y, z, x1, y1, z1;

  // correction to accomodate multiple instances
  if(config.torusMap == 1) {
    int dimNX, dimNY, dimNZ;
    int longDim = -1, maxD;

    x1 = dimNX = topoMgr->getDimNX();
    y1 = dimNY = topoMgr->getDimNY();
    z1 = dimNZ = topoMgr->getDimNZ();

    maxD = dimNX;
    longDim = 1;

    if (dimNY > maxD) { maxD = dimNY; longDim = 2; }
    if (dimNZ > maxD) { maxD = dimNZ; longDim = 3; }

    switch(longDim) {
      case 1:
	x = numInst*(maxD/config.numInstances); y = 0; z = 0;
	x1 = dimNX / config.numInstances;
	break;
      case 2:
	x = 0; y = numInst*(maxD/config.numInstances); z = 0;
	y1 = dimNY / config.numInstances;
	break;
      case 3:
	x = 0; y = 0; z = numInst*(maxD/config.numInstances);
	z1 = dimNZ / config.numInstances;
	break;
    }
  }

  GSImaptable[numInst].buildMap(nstates, nchareG);
  if (firstInstance) {
    int success = 0;
    if(config.loadMapFiles) {
      int size[2];
      size[0] = nstates; size[1] = nchareG;
      MapFile *mf = new MapFile("GSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      success = mf->loadMap("GSMap", &GSImaptable[0]);
#else
      success = mf->loadMap("GSMap", &GSmaptable);
#endif
    delete mf;
    }
    
    if(success == 0) {
#ifdef USE_INT_MAP
      GSMapTable gsTable = GSMapTable(&GSImaptable[0], &GSImaptable[numInst], availGlobG, 
		nchareG, nstates, Gstates_per_pe, config.useCuboidMap, numInst, x, y, z);
#else
      GSMapTable gsTable = GSMapTable(&GSmaptable, availGlobG, nchareG,
				      nstates, Gstates_per_pe, config.useCuboidMap);
#endif
    }

    CkPrintf("GSMap created in %g\n", newtime-Timer);
  } else {
      GSMapTable gsTable = GSMapTable(&GSImaptable[0], &GSImaptable[numInst], availGlobG, 
		nchareG, nstates, Gstates_per_pe, config.useCuboidMap, numInst, x, y, z);
  }
  
  // there is only one IntMap per chare type, but each instance has
  // its own map group
  CProxy_GSMap gsMap = CProxy_GSMap::ckNew(thisInstance);


  //  CkArrayOptions gSpaceOpts(nstates,nchareG);
  CkArrayOptions gSpaceOpts(nstates,nchareG);
  string forwardname("GSpaceForward");
  std::ostringstream fwdstrm;
  fwdstrm << forwardname << "." << thisInstance.idxU.x << "." << thisInstance.idxU.y << "." << thisInstance.idxU.z; 
  int gforward=keeperRegister(fwdstrm.str());
  string backwardname("GSpaceBackward");
  std::ostringstream bwdstrm;
  bwdstrm << backwardname << "." << thisInstance.idxU.x << "." << thisInstance.idxU.y << "." << thisInstance.idxU.z; 
  int gbackward=keeperRegister(bwdstrm.str());
  gSpaceOpts.setMap(gsMap);
  UgSpacePlaneProxy.push_back(CProxy_CP_State_GSpacePlane::ckNew(
                     sizeX,  1, 1, sGrainSize,  gforward, 
		     gbackward, thisInstance, gSpaceOpts));
  UgSpacePlaneProxy[thisInstance.proxyOffset].doneInserting();
  // CkPrintf("{%d} main uGSpacePlaneProxy[%d] is %d\n",thisInstance.proxyOffset,thisInstance.proxyOffset,CkGroupID(UgSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID()).idx);

 //--------------------------------------------------------------------------------
 // Bind the GSpaceDriver array to the GSpacePlane array so that they migrate together
 CkArrayOptions gspDriverOpts(nstates,nchareG);
 gspDriverOpts.bindTo(UgSpacePlaneProxy[thisInstance.proxyOffset]);
 UgSpaceDriverProxy.push_back( CProxy_GSpaceDriver::ckNew(thisInstance,gspDriverOpts) );
 UgSpaceDriverProxy[thisInstance.proxyOffset].doneInserting();
 //--------------------------------------------------------------------------------
 // We bind the particlePlane array to the gSpacePlane array migrate together

  //  CkArrayOptions particleOpts(nstates,nchareG);
  CkArrayOptions particleOpts(nstates,nchareG);
  particleOpts.setMap(gsMap); // the maps for both the arrays are the same
  particleOpts.bindTo(UgSpacePlaneProxy[thisInstance.proxyOffset]);
  UparticlePlaneProxy.push_back(CProxy_CP_State_ParticlePlane::ckNew(
                       nchareG, sim->sizeY, sim->sizeZ,ngridaNl,ngridbNl,ngridcNl,
                       1, numSfGrps, natm_nl, natm_nl_grp_max, nstates, 
                       nchareG, Gstates_per_pe, numIterNL, ees_nonlocal_on, 
		       thisInstance,  particleOpts));
  UparticlePlaneProxy[thisInstance.proxyOffset].doneInserting();

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareG;
    MapFile *mf = new MapFile("GSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&GSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&GSmaptable);
#endif
    delete mf;
  }

  if(config.dumpMapCoordFiles) {
    // if someone wants to dump nontopo maps, they'll need a manager
    if(topoMgr == NULL)
      topoMgr = new TopoManager();

    int size[2];
    size[0] = nstates; size[1] = nchareG;
    MapFile *mf = new MapFile("GSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&GSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&GSmaptable);
#endif
    delete mf;
  }

  //---------------------------------------------------------------------------
  // state r-space

  // correction to accomodate multiple instances
  RSImaptable[numInst].buildMap(nstates, nchareR);

  if(firstInstance) {
    Timer=CmiWallTimer();
    availGlobR->reset();

    int success = 0;
    if(config.loadMapFiles) {
      int size[2];
      size[0] = nstates; size[1] = nchareR;
      MapFile *mf = new MapFile("RSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      success = mf->loadMap("RSMap", &RSImaptable[0]);
#else
      success = mf->loadMap("RSMap", &RSmaptable);
#endif
      delete mf;
    }
    
    if(success == 0) {
#ifdef USE_INT_MAP
      RSMapTable RStable= RSMapTable(&RSImaptable[0], &RSImaptable[numInst], availGlobR, 
			    nstates, nchareR, Rstates_per_pe, config.useCuboidMapRS, 
			    &GSImaptable[0], config.nchareG, numInst, x, y, z);
#else
      RSMapTable RStable= RSMapTable(&RSmaptable, availGlobR, nstates, nchareR, 
				     Rstates_per_pe, config.useCuboidMapRS, &GSmaptable, config.nchareG);
#endif
    }
    newtime=CmiWallTimer();
    CkPrintf("RSMap created in %g\n", newtime-Timer);
    Timer=newtime;
  } else {
    RSMapTable RStable= RSMapTable(&RSImaptable[0], &RSImaptable[numInst], availGlobR, 
			    nstates, nchareR, Rstates_per_pe, config.useCuboidMapRS, 
			    &GSImaptable[0], config.nchareG, numInst, x, y, z);
    CkPrintf("RSMap instance %d created in %g\n", numInst, newtime-Timer);
  }

  CProxy_RSMap rsMap= CProxy_RSMap::ckNew(thisInstance);
  //  CkArrayOptions realSpaceOpts(nstates,nchareR);
  CkArrayOptions realSpaceOpts(nstates,nchareR);
  realSpaceOpts.setMap(rsMap);
  int rforward=keeperRegister(string("RealSpaceForward"));
  int rbackward=keeperRegister(string("RealSpaceBackward"));

  UrealSpacePlaneProxy.push_back( CProxy_CP_State_RealSpacePlane::ckNew(1, 1, ngrida, ngridb, ngridc, rforward, rbackward, thisInstance, realSpaceOpts));
    UrealSpacePlaneProxy[thisInstance.proxyOffset].doneInserting();  

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareR;
    MapFile *mf = new MapFile("RSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RSmaptable);
#endif
    delete mf;
  }

  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareR;
    MapFile *mf = new MapFile("RSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RSmaptable);
#endif
    delete mf;

  }

 //--------------------------------------------------------------------------------
 // state r-particleplane

  availGlobR->reset();

  //--------------------------------------------------------------------------------
  // Do a fancy dance to determine placement of structure factors
    CkPrintf("Making SF non-EES\n");
    if(!ees_nonlocal_on)
      {
	CkPrintf("Making SF non-EES\n");
	int *nsend       = new int[nchareG];
	int **listpe     = new int * [nchareG];
	int numproc      = config.numPes;
	int *gspace_proc = new int [numproc];
   
	for(int i =0;i<numproc;i++){gspace_proc[i]=0;}
	for(int j=0;j<nchareG;j++){   
	  listpe[j]= new int[nstates];
	  nsend[j]=0;
	  for(int i=0;i<nstates;i++){
	    listpe[j][i] = gsprocNum(sim, i, j, thisInstance.getPO());
	    gspace_proc[listpe[j][i]]+=1;
	  }//endfor
	  lst_sort_clean(nstates, &nsend[j], listpe[j]);
	}//endfor
    
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
        PRINT_LINE_STAR;
	CkPrintf("Structure factor chareG dests\n");
	CkPrintf("Number of g-space chares : %d\n",nchareG);    
        PRINT_LINE_DASH;
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
        PRINT_LINE_DASH;
	CkPrintf("SFSends min %d max %d avg %g\n",minsend,maxsend,avgsend/(double)chareG_use);
        PRINT_LINE_STAR;

	//--------------------------------------------------------------------------------
	// Insert the objects into the StructureFactor array
	
	int dupmax=maxsend;  // there is no point in ever having more than that
	if(config.numSfDups<dupmax){dupmax=config.numSfDups;}
	config.numSfDups = dupmax;
	int numSfDups    = dupmax;
	CkPrintf("real numSfdups is %d based on maxsend of %d\n",numSfDups, maxsend);
	CkVec  <int>  peUsedBySF;
	for (int dup=0; dup<dupmax; dup++){
	  for (int x = 0; x < nchareG; x += 1){
	    int num_dup, istart, iend;
	    get_grp_params( nsend[x],  numSfDups,  dup, x ,&num_dup,  &istart, &iend);
	    int pe_ind=istart;
	    if(x%2==0)
	      pe_ind=iend;
	    for (int AtmGrp=0; AtmGrp<numSfGrps; AtmGrp++){
	      UsfCompProxy[thisInstance.proxyOffset](AtmGrp, x, dup).insert(numSfGrps,numSfDups, 
						 natm_nl_grp_max,  num_dup, &(listpe[x][istart]),
						  thisInstance,atmGrpMap(istart, num_dup, nsend[x], listpe[x],AtmGrp,dup,x));
	      peUsedBySF.push_back(atmGrpMap(istart, num_dup, nsend[x], listpe[x],   AtmGrp, dup,x));	      

	      pe_ind++;
	      if(pe_ind>nsend[x]){ pe_ind=0;}
	    }//endfor : AtmGrp
	  }//endfor : chareG
	}//endfor : Dups
	UsfCompProxy[thisInstance.proxyOffset].doneInserting();
	UpeUsedBySF.push_back(peUsedBySF);
	for(int j=0;j<nchareG;j++){delete [] listpe[j];}
	delete [] listpe;
	delete [] nsend;
      }

  //============================================================================
  // Set some com strategy of Sameer
  if(firstInstance) {

    if(config.useCommlib) {
      // TODO: do we need to do this for each Uber instance?
      // technically they are operating on mutually exclusive process
      // sets, so modulo insane commlib global state, it should be
      // safe to reuse these across ubers
      
      CkPrintf("Making State streaming strats\n");
      //mstrat->enableShortArrayMessagePacking();
      //rspaceState to gspaceState : gspaceState to rspaceState 

      if(config.useGssInsRealP)
	{
	  StreamingStrategy *gmstrat = new StreamingStrategy(config.gStreamPeriod,
							     config.gBucketSize);
	  gssInstance= ComlibRegister(gmstrat);    
	}
      //mstrat->enableShortArrayMessagePacking();
      //rPPState to gPPState : gPPState to rPPState 
      if(config.useMssInsGP){
        StreamingStrategy *mstrat = new StreamingStrategy(config.rStreamPeriod,
                                                          config.rBucketSize);

        mssInstance= ComlibRegister(mstrat);    
      }
      if (config.useMssInsGPP){
	StreamingStrategy *rpmstrat = new StreamingStrategy(config.rStreamPeriod,
							    config.rBucketSize);

	mssPInstance= ComlibRegister(rpmstrat);    
      }
      if (config.useGssInsRealPP){
	StreamingStrategy *gpmstrat = new StreamingStrategy(config.gStreamPeriod,
							    config.gBucketSize);
	gssPInstance= ComlibRegister(gpmstrat);    
      }
    }//endif
  }
//============================================================================

  printf("\n");
  
  // set firstInstance to false now
  firstInstance = false;

  PRINT_LINE_DASH;
  PRINTF("Completed G-space/R-space state chare array build\n");
  PRINT_LINE_STAR;printf("\n");

//----------------------------------------------------------------------------
    }//end routine
//============================================================================

//============================================================================
// Creating arrays CP_StateRealParticlePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void init_eesNL_chares(int natm_nl,int natm_nl_grp_max,
                       int doublePack, PeList *exclusion, CPcharmParaInfo *sim,
		       UberCollection thisInstance)
//============================================================================
   { //begin routine 
//============================================================================
/*
 * Set up the map variables used to control the location of the
 * 2d-array elements over processors   
*/
//============================================================================
// Useful Local variables


  int ngridaNl        = sim->ngrid_nloc_a;
  int ngridbNl        = sim->ngrid_nloc_b;
  int ngridcNl        = sim->ngrid_nloc_c;
  int ees_nonlocal_on = sim->ees_nloc_on;

  int nchareG         = sim->nchareG;
  int nchareRPP       = ngridcNl;

  int numIterNL       = sim->nlIters;
  int zmatSizeMax     = sim->nmem_zmat_max;

  PeList *nlexcludePes;
  if(config.useRhoExclusionMap)
    nlexcludePes=exclusion;
  else if(config.useReductionExclusionMap)
    nlexcludePes= new PeList(UpeUsedByNLZ[thisInstance.proxyOffset]);   
  else
    nlexcludePes= new PeList(0);
  if(config.excludePE0 && ! config.loadMapFiles)
      nlexcludePes->appendOne(0);
  int Rstates_per_pe  = config.Rstates_per_pe;
  availGlobG->reset();
  double newtime=CmiWallTimer();
#ifdef USE_INT_MAP
  RPPImaptable[thisInstance.getPO()].buildMap(nstates, nchareRPP);
#endif

  int success = 0;
  if(config.loadMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareRPP;
    MapFile *mf = new MapFile("RPPMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    success = mf->loadMap("RPPMap", &RPPImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("RPPMap", &RPPmaptable);
#endif
    delete mf;
  }

  if(success == 0) {
#ifdef USE_INT_MAP
    RPPMapTable RPPtable= RPPMapTable(&RPPImaptable[thisInstance.getPO()], availGlobG, nlexcludePes, 
				    nstates,  nchareRPP, Rstates_per_pe,
				    boxSize, config.useCuboidMap, 
				    config.nchareG, &GSImaptable[thisInstance.getPO()]);
#else
    RPPMapTable RPPtable= RPPMapTable(&RPPmaptable, availGlobG, nlexcludePes, 
				    nstates,  nchareRPP, Rstates_per_pe,
				    boxSize, config.useCuboidMap, 
				    config.nchareG, &GSmaptable);
#endif
  }
  CProxy_RPPMap rspMap= CProxy_RPPMap::ckNew(thisInstance);
  newtime=CmiWallTimer();
  CmiNetworkProgressAfter(0);
  CkPrintf("RPPMap created in %g\n",newtime-Timer);

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareRPP;
    MapFile *mf = new MapFile("RPPMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RPPImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RPPmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nstates; size[1] = nchareRPP;
    MapFile *mf = new MapFile("RPPMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RPPImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RPPmaptable);
#endif
    delete mf;
  }

  if(config.useRhoExclusionMap)
    {
      nlexcludePes=NULL;
    }
    else
    {
	if(nlexcludePes!=NULL)
	    delete nlexcludePes;
    }
  if(config.excludePE0 &&  !config.loadMapFiles)
  {
      if(nlexcludePes==NULL)
	  nlexcludePes=new PeList(0);
      nlexcludePes->appendOne(0);
  }
  Timer=newtime;
  CkArrayOptions pRealSpaceOpts(nstates,ngridcNl);
  pRealSpaceOpts.setMap(rspMap);
  UrealParticlePlaneProxy.push_back(CProxy_CP_State_RealParticlePlane::ckNew(
                                ngridaNl,ngridbNl,ngridcNl,
                                numIterNL,zmatSizeMax,Rstates_per_pe,
		                nchareG,ees_nonlocal_on, thisInstance,
				pRealSpaceOpts));
  UrealParticlePlaneProxy[thisInstance.proxyOffset].doneInserting();
  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed RealParticle chare array build\n");
  PRINT_LINE_STAR;printf("\n");

}

//============================================================================
// Creating arrays CP_Rho_GSpacePlane, CP_Rho_GSpacePlaneHelper 
// and CP_Rho_RealSpacePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_rho_chares(CPcharmParaInfo *sim, UberCollection thisInstance)
//============================================================================
{//begin routine
  //============================================================================
  /*
   * create the array for real-space densities (two-dimensional chare array)
   */    
  //============================================================================
  //  Chare array sizes and offsets 

  int ngrid_eext_a    = sim->ngrid_eext_a;
  int ngrid_eext_b    = sim->ngrid_eext_b;
  int ngrid_eext_c    = sim->ngrid_eext_c;
  int ees_eext_on     = sim->ees_eext_on;
  int ees_nonlocal_on = sim->ees_nloc_on;
  int natmTyp         = sim->natm_typ;
  int nchareRhoG      = sim->nchareRhoG;
  int nchareRhoR      = sim->sizeZ;
  int rhoGHelpers     = config.rhoGHelpers;
  int nchareHartAtmT  = config.nchareHartAtmT;
  int nchareRhoGHart  = rhoGHelpers*nchareRhoG;
  int nchareRhoRHart  = ngrid_eext_c;


  //============================================================================
  // Output to the screen

  PRINT_LINE_STAR;
  CkPrintf("Building RhoR, RhoG, RhoGHartExt, RhoRHartExt Chares %d %d %d %d natmtyp:%d\n",
	   nchareRhoR,nchareRhoG,nchareRhoGHart,nchareRhoRHart,natmTyp);
  PRINT_LINE_DASH;printf("\n");

  //============================================================================
  // Nuke some procs from the list : reset, nuke, reset if you run out

  
  PeList *RhoAvail=NULL;
  if(!config.loadMapFiles)
  {
      availGlobR->reset();
      RhoAvail=new PeList(*availGlobR);
  }
  //------------------------------------------------------------------------
  // subtract processors used by other nonscaling chares (non local reduceZ)
  excludePes= new PeList(0);   
  if(config.excludePE0 && !config.loadMapFiles)
      excludePes->appendOne(0);

  if(config.useReductionExclusionMap && !config.loadMapFiles)
    {
      if( nchareRhoR*config.rhoRsubplanes+UpeUsedByNLZ[thisInstance.proxyOffset].size() <
	  RhoAvail->count()){

	CkPrintf("subtracting %d NLZ nodes from %d for RhoR Map\n",
		 UpeUsedByNLZ[thisInstance.proxyOffset].size(),RhoAvail->count());
	//       nlz.dump();
	*RhoAvail-*excludePes; //unary minus
	RhoAvail->reindex();
	CkPrintf("Leaving %d for RhoR Map\n",RhoAvail->count());
      }//endif

      //------------------------------------------------------------------------
      // subtract processors used by other nonscaling chares

      if(ees_nonlocal_on==0){
	if( nchareRhoR*config.rhoRsubplanes+UpeUsedBySF[thisInstance.proxyOffset].size()<RhoAvail->count()){
	  CkPrintf("subtracting %d SF nodes from %d for RhoR Map\n",
		   UpeUsedBySF[thisInstance.proxyOffset].size(),RhoAvail->count());
	  PeList sf(UpeUsedBySF[thisInstance.proxyOffset]);
	  *RhoAvail-sf;
	  RhoAvail->reindex();
	}//endif
      }//endif
    }
  if(!config.loadMapFiles && RhoAvail->count()>2 ) { RhoAvail->reindex(); }

  //============================================================================
  // Maps and options
  //CkPrintf("RhoR map for %d x %d=%d chares, using %d procs\n",nchareRhoR, config.rhoRsubplanes, nchareRhoR*config.rhoRsubplanes, RhoAvail->count());

  //---------------------------------------------------------------------------
  // rho RS 
#ifdef USE_INT_MAP
  RhoRSImaptable[thisInstance.getPO()].buildMap(nchareRhoR, config.rhoRsubplanes);
#endif

  int success = 0;
  if(config.loadMapFiles) {
    int size[2];
    size[0] = nchareRhoR; size[1] = config.rhoRsubplanes;
    MapFile *mf = new MapFile("RhoRSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    success = mf->loadMap("RhoRSMap", &RhoRSImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("RhoRSMap", &RhoRSmaptable);
#endif
    delete mf;
  }

  if(success == 0) {
#ifdef USE_INT_MAP
    RhoRSMapTable RhoRStable(&RhoRSImaptable[thisInstance.getPO()], RhoAvail, nchareRhoR, config.rhoRsubplanes, config.nstates, config.useCentroidMapRho, &RSImaptable[thisInstance.getPO()], excludePes);
#else
    RhoRSMapTable RhoRStable(&RhoRSmaptable, RhoAvail, nchareRhoR,  config.rhoRsubplanes, config.nstates, config.useCentroidMapRho, &RSmaptable, excludePes);
#endif
  }

  CProxy_RhoRSMap rhorsMap = CProxy_RhoRSMap::ckNew(thisInstance);
  CkArrayOptions rhorsOpts(nchareRhoR, config.rhoRsubplanes);
  //CkArrayOptions rhorsOpts;
  rhorsOpts.setMap(rhorsMap);


  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nchareRhoR; size[1] = config.rhoRsubplanes;
    MapFile *mf = new MapFile("RhoRSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RhoRSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RhoRSmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nchareRhoR; size[1] = config.rhoRsubplanes;
    MapFile *mf = new MapFile("RhoRSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RhoRSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RhoRSmaptable);
#endif
    delete mf;
  }
  CmiNetworkProgressAfter(0);
  //---------------------------------------------------------------------------
  // rho GS 
  // if there aren't enough free procs refresh the RhoAvail list;
  if(!config.loadMapFiles && nchareRhoG>RhoAvail->count()) 
    {
      CkPrintf("refreshing avail list count %d less than rhog %d\n",RhoAvail->count(), nchareRhoG);
      RhoAvail->reset();
    }
#ifdef USE_INT_MAP
  RhoGSImaptable[thisInstance.getPO()].buildMap(nchareRhoG, 1);
#endif

  success = 0;
  if(config.loadMapFiles) {
    int size[2];
    size[0] = nchareRhoG; size[1] = 1;
    MapFile *mf = new MapFile("RhoGSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    success = mf->loadMap("RhoGSMap", &RhoGSImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("RhoGSMap", &RhoGSmaptable);
#endif
    delete mf;
  }

  if(success == 0) {
#ifdef USE_INT_MAP
    RhoGSMapTable RhoGStable(&RhoGSImaptable[thisInstance.getPO()], RhoAvail,nchareRhoG, config.useCentroidMapRho, &RhoRSImaptable[thisInstance.getPO()], excludePes);
#else
    RhoGSMapTable RhoGStable(&RhoGSmaptable, RhoAvail,nchareRhoG, config.useCentroidMapRho, &RhoRSmaptable, excludePes);
#endif
  }

  CProxy_RhoGSMap rhogsMap = CProxy_RhoGSMap::ckNew(thisInstance);
  CkArrayOptions rhogsOpts(nchareRhoG,1);
  //CkArrayOptions rhogsOpts;
  rhogsOpts.setMap(rhogsMap);

  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nchareRhoG; size[1] = 1;
    MapFile *mf = new MapFile("RhoGSMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RhoGSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RhoGSmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nchareRhoG; size[1] = 1;
    MapFile *mf = new MapFile("RhoGSMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RhoGSImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RhoGSmaptable);
#endif
    delete mf;
  }


  //---------------------------------------------------------------------------
  // rho RHart 
  // if there aren't enough free procs refresh the avail list;
  if(!config.loadMapFiles && nchareRhoRHart*nchareHartAtmT > RhoAvail->count())
    RhoAvail->reset();
  CkArrayOptions rhorhartOpts(nchareRhoRHart, config.rhoRsubplanes, nchareHartAtmT);
  //CkArrayOptions rhorhartOpts;
    
  if(ees_eext_on) {
#ifdef USE_INT_MAP
    RhoRHartImaptable[thisInstance.getPO()].buildMap(nchareRhoRHart, config.rhoRsubplanes, nchareHartAtmT);
#endif

    success = 0;
    if(config.loadMapFiles) {
      int size[3];
      size[0] = nchareRhoRHart; size[1] = config.rhoRsubplanes;
      size[2] = nchareHartAtmT;
      MapFile *mf = new MapFile("RhoRHartMap", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      success = mf->loadMap("RhoRHartMap", &RhoRHartImaptable[thisInstance.getPO()]);
#else
      success = mf->loadMap("RhoRHartMap", &RhoRHartmaptable);
#endif
      delete mf;
    }
    if(success == 0) {
#ifdef USE_INT_MAP
      RhoRHartMapTable RhoRHarttable(&RhoRHartImaptable[thisInstance.getPO()], RhoAvail, 
				     nchareRhoRHart, config.rhoRsubplanes, 
				     config.nchareHartAtmT, excludePes);
#else
      RhoRHartMapTable RhoRHarttable(&RhoRHartmaptable, RhoAvail,
				     nchareRhoRHart, config.rhoRsubplanes, 
				     config.nchareHartAtmT, excludePes);
#endif
    }

    CProxy_RhoRHartMap rhorHartMap = CProxy_RhoRHartMap::ckNew(thisInstance);
    rhorhartOpts.setMap(rhorHartMap);

    if(config.dumpMapFiles) {
      int size[3];
      size[0] = nchareRhoRHart; size[1] = config.rhoRsubplanes,
      size[2] = nchareHartAtmT;				  
      MapFile *mf = new MapFile("RhoRHartMap", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      mf->dumpMap(&RhoRHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
      mf->dumpMap(&RhoRHartmaptable);
#endif
      delete mf;
    }
    if(config.dumpMapCoordFiles) {
      int size[3];
      size[0] = nchareRhoRHart; size[1] = config.rhoRsubplanes,
      size[2] = nchareHartAtmT;				  
      MapFile *mf = new MapFile("RhoRHartMap_coord", 3, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
      mf->dumpMapCoords(&RhoRHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
      mf->dumpMapCoords(&RhoRHartmaptable);
#endif
      delete mf;
    }
  } //endif : ees_ext_on
  CkPrintf("RhoRHartMap built %d x %d x %d\n",nchareRhoRHart, config.rhoRsubplanes, config.nchareHartAtmT);
  //---------------------------------------------------------------------------
  // rho GHart 
  // if there aren't enough free procs refresh the avail list;
  if(!config.loadMapFiles && nchareRhoGHart>RhoAvail->count())
    RhoAvail->reset();
#ifdef USE_INT_MAP
  RhoGHartImaptable[thisInstance.getPO()].buildMap(nchareRhoGHart, nchareHartAtmT);
#endif

  success = 0;
  if(config.loadMapFiles) {
    int size[2];
    size[0] = nchareRhoGHart; size[1] = nchareHartAtmT;
    MapFile *mf = new MapFile("RhoGHartMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    success = mf->loadMap("RhoGHartMap", &RhoGHartImaptable[thisInstance.getPO()]);
#else
    success = mf->loadMap("RhoGHartMap", &RhoGHartmaptable);
#endif
    delete mf;
  }

  if(success == 0) {
#ifdef USE_INT_MAP
    MapType3 *RhoRHartImaptablep=NULL;
    if(ees_eext_on)
      RhoRHartImaptablep=&RhoRHartImaptable[thisInstance.getPO()];
    RhoGHartMapTable RhoGHarttable(&RhoGHartImaptable[thisInstance.getPO()], RhoAvail, 
				   nchareRhoGHart, config.nchareHartAtmT,
				   config.useCentroidMapRho, 
				   RhoRHartImaptablep, excludePes);
#else
    RhoGHartMapTable RhoGHarttable(&RhoGHartmaptable, RhoAvail, nchareRhoGHart,
				   config.useCentroidMapRho, &RhoRHartmaptable,
				   excludePes);
#endif
  }

  CProxy_RhoGHartMap rhogHartMap = CProxy_RhoGHartMap::ckNew(thisInstance);
  CkArrayOptions rhoghartOpts(nchareRhoGHart, nchareHartAtmT);
  //  CkArrayOptions rhoghartOpts;
  rhoghartOpts.setMap(rhogHartMap);
  CmiNetworkProgressAfter(0);
  if(config.dumpMapFiles) {
    int size[2];
    size[0] = nchareRhoGHart; size[1] =  nchareHartAtmT;
    MapFile *mf = new MapFile("RhoGHartMap", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMap(&RhoGHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMap(&RhoGHartmaptable);
#endif
    delete mf;
  }
  if(config.dumpMapCoordFiles) {
    int size[2];
    size[0] = nchareRhoGHart; size[1] =  nchareHartAtmT;
    MapFile *mf = new MapFile("RhoGHartMap_coord", 2, size, config.numPes, "TXYZ", 2, 1, 1, 1);
#ifdef USE_INT_MAP
    mf->dumpMapCoords(&RhoGHartImaptable[thisInstance.getPO()], thisInstance.getPO());
#else
    mf->dumpMapCoords(&RhoGHartmaptable);
#endif
    delete mf;
  }

  //============================================================================
  // Instantiate the chares

  bool dummy = true;

  //--------------------------------------------------------------------------
    
  // insert rhoreal
  int rhokeeper= keeperRegister(string("Density"));
  UrhoRealProxy.push_back(CProxy_CP_Rho_RealSpacePlane::ckNew(sizeX,dummy, 
						     ees_eext_on, ngrid_eext_c,
						     rhokeeper,
						     thisInstance,
						     rhorsOpts));
  /*    for (int i = 0; i < nchareRhoR; i++){
      for (int j = 0; j < config.rhoRsubplanes; j++){
	UrhoRealProxy[thisInstance.proxyOffset](i,j).insert(sizeX,dummy, ees_eext_on, ngrid_eext_c, rhokeeper);
      } //endfor
    } //endfor
    */
  UrhoRealProxy[thisInstance.proxyOffset].doneInserting();
  /// @todo: valgrind complains of a tiny memleak here. Check if callbacks get destroyed properly. 
  UrhoRealProxy[thisInstance.proxyOffset].ckSetReductionClient( new CkCallback(CkIndex_InstanceController::printEnergyEexc(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
  CmiNetworkProgressAfter(0);
  //--------------------------------------------------------------------------
  // insert rhog
  UrhoGProxy.push_back(CProxy_CP_Rho_GSpacePlane::ckNew(sizeX, 1, 
					       1, dummy, thisInstance,
					       rhogsOpts));
/*  for (int i = 0; i < nchareRhoG; i++){
    UrhoGProxy[thisInstance.proxyOffset](i,0).insert(sizeX, 1,1,dummy );
  }//endfor
*/
  UrhoGProxy[thisInstance.proxyOffset].doneInserting();
  CmiNetworkProgressAfter(0);
  //--------------------------------------------------------------------------
  // insert rhoghart
  UrhoGHartExtProxy.push_back(CProxy_CP_Rho_GHartExt::ckNew(ngrid_eext_a,ngrid_eext_b,
						   ngrid_eext_c,ees_eext_on,
						   natmTyp, thisInstance,
						   rhoghartOpts));
  /*
  for (int k = 0; k < nchareHartAtmT; k++){
    for (int i = 0; i < nchareRhoGHart; i++){
      UrhoGHartExtProxy[thisInstance.proxyOffset](i,k).insert(ngrid_eext_a,ngrid_eext_b,
				   ngrid_eext_c,ees_eext_on,natmTyp);
    }//endfor
  }//endfor
  */
  /// @todo: valgrind complains of a tiny memleak here. Check if callbacks get destroyed properly. 
  UrhoGHartExtProxy[thisInstance.proxyOffset].ckSetReductionClient(new CkCallback(CkIndex_InstanceController::printEnergyHart(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
  UrhoGHartExtProxy[thisInstance.proxyOffset].doneInserting();
  CmiNetworkProgressAfter(0);
  //--------------------------------------------------------------------------
  // insert rhoRhart
  if(ees_eext_on){
    UrhoRHartExtProxy.push_back(CProxy_CP_Rho_RHartExt::ckNew(ngrid_eext_a,ngrid_eext_b,
						     ngrid_eext_c,ees_eext_on,
						     natmTyp,thisInstance, 
						     rhorhartOpts));
    /*
    for (int k = 0; k < nchareHartAtmT; k++){
      for (int i = 0; i < nchareRhoRHart; i++){
	for (int j = 0; j < config.rhoRsubplanes; j++){
	  UrhoRHartExtProxy[thisInstance.proxyOffset](i,j,k).insert(ngrid_eext_a,ngrid_eext_b,ngrid_eext_c,
					 ees_eext_on,natmTyp);
	}//endfor
      }//endfor
    }//endfor
    */
    UrhoRHartExtProxy[thisInstance.proxyOffset].doneInserting();
  }//endif
  CmiNetworkProgressAfter(0);
  //===========================================================================
  // Output to the screen
  // need to add maps for these, for now just let em default
  // IF some condition which triggers QMMM
  int nchareRhoGLSP=1;
  int nchareRhoRealLSP=1;
  int nchareRhoRealLSPsubplanes=1;
  CkArrayOptions lspgspOpts(nchareRhoGLSP);
  CkArrayOptions lsprealOpts(nchareRhoRealLSP, nchareRhoRealLSPsubplanes);
  UlsRhoGProxy.push_back(CProxy_CP_LargeSP_RhoGSpacePlane::ckNew(thisInstance,lspgspOpts));
  UlsRhoRealProxy.push_back(CProxy_CP_LargeSP_RhoRealSpacePlane::ckNew(thisInstance,lsprealOpts));
  printf("\n");
  PRINT_LINE_DASH;
  PRINTF("Completed G-space/R-space Rho chare array build\n");
  PRINT_LINE_STAR;printf("\n");
  if(RhoAvail!=NULL)
    delete RhoAvail;

  //===========================================================================
}//end routine
//============================================================================



//============================================================================
// Get the atoms and the parainfo
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void control_physics_to_driver(UberCollection thisInstance){
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
    UatomsGrpProxy.push_back( CProxy_AtomsGrp::ckNew(natm,natm_nl,len_nhc,iextended_on,
                                           cp_min_opt,cp_wave_opt,isokin_opt,
                                           kT,atoms,atomsNHC,thisInstance));
    delete [] atoms;
    delete [] atomsNHC;
    delete PhysicsAtom;

//=====================================================================
// Make a group for the energies

    UegroupProxy.push_back(CProxy_EnergyGroup::ckNew(thisInstance)); 

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
static CkReductionMsg *complexSum(int nMsg, CkReductionMsg **msgs){
//============================================================================
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
//============================================================================
 {// begin routine
//============================================================================

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
	 //         CkExit();
       }
  
   (*n_ret)     = n;
   (*istrt_ret) = istrt;
   (*iend_ret)  = iend;

//---------------------------------------------------------------------------
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int atmGrpMap(int istart, int nsend, int listsize, int *listpe, int AtmGrp, 
              int dup, int planeIndex)
//============================================================================
  {// begin routine
//============================================================================

  int numSfDups=config.numSfDups;

  if(listsize <= numSfDups){ //dup list not unique
      return listpe[(dup % listsize)];
  }//endif

  if(nsend <= config.numSfGrps){  //atom group list not unique
      return listpe[(istart + (AtmGrp % nsend))];
  }//endif

  // nsend must be > numSfGrps and its list is unique
  return listpe[(istart + ((AtmGrp + config.numSfGrps * planeIndex)%nsend))];

//---------------------------------------------------------------------------
  }// end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int gsprocNum(CPcharmParaInfo *sim, int state, int plane, int numInst) {
  int proc;
#ifdef USE_INT_MAP
  proc = GSImaptable[numInst].get(state, plane);
#else
  proc = GSmaptable.get(intdual(state, plane));
#endif
  if(config.fakeTorus)
    proc=proc%CkNumPes();
  return(proc);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void mapOutput()
//============================================================================
 {//begin routine
//============================================================================

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// return the cuboid x,y,z of a subpartition exactly matching that volume
//============================================================================
bool findCuboid(int &x, int &y, int &z, int &order, int maxX, int maxY, int maxZ, int maxT, int volume, int vn){
//============================================================================
  int maxD=maxX;
  int minD=maxX;
  int minX=maxX;
  int minY=maxY;
  int minZ=maxZ;
  //  vn=0;
  if(config.forceMappingAxis>-1)
    {
      order=config.forceMappingAxis;
      if (config.forceMappingAxis==0) 
	minD=minX;
      if (config.forceMappingAxis==1) 
	minD=maxY;
      if (config.forceMappingAxis==2) 
	minD=maxZ;
    }
  else
    { 
      order=0;
      maxD = (maxY>maxD) ? maxY : maxD;
      maxD = (maxZ>maxD) ? maxZ : maxD;
      if(vn)
	{  // using Y as the prism axis seems to suck
	  maxD = (maxY>maxD) ? maxY : maxD;
	  maxD = (maxZ>maxD) ? maxZ : maxD;
	  minD = (maxY<minD) ? maxY : minD;
	  minD = (maxZ<maxD) ? maxZ : minD;

	}
      else
	{
	  maxD = (maxY>maxD) ? maxY : maxD;
	  maxD = (maxZ>maxD) ? maxZ : maxD;
	  minD = (maxY<minD) ? maxY : minD;
	  minD = (maxZ<minD) ? maxZ : minD;

	}
    }
  // CkPrintf("minD %d maxD %d\n",minD, maxD);
  if(config.useCuboidMapRS)
    {
      CkPrintf("Using long prisms for useCuboidMapRS\n");
    }

  // We were reducing the volume by maxT and then finding the dimensions of the
  // box in terms of the no. of nodes and not processors

  // That worked ok on BG/L, but BG/P maxT of 4 distorts the long axis logic.
  // in TXYZ you have a 32x8x16 for 1 Rack or 32x8x32 for 2 Rack
  // in XYZT you have 8x8x64 for 1 Rack or 8x8x128 for 2 Rack

  // The latter case makes clear the bizarre and unreal distortion created by
  // pretending that maxT is just a multiplier along one dimension when
  // considering actual box size.  What was merely odd on BG/L is now
  // seriously peculiar.  Also default TXYZ mappings and default Order 0
  // mappings to use the X as long axis for prisms fail to work right because
  // we aren't really spanning the X*maxT dimension with these long prisms.
  
  int redVol = volume / maxT;
  double cubert= cbrt((double) redVol);
  int cubetrunc= (int) cubert;
  x=y=z=cubetrunc;
  if(cubetrunc>minD)
    cubetrunc=minD;
  if(cubetrunc>maxY)
    cubetrunc=maxY;
  if(cubetrunc>maxZ)
    cubetrunc=maxZ;
  if(redVol==x*y*z && !config.useCuboidMapRS)
    return true;
  bool switchSet=false;
  CkAssert(redVol>0);
  switch (redVol) // for the common values we just pick cuboids we like
    {
    case 1: 
      x=1; y=1; z=1; switchSet=true; break;
    case 2: 
      x=2; y=1; z=1; switchSet=true; break;
    case 3:
      x=3; y=1; z=1; switchSet=true; break;
    case 4:
      x=2; y=2; z=1; switchSet=true; break;
    case 5:
      x=5; y=1; z=1; switchSet=true; break;
    case 6:
      x=3; y=2; z=1; switchSet=true; break;
    case 7:
      x=7; y=1; z=1; switchSet=true; break;
    case 8:
      if(config.useCuboidMapRS)
      {
	  if(minD>=2)

	  { x=2; y=2; z=2; switchSet=true; break;}
/* no evidence that these other schemes help
	  if(minD==4)
	  { x=4; y=2; z=1; switchSet=true; break;}
	  if(minD>=8)
	  { x=8; y=1; z=1; switchSet=true; break;}
*/
	}
    case 9:
      x=3; y=3; z=1; switchSet=true; break;
    case 10:
      x=5; y=2; z=1; switchSet=true; break;
    case 12:
      x=2; y=3; z=2; switchSet=true; break;
    case 14:
      x=7; y=2; z=1; switchSet=true; break;
    case 15:
      x=5; y=3; z=1; switchSet=true; break;
    case 16:
      if(config.useCuboidMapRS)
	{
	  if(minD>=8)
	    { x=8; y=2; z=1; switchSet=true; break;}
	}
      x=4; y=2; z=2; switchSet=true; break;
    case 18:
      x=3; y=3; z=2; switchSet=true; break;
    case 20:
      x=5; y=2; z=2; switchSet=true; break;
    case 21:
      x=7; y=3; z=1; switchSet=true; break;
    case 24:
      x=4; y=3; z=2; switchSet=true; break;
    case 25:
      x=5; y=5; z=1; switchSet=true; break;
    case 27:
      x=3; y=3; z=3; switchSet=true; break;
    case 28:
      x=7; y=2; z=2; switchSet=true; break;
    case 30:
      x=5; y=2; z=2; switchSet=true; break;
    case 32:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=2; z=2; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=2; z=1; switchSet=true; break;}
	  if(minD>=32)
	    { x=32; y=1; z=1; switchSet=true; break;}

	}
      x=4; y=2; z=4; switchSet=true; break;
    case 35:
      x=7; y=5; z=1; switchSet=true; break;
    case 36:
      x=4; y=3; z=3; switchSet=true; break;
    case 40:
      x=5; y=4; z=2; switchSet=true; break;
    case 42:
      x=7; y=3; z=2; switchSet=true; break;
    case 43:
      x=7; y=3; z=2; switchSet=true; break;
    case 45:
      x=5; y=3; y=3; switchSet=true; break;
    case 48:
      x=4; y=3; z=4; switchSet=true; break;
    case 50:
      x=5; y=5; z=2; switchSet=true; break;
    case 54:
      x=6; y=3; z=3; switchSet=true; break;
    case 56:
      x=7; y=4; z=2; switchSet=true; break;
    case 60:
      x=5; y= 4; z=3; switchSet=true; break;
    case 64:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=4; z=2; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=2; z=2; switchSet=true; break;}
	  if(minD==32)
	    { x=32; y=2; z=1; switchSet=true; break;}
	  if(minD>=64)
	    { x=64; y=1; z=1; switchSet=true; break;}
	}
      x=4; y=4; z=4; switchSet=true; break;
    case 128:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=4; z=4; switchSet=true; break;}
	  if(minD==16)
	    {x=16; y=4; z=2; switchSet=true; break;}
	  if(minD==32)
	    {  x=32; y=2; z=2; switchSet=true; break;	}
	  if(minD==64)
	    {  x=64; y=2; z=1; switchSet=true; break;	}
	  if(minD>=128)
	    {  x=128; y=1; z=1; switchSet=true; break;	}
	}
      x=8; y=4; z=4; switchSet=true; break;
    case 256:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=8; z=4; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=4; z=4; switchSet=true; break;}
	  if(minD==32)
	    { x=32; y=4; z=2; switchSet=true; break;}
	  if(minD==64)
	    { x=64; y=2; z=2; switchSet=true; break;}
	  if(minD>=128)
	    { x=64; y=2; z=1; switchSet=true; break;}
	}
      x=8; y=8; z=4; switchSet=true; break;
    case 512:
      if(config.useCuboidMapRS)
	{
	  if(minD==8)
	    { x=8; y=8; z=8; switchSet=true; break;}
	  if(minD==16)
	    { x=16; y=4; z=8; switchSet=true; break;}
	  if(minD==32)
	    { x=32; y=4; z=4; switchSet=true; break;}
	  if(minD==64)
	    { x=64; y=4; z=2; switchSet=true; break;}
	  if(minD==128)
	    { x=128; y=2; z=2; switchSet=true; break;}
	  if(minD>=256)
	    { x=256; y=2; z=1; switchSet=true; break;}
	}
      x=8; y=8; z=8; switchSet=true; break;

    default:
      break;
    }
  if(switchSet && config.forceMappingAxis>-1)
    {
      
      switch(config.forceMappingAxis)
	{
	case 0:
	  return true; 	  //no change
	case 1:
	  {
	  order=1; //YXZ
	  int swap=x;
	  x=y;
	  y=swap;
	  return true;
	  }
	case 2:
	  {
	  order=2; //ZXY
	  int swap=x;
	  x=z;
	  z=swap;
	  return true;
	  }
	}
    }
  else if (switchSet)
    {
      // now correct the x,y,z to put long prism axis along the
      // smallest torus dimension which will fit.
      if(x==maxX)
	return true;
      if(x==maxY)
	{ // change to Y
	  order=1; //YXZ
	  int swap=x;
	  x=y;
	  y=swap;
	  return true;
	}
      if(x==maxZ)
	{ // change to Z
	  order=2; //ZXY
	  int swap=x;
	  x=z;
	  z=swap;
	  return true;
	}
      // if we're here then we don't have a spanning prism
      // just pick the smallest which will fit.
      if(x<maxX)
	return true;
      if(x<maxY)
	{ // change to Y
	  order=1; //YXZ
	  int swap=x;
	  x=y;
	  y=swap;
	  return true;
	}
      if(x<maxZ)
	{ // change to Z
	  order=2; //ZXY
	  int swap=x;
	  x=z;
	  z=swap;
	  return true;
	}
    }
  // its something weird so try a best fit
  int start=cubetrunc-1;
  if(config.useCuboidMapRS)
    x = (redVol>=maxX) ? maxX : cubetrunc;
  else
    x=cubetrunc;
  for(; x<=maxX;x++)
    {
      for(y=start; y<=maxY;y++)
	{
	  for(z=start; z<=maxZ;z++)
	    {
	      if(redVol==x*y*z)
		return true;
	    }
	}
    }
  return false;


}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CkReduction::reducerType sumFastDoubleType;
void registersumFastDouble(void){ 
  sumFastDoubleType=CkReduction::addReducer(sumFastDouble);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// sum together matrices of doubles
// possibly faster than sum_double due to minimizing copies
// and using sameer's fastadd
//============================================================================
inline CkReductionMsg *sumFastDouble(int nMsg, CkReductionMsg **msgs){

  int size0=msgs[0]->getSize();
  int size=size0/sizeof(double);

  double *inmatrix;
  //  int progcount=0;

  double *ret=(double *)msgs[0]->getData();

  for(int i=1; i<nMsg;i++)
    {
#ifdef CMK_BLUEGENEL
      fastAdd(ret, (double *) msgs[i]->getData(), size);
#else
      inmatrix=(double *) msgs[i]->getData();
	for(int d=0;d<size;d++)
	  ret[d]+=inmatrix[d];
#endif
    }
  return CkReductionMsg::buildNew(size0,ret);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void create_Rho_fft_numbers(int nchareR,int nchareRHart,int rhoRsubplanes,
                            int nplane, int sizeY, int ngridbHart,
                            int *numRXRho,int *numRYRho,int *numRXEext,int *numRYEext,
                            int *numSubGx)
//============================================================================
    {//begin routine
//============================================================================
// Simple case  : everyone is the same size

  if(rhoRsubplanes==1){

    for(int i=0;i<nchareR;i++){
       numRXRho[i]  = sizeY;
       numRYRho[i]  = nplane;
    }//endfor

    for(int i=0;i<nchareRHart;i++){
       numRXEext[i] = ngridbHart;
       numRYEext[i]   = nplane;
    }//endfor

  }//endif

//============================================================================
// Subplane case : different sizes with subplanes label

  if(rhoRsubplanes>1){
    int div,rem;

    //---------------------------------------
    // how many x FFTs for Rho
    div  = (sizeY / rhoRsubplanes); 
    rem  = (sizeY % rhoRsubplanes);
    for(int j=0;j<rhoRsubplanes;j++){
      int mySizeY = (j < rem ? div+1 : div);
      for(int i=0;i<nchareR;i++){
        int ind       = j*nchareR + i;
        numRXRho[ind] = mySizeY;
     }//endfor
    }//endfor

    //---------------------------------------
    // how many y FFTs for Rho
    for(int j=0;j<rhoRsubplanes;j++){
      int myNplane = numSubGx[j];
      for(int i=0;i<nchareR;i++){
        int ind       = j*nchareR + i;
        numRYRho[ind] = myNplane;
      }//endfor
    }//endfor

    //---------------------------------------
    // how many x FFTs for HartEExt EES
    div  = (ngridbHart / rhoRsubplanes); 
    rem  = (ngridbHart % rhoRsubplanes);
    for(int j=0;j<rhoRsubplanes;j++){
      int myNgridb   = (j < rem ? div+1 : div);
      for(int i=0;i<nchareRHart;i++){
        int ind        = j*nchareRHart + i;
        numRXEext[ind] = myNgridb;
     }//endfor
    }//endfor
    
    //---------------------------------------
    // how many y FFTs for HartEExt EES
    for(int j=0;j<rhoRsubplanes;j++){
      int myNplane = numSubGx[j];
      for(int i=0;i<nchareRHart;i++){
        int ind        = j*nchareRHart + i;
        numRYEext[ind] = myNplane;
      }//endfor
    }//endfor

  }//endif

//---------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
#include "CPcharmParaInfo.def.h"
#include "cpaimd.def.h"

//============================================================================


