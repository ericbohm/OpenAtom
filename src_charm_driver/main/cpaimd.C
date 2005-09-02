//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//                         cpaimd-charm-driver
//     Software developed by the Parallel Programing Laboratory, UIUC.
//     in collaboration with IBM and NYU.
//    
//     This file contains cpaimd-charm-driver main. It creates and 
//     initializes all the arrays and libraries. 
//      
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

#include <math.h>
#include "charm++.h"
#include "ckarray.h"
#include "util.h"

#include "groups.h"
#include "cpaimd.h"
#include "ortho.h"
#include "matmul.h"
#include "sim_subroutines.h"
#include "StructFactorCache.h"
#include "StructureFactor.h"
#include "CP_State_Plane.h"
#include "MeshStreamingStrategy.h"

#include "../include/CPcharmParaInfo.h"
#include "../../src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "../../src_piny_physics_v1.0/include/charm_defs/Interface_ctrl.decl.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h"

#ifdef MACOSX_BLAS
#include   <Accelerate/Accelerate.h>
#endif

//============================================================================
// Defining all Charm++ readonly variables for PINY physics 

extern MDINTEGRATE  readonly_mdintegrate;
extern MDATOMS      readonly_mdatoms;
extern MDINTER      readonly_mdinter;
extern MDINTRA      readonly_mdintra;
extern GENERAL_DATA readonly_general_data;
extern CP           readonly_cp; 

#include "MultiRingMulticast.h"

//============================================================================
/*
 * Defining all the Charm++ readonly variables, which include proxies
 * to access the arrays and groups and the Communication Library
 * handles.
 */

Config config;
PairCalcID pairCalcID1;
PairCalcID pairCalcID2;
CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
CProxy_CP_State_ParticlePlane particlePlaneProxy;
CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
CProxy_CP_Rho_GSpacePlane rhoGProxy;
CProxy_CP_Rho_GSpacePlaneHelper rhoGHelperProxy;
CProxy_Ortho orthoProxy;
CProxy_matmul matmulProxy1;
CProxy_matmul matmulProxy2;
CProxy_matmul matmulProxy3;
CProxy_CPcharmParaInfoGrp scProxy;
CProxy_AtomsGrp atomsGrpProxy;
CProxy_EnergyGroup egroupProxy;
CProxy_FFTcache fftCacheProxy;
CProxy_StructFactCache sfCacheProxy;
CProxy_StructureFactor sfCompProxy;

int atom_integrate_done;  // not a real global : more like a group of size 1

int nstates;  // readonly globals
int sizeX;
int Ortho_UE_step2;
int Ortho_UE_step3;
int Ortho_UE_error;
bool Ortho_use_local_cb;


//============================================================================
// For using the multicast library :  Set some reduction clients

CkGroupID mCastGrpId; 
// For using the communication library
ComlibInstanceHandle commInstance;
ComlibInstanceHandle commRealInstance;
ComlibInstanceHandle mcastInstance;
ComlibInstanceHandle ssInstance;
ComlibInstanceHandle mssInstance;
ComlibInstanceHandle mcastInstancePP;
CkReduction::reducerType complexVectorAdderType;

#include "ReductionClients.h"


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/* dummy function for using the PairCalculator library*/
void myFunc(complex a, complex b) {}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 * The Main of CPAIMD. It calls all the init functions
 */
//============================================================================

main::main(CkArgMsg *m) {

//============================================================================
/* Check arguments : Tell people what we are doing*/

    if (m->argc < 3) {
      CkAbort("Usage: pgm config_file_charm config_file_piny");
    }//endif

    CkPrintf("\n================================================\n");
    CkPrintf("Starting Cpaimd-Charm-Driver Setup Phase \n");
    CkPrintf("---------------------------------------------------\n");
    CkPrintf("  Cpaimd-Charm-Driver running on %d processors. \n", CkNumPes());
    CkPrintf("  Reading Physics input from %s\n",m->argv[2]);
    CkPrintf("  Reading Driver  input from %s\n",m->argv[1]);
    CkPrintf("---------------------------------------------------\n\n");

//============================================================================    
/* Invoke PINY input class */

    CkCallback piny_callback (CkCallback::ignore);
    Interface_ctrl piny_interface (m->argv[2],piny_callback);

    CPcharmParaInfo *sim  = new CPcharmParaInfo();
    PhysicsParamTransfer::ParaInfoInit(sim);
    int ibinary_opt = sim->ibinary_opt;
    int natm_nl     = sim->natm_nl;

//============================================================================    
/* Invoke Ramkumar input class */

    CkPrintf("\n================================================\n");
    CkPrintf("Cpaimd-Charm-Driver input started \n");
    CkPrintf("------------------------------------------------\n");

    Config::readConfig(m->argv[1],config,sim->nstates,
                       sim->sizeX,sim->sizeY,sim->sizeZ,
                       sim->ntime,ibinary_opt,natm_nl);

    int numSfGrps    = config.numSfGrps;  // local copies are nice
    int doublePack   = config.doublePack;
    size2d sizeYZ    = size2d(sim->sizeY,sim->sizeZ);
    int gSpacePPC    = config.gSpacePPC;  
    int realSpacePPC = config.realSpacePPC;
    int rhoGPPC      = config.rhoGPPC;

    nstates          = config.nstates;    // globals on all procs
    sizeX            = sim->sizeX;

    config.print();

    CkPrintf("------------------------------------------------\n");
    CkPrintf("Cpaimd-Charm-Driver input completed \n");
    CkPrintf("================================================\n\n");

//============================================================================
// Set user trace events

     traceRegisterUserEvent("doRealFwFFT", doRealFwFFT_);
     traceRegisterUserEvent("doRealBwFFT", doRealBwFFT_);
     traceRegisterUserEvent("GspaceFwFFT", GspaceFwFFT_);
     traceRegisterUserEvent("GspaceBwFFT", GspaceBwFFT_);
     traceRegisterUserEvent("RhoRtoGxyFFT", RhoRtoGxzFFT_);
     traceRegisterUserEvent("RhoRtoGyFFT", RhoRtoGyFFT_);
     traceRegisterUserEvent("RhoDivRhoXFFT", RhoDivRhoXFFT_);
     traceRegisterUserEvent("RhoDivRhoYFFT", RhoDivRhoYFFT_);
     traceRegisterUserEvent("RhoDivRhoZFFT", RhoDivRhoZFFT_);
     traceRegisterUserEvent("VksofGFFT", VksofGFFT_);
     traceRegisterUserEvent("VksofRFFT", VksofRFFT_);
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

    scProxy = CProxy_CPcharmParaInfoGrp::ckNew(*sim);

//============================================================================    
// Create the multicast/reduction manager for array sections
// Create the parainfo group from sim
// Initialize chare arrays for real and g-space of states 

    mCastGrpId = CProxy_CkMulticastMgr::ckNew();

    init_planes(sizeYZ,natm_nl,natm_nl_grp_max,numSfGrps,doublePack,
                gSpacePPC,realSpacePPC,sim);

//============================================================================    
// Transfer parameters from physics to driver
//    read in atoms/states : create atoms group 
    
    control_physics_to_driver();

//============================================================================ 
// Some intense initialization of paircalculator : Comments anyone?
//   Can we hide all this in a function call?

  //-------------------------------------------------------------
  // Create mapping classes for Paircalcular
    int indexSize = 0;
    if(doublePack){
    	indexSize = config.low_x_size;
    }else{
    	indexSize = config.low_x_size + (sizeX-config.high_x_size-1);
    }

    int* indexZ = new int[indexSize];
    for(int i=0, count=0; i<sizeX; i++){
        if(i < config.low_x_size || ((i > config.high_x_size)&&(!doublePack))){
        indexZ[count] = i;
        count++;
      }
    }

    CProxy_SCalcMap scMap_sym = CProxy_SCalcMap::ckNew(config.nstates,
                   sizeX / gSpacePPC,config.sGrainSize,CmiTrue,sim->nplane_x, 
                   sim->lines_per_plane, sim->pts_per_plane) ;
    
    CProxy_SCalcMap scMap_asym = CProxy_SCalcMap::ckNew(config.nstates,
  	           sizeX / gSpacePPC,config.sGrainSize, CmiFalse,sim->nplane_x, 
                   sim->lines_per_plane, sim->pts_per_plane);
    
    CkGroupID scalc_sym_id = scMap_sym.ckGetGroupID();
    CkGroupID scalc_asym_id = scMap_asym.ckGetGroupID();


  //-------------------------------------------------------------
  // Register the PCs

    int gsp_ep =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_mySendMsg;
    if(config.gspacesum)
      gsp_ep =  CkIndex_CP_State_GSpacePlane::__idx_acceptNewPsi_partialResultMsg;



    createPairCalculator(true, nstates, config.sGrainSize, indexSize, indexZ, 0, myFunc, 0, myFunc, CkCallback(CkIndex_Ortho::start_calc(NULL), orthoProxy), &pairCalcID1, gsp_ep, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_sym_id, doublePack, config.conserveMemory,config.lbpaircalc, config.psipriority, mCastGrpId, config.gspacesum );

    CkArrayIndex2D myindex(0, 0);

    gsp_ep = CkIndex_CP_State_GSpacePlane::__idx_acceptLambda_mySendMsg;
    int myPack=0;

    bool asymm_gspacesum=false; // not supported yet

//    createPairCalculator(false, nstates,  config.sGrainSize, indexSize, indexZ, 0, myFunc, 0, myFunc,CkCallback(CkIndex_Ortho::acceptAllLambda(NULL), myindex, gSpacePlaneProxy.ckGetArrayID()), &pairCalcID2, gsp_ep, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_asym_id, myPack, config.conserveMemory,config.lbpaircalc, config.lambdapriority, mCastGrpId, asymm_gspacesum );
    createPairCalculator(false, nstates,  config.sGrainSize, indexSize, indexZ, 0, myFunc, 0, myFunc,CkCallback(CkIndex_CP_State_GSpacePlane::acceptAllLambda(NULL), myindex, gSpacePlaneProxy.ckGetArrayID()), &pairCalcID2, gsp_ep, gSpacePlaneProxy.ckGetArrayID(), 1, &scalc_asym_id, myPack, config.conserveMemory,config.lbpaircalc, config.lambdapriority, mCastGrpId, asymm_gspacesum );
    
  //-------------------------------------------------------------
  // Create stuff for ortho which PC invokes by section reduction

    orthoProxy = CProxy_Ortho::ckNew();
    CkArrayOptions opts(0);
    opts.bindTo(orthoProxy);
    matmulProxy1 = CProxy_matmul::ckNew(opts);
    matmulProxy2 = CProxy_matmul::ckNew(opts);
    matmulProxy3 = CProxy_matmul::ckNew(opts);
    int init_pe = 0;
    CkCallback *orthoReduction = new CkCallback(CkIndex_Ortho::collect_error(NULL), orthoProxy(0, 0));
    orthoProxy.ckSetReductionClient(orthoReduction);
    
    // extra triangle ortho elements are really a waste of our time
    // and resources, but we don't have a triangular solver for
    // inv_square, so we'll just make do.

    // They need to exist solely so that the inv_sq method can work.
    // So we need to copy their mirror elements data into them.
    // then when complete they need to know not to call finishpaircalc.
    // Because their redundant data has nowhere to go.

    // punch in an actual map for this, 3d cube for BG/L  spread them out to maximize bi-section band
    CProxySection_PairCalculator lambdaSectProxy;
    for (int s1 = 0; s1 < nstates; s1 += config.sGrainSize) {
      for (int s2 = s1; s2 < nstates; s2 += config.sGrainSize) {

	int indX = s1 / config.sGrainSize;
	int indY = s2 / config.sGrainSize;
	CkCallback cb = CkCallback(CkIndex_Ortho::start_calc(NULL), orthoProxy(indX, indY));
	CkCallback cbl = CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), orthoProxy(indX, indY));
	if(config.parlambda)
	  lambdaSectProxy=initOneRedSect(indexSize, indexZ, 1, &pairCalcID2, cbl, s1, s2, 0);	    
	CProxySection_PairCalculator sectProxy=initOneRedSect(indexSize, indexZ, 1, &pairCalcID1, cb, s1, s2, 0);	    

	orthoProxy(indX, indY).insert(config.sGrainSize, sectProxy, lambdaSectProxy, init_pe);
	int chunks = nstates / config.sGrainSize;
	matmulProxy1(indX, indY).insert(chunks, config.sGrainSize, CkCallback(CkIndex_Ortho::ready(), orthoProxy(indX, indY)), init_pe);
	matmulProxy2(indX, indY).insert(chunks, config.sGrainSize, CkCallback(CkIndex_Ortho::ready(), orthoProxy(indX, indY)), init_pe);
	matmulProxy3(indX, indY).insert(chunks, config.sGrainSize, CkCallback(CkIndex_Ortho::ready(), orthoProxy(indX, indY)), init_pe);
	if(s2>s1) // non diagonal 
	  { 
	    init_pe = (init_pe + 1) % CkNumPes();
	    CkCallback cbl = CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), orthoProxy(indY, indX));
	    if(config.parlambda)
	      lambdaSectProxy=initOneRedSect(indexSize, indexZ, 1, &pairCalcID2, cbl, s2, s1, 0);	    
	    orthoProxy(indY, indX).insert(config.sGrainSize, sectProxy, lambdaSectProxy, init_pe);
	    int chunks = nstates / config.sGrainSize;
	    matmulProxy1(indY, indX).insert(chunks, config.sGrainSize, CkCallback(CkIndex_Ortho::ready(), orthoProxy(indY, indX)), init_pe);
	    matmulProxy2(indY, indX).insert(chunks, config.sGrainSize, CkCallback(CkIndex_Ortho::ready(), orthoProxy(indY, indX)), init_pe);
	    matmulProxy3(indY, indX).insert(chunks, config.sGrainSize, CkCallback(CkIndex_Ortho::ready(), orthoProxy(indY, indX)), init_pe);
	  }
	init_pe = (init_pe + 1) % CkNumPes();
	
      }
    }
    orthoProxy.doneInserting();
    matmulProxy1.doneInserting();
    matmulProxy2.doneInserting();
    matmulProxy3.doneInserting();

    delete [] indexZ;

//============================================================================ 
// Initialize the density chare arrays
    
    init_rho(sizeYZ,gSpacePPC,realSpacePPC,rhoGPPC);

//============================================================================ 
// Sameer's new communication strategies  : function call anyone?
//                                          comments?

    int i = 0;
    //Initialize the communication library strategies.
    if (config.useCommlib) {        
        int nsrcelements = sizeYZ[1]/rhoGPPC;
        int ndestelements = sizeYZ[0]/realSpacePPC;

        CkArrayIndexMax *srcelements = new CkArrayIndexMax[sizeYZ[1]/rhoGPPC];
        for (i = 0; i < sizeYZ[1]/rhoGPPC; i++) {
            srcelements[i] = CkArrayIndex1D(i);
        }

        CkArrayIndexMax *destelements = new 
            CkArrayIndexMax[sizeYZ[0]/realSpacePPC];
        for(i = 0; i < sizeYZ[0]/realSpacePPC; i++) {
            destelements[i] = CkArrayIndex1D(i);
        }

        CharmStrategy *strat = new EachToManyMulticastStrategy
            (USE_MESH, rhoGProxy.ckGetArrayID(), rhoRealProxy.ckGetArrayID(), 
             nsrcelements, srcelements, ndestelements, destelements);
        
        srcelements = new CkArrayIndexMax[sizeYZ[1]/rhoGPPC];
        for (i = 0; i < sizeYZ[1]/rhoGPPC; i++) 
            srcelements[i] = CkArrayIndex1D(i);
        
        destelements = new CkArrayIndexMax[sizeYZ[0]/realSpacePPC];
        for(i = 0; i < sizeYZ[0]/realSpacePPC; i++)
            destelements[i] = CkArrayIndex1D(i);
        
        CharmStrategy *real_strat = new EachToManyMulticastStrategy
            (USE_MESH, rhoRealProxy.ckGetArrayID(), rhoGProxy.ckGetArrayID(),
             ndestelements, destelements, nsrcelements, srcelements);
        
        commInstance = CkGetComlibInstance();
        commRealInstance = CkGetComlibInstance();
        commInstance.setStrategy(strat);
        commRealInstance.setStrategy(real_strat);
    }

    if (config.useCommlibMulticast) {
        mcastInstance = CkGetComlibInstance();         
        mcastInstancePP = CkGetComlibInstance();         
        DirectMulticastStrategy *dstrat = new DirectMulticastStrategy
	    (realSpacePlaneProxy.ckGetArrayID(), 1);
        
        RingMulticastStrategy *rstrat = new RingMulticastStrategy
            (realSpacePlaneProxy.ckGetArrayID(), 1);
        
        MultiRingMulticast *mrstrat = new MultiRingMulticast
            (realSpacePlaneProxy.ckGetArrayID(), 1);
        
        DirectMulticastStrategy *d1strat = new DirectMulticastStrategy
            (particlePlaneProxy.ckGetArrayID(), 1);

        RingMulticastStrategy *r1strat = new RingMulticastStrategy
            (particlePlaneProxy.ckGetArrayID(), 1);

        MultiRingMulticast *mr1strat = new MultiRingMulticast
            (particlePlaneProxy.ckGetArrayID(), 1);

	if(CkNumNodes()>64) //multiring should be good on large runs, but not on BG/L
	  {
	      mcastInstance.setStrategy(dstrat);
//	      mcastInstancePP.setStrategy(d1strat);
	      mcastInstancePP.setStrategy(mr1strat);
	    }
	else
	  {
	      mcastInstance.setStrategy(rstrat);
	      mcastInstancePP.setStrategy(r1strat);
	  }
    }        
    // end Sameer's new communication strategies 

//============================================================================
// clean up

    delete m;
    delete sim;

//============================================================================


    CkPrintf("\n------------------------------------------------\n");
    CkPrintf("Cpaimd-Charm-Driver setup phase complete\n");
    CkPrintf("================================================\n\n");

//--------------------------------------------------------------------------
   }// end Main
//============================================================================


//============================================================================
//Create the array elements for the GSpace, Particle and Real Space planes
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void init_planes(size2d sizeYZ, int natm_nl,int natm_nl_grp_max,int numSfGrps,
                 int doublePack,int gSpacePPC,int realSpacePPC,
                 CPcharmParaInfo *sim)

//============================================================================
   { //begin routine 
//============================================================================
    /*
     * Set up the map variables used to control the location of the
     * 2d-array elements over processors   
     */

  PRINT_LINE_STAR;
  PRINTF("Building G-space Chares\n");
  PRINT_LINE_DASH;printf("\n");

    CProxy_GSMap gsMap;
    CProxy_RSMap rsMap;

    rsMap = CProxy_RSMap::ckNew();

    CkArrayOptions realSpaceOpts;
    realSpaceOpts.setMap(rsMap);
    size2d sizeRealPlane(sizeYZ[1], sizeX);

    realSpacePlaneProxy = CProxy_CP_State_RealSpacePlane::ckNew(sizeRealPlane,gSpacePPC,
                                                                realSpacePPC, 
                                                                realSpaceOpts);

								
    fftCacheProxy = CProxy_FFTcache::ckNew(sizeRealPlane, realSpacePPC);
    sfCacheProxy = CProxy_StructFactCache::ckNew(numSfGrps,natm_nl,natm_nl_grp_max);
    sfCompProxy = CProxy_StructureFactor::ckNew();
    
    gsMap = CProxy_GSMap::ckNew(sim->nplane_x, sim->lines_per_plane, sim->pts_per_plane);

    CkArrayOptions gSpaceOpts;
    gSpaceOpts.setMap(gsMap);

    gSpacePlaneProxy = CProxy_CP_State_GSpacePlane::ckNew(sizeX, sizeYZ, gSpacePPC, 
                                                          realSpacePPC, 
                                                          config.sGrainSize, gSpaceOpts);

    // We bind the particlePlane array to the gSpacePlane array migrate together
    CkArrayOptions particleOpts;
    particleOpts.setMap(gsMap); // the maps for both the arrays are the same
    particleOpts.bindTo(gSpacePlaneProxy);
    particlePlaneProxy = CProxy_CP_State_ParticlePlane::ckNew(sizeX, sizeYZ[0], sizeYZ[1],   
		      gSpacePPC,numSfGrps,natm_nl,natm_nl_grp_max,particleOpts);

    /*
     * Insert the planes in the particle plane array, gSpacePlane array
     */
    gSpacePlaneProxy.setReductionClient(doneInit, (void *) NULL);
    realSpacePlaneProxy.setReductionClient(doneInit, (void *) NULL);
    int s,x;
    for (s = 0; s < nstates; s++){
      for (x = 0; x < sizeX; x += gSpacePPC){
        if (x < config.low_x_size || ((x > config.high_x_size) &&(!doublePack) ) ) {
             gSpacePlaneProxy(s, x).insert(sizeX, sizeYZ, gSpacePPC, 
                                    realSpacePPC,config.sGrainSize);
             particlePlaneProxy(s, x).insert(sizeX, sizeYZ[0], sizeYZ[1],   
				      gSpacePPC,numSfGrps,natm_nl,natm_nl_grp_max);
	}
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
        for (y = 0; y < sizeYZ[0]; y += realSpacePPC)
            realSpacePlaneProxy(s, y).insert(sizeRealPlane, gSpacePPC, 
                                             realSpacePPC);
    realSpacePlaneProxy.doneInserting();
    // 
    int nplanes=sizeX;
    int *nsend= new int[nplanes];
    int **listpe = new int * [nplanes];

    for(int j=0;j<nplanes;j+=gSpacePPC){   
      listpe[j]= new int[nstates];
      nsend[j]=0;
//      GSMap gmap(sim->nplane_x, sim->lines_per_plane, sim->pts_per_plane);
      if (j < config.low_x_size || ((j > config.high_x_size) &&(!doublePack))) 
	{
	  for(int i=0;i<nstates;i++)    
	    {

	      listpe[j][i]=cheesyhackgsprocNum(sim, i,j);
//	      listpe[j][i]=gmap.slowprocNum(0, CkArrayIndex2D(i,j));
//	      listpe[j][i]=gmap.procNum(0, CkArrayIndex2D(i,j));
//	      CkPrintf("[%d %d] pe %d\n",j,i,listpe[j][i]);
	    }

	  lst_sort_clean(nstates, &nsend[j], listpe[j]);
	}
    }
    int minsend=nstates;
    int maxsend=0;
    double avgsend=0.0;
    int nplane_use=0;
    CkPrintf("============================\n");
    CkPrintf("Structure factor plane dests\n");
    CkPrintf("---------------------------\n");
    for(int lsi=0;lsi<nplanes;lsi++){
      if (lsi < config.low_x_size || ((lsi > config.high_x_size) &&(!doublePack) ) ) {
        nplane_use++;
	CkPrintf("plane [%d] nsend %d\n",lsi,nsend[lsi]);
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
      }//endif
    }//endfor
    CkPrintf("---------------------------\n");
    CkPrintf("SFSends min %d max %d avg %g\n",minsend,maxsend,avgsend/(double)nplane_use);
    CkPrintf("============================\n");
    // Insert the objects into the StructureFactor array
    for (int dup=0; dup<config.numSfDups; dup++)
      for (x = 0; x < sizeX; x += gSpacePPC){
	if (x < config.low_x_size || ((x > config.high_x_size) &&(!doublePack) ) ) 
	  {
	    int num_dup, istart, iend;
	    get_grp_params( nsend[x],  config.numSfDups,  dup, x ,&num_dup,  &istart, &iend);
	    int pe_ind=istart;
	    if(x%2==0)
	      pe_ind=iend;
	    for (int AtmGrp=0; AtmGrp<numSfGrps; AtmGrp++)
	      {
//		CkPrintf("inserting SF[%d %d %d] to pe %d\n",AtmGrp,x,dup,listpe[x][pe_ind]);
		sfCompProxy(AtmGrp, x, dup).insert(numSfGrps, config.numSfDups, natm_nl_grp_max,  num_dup, &(listpe[x][istart]),atmGrpMap(istart, num_dup, nsend[x], listpe[x], AtmGrp, dup,x));	      
		pe_ind++;
		if(pe_ind>nsend[x])
		  pe_ind=0;
	      }
	  }
      }
    sfCompProxy.doneInserting();

    for(int j=0;j<nplanes;j++)
      delete [] listpe[j];
    delete [] listpe;
    if(config.useCommlib) {
        // Set some com strategy of Sameer
        ssInstance = CkGetComlibInstance();
        StreamingStrategy *strat = new StreamingStrategy(0.2,10);
        //strat->enableShortArrayMessagePacking();
        ssInstance.setStrategy(strat);
        
        // Set some com strategy of Sameer
        mssInstance = CkGetComlibInstance();
        StreamingStrategy *mstrat = new StreamingStrategy(0.2,5);
        //mstrat->enableShortArrayMessagePacking();
        mssInstance.setStrategy(mstrat);    
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

void init_rho(size2d sizeYZ, int gSpacePPC, int realSpacePPC, int rhoGPPC)

//============================================================================
    {//begin routine
//============================================================================
/*
 * create the array for real-space densities (one-dimensional chare array)
 * and hartree energy computation
 */    
//============================================================================
    
        int pe = 0;
        CkAssert(config.rhoGPPC == 1);
        rhoGHelperProxy = CProxy_CP_Rho_GSpacePlaneHelper::ckNew();
	rhoGHelperProxy.setReductionClient(printEnergyHart, NULL);

        int z, y;
        int helperSize = sizeYZ[0]/config.rhoGHelpers;

        if( (sizeYZ[0] % config.rhoGHelpers) !=0 ){
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         CkPrintf("Helper size must be a mod of %d.\n",sizeYZ[0]);
         CkPrintf("Please fix your cpaimd_config.\n");
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
         CkExit();
	}//endif

        for (z = 0; z < sizeYZ[1]; z++){
            for (y = 0; y < sizeYZ[0]; y += helperSize) {
                rhoGHelperProxy(z,y).insert(sizeX, sizeYZ, y, pe, 0);
                pe = (pe + 1)%CkNumPes();
            }
        }
        rhoGHelperProxy.doneInserting();

    
    int pestride =  CkNumPes() / (sizeYZ[0]/realSpacePPC);
    if(pestride < 1)
      pestride = 1;
    
    if (sizeYZ[0]/realSpacePPC == sizeYZ[1]/rhoGPPC ) {
      rhoGProxy = CProxy_CP_Rho_GSpacePlane::ckNew();
      rhoRealProxy = CProxy_CP_Rho_RealSpacePlane::ckNew();
     
      int peg = 0;
      int per = pestride/2;
      int i;

      bool fftuseCommlib = config.fftuseCommlib;
      ComlibInstanceHandle fftcommInstance;
      if (fftuseCommlib) {        
	  int period_in_ms = 1, nmsgs = 1000;
	  StreamingStrategy * strat = new StreamingStrategy(period_in_ms, nmsgs);
	  fftcommInstance = CkGetComlibInstance();
	  fftcommInstance.setStrategy(strat);
      }
     for (i = 0; i < sizeYZ[1]/rhoGPPC; i++) {
	rhoGProxy[i].insert(sizeX, sizeYZ, realSpacePPC, rhoGPPC, fftuseCommlib, 
                            fftcommInstance, peg);
	rhoRealProxy[i].insert(sizeX, sizeYZ, realSpacePPC,
			       rhoGPPC, fftuseCommlib, fftcommInstance, per);
	
	peg += pestride;
	per += pestride;
                                                                                   
	if (peg >= CkNumPes())
	  peg = 0;
                                                                                   
	if (per >= CkNumPes())
	  per = pestride/2;
      }
                                                                                   
      rhoGProxy.doneInserting();
      rhoRealProxy.doneInserting();
    }else{ 
      CkAbort("sizeYZ[0]/realSpacePPC != sizeYZ[1]/rhoGPPC");
    }
    
    rhoRealProxy.setReductionClient(printEnergyEexc, 0);

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
    int natm    = PhysicsAtom->natm_tot;
    int natm_nl = PhysicsAtom->natm_nl;
    Atom *atoms = new Atom[natm];
    PhysicsAtom->DriverAtomInit(atoms);

    atomsGrpProxy = CProxy_AtomsGrp::ckNew(natm,natm_nl,atoms);
    delete [] atoms;
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
void get_grp_params(int natm_nl, int numSfGrps, int indexSfGrp, int planeIndex,
		    int *n_ret, int *istrt_ret, int *iend_ret)
{

   int n     = (natm_nl/numSfGrps);
   int m     = (natm_nl % numSfGrps);

   int istrt = n*indexSfGrp;
   if(indexSfGrp>=m){istrt += m;}
   if(indexSfGrp<m) {istrt += indexSfGrp;}
   if(indexSfGrp<m) {n++;}
   int iend  = n+istrt;
   if(numSfGrps>natm_nl)
     if(m>=indexSfGrp)
       {
	 n=0;
	 istrt=natm_nl+1;
	 iend=natm_nl;
         CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@\n");
	 CkPrintf("Redundant DupSF for plane %d\n",planeIndex);
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


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
int cheesyhackgsprocNum(CPcharmParaInfo *sim,int state, int plane)
{
  CkArrayIndex2D idx2d(state,plane);

  int pe = 0;
  int numPlanes = 0;
  
  if(config.doublePack) 
      numPlanes = sizeX/4;
  else
      numPlanes = sizeX/2;
  
  //pe = basicMap(idx2d, numPlanes);
  //return pe;
  double state_load=0.0;
  if(state_load <= 0.0 ) {
    for(int x = 0; x  < numPlanes; x++) {
      double curload = 0.0;
      double sload = 0.0;
      
      hackGSpacePlaneLoad(sim, x, &curload, &sload);
      
      state_load += curload;
    }
  }


  int pes_per_state = config.GpesPerState;
  int np = CkNumPes()/pes_per_state;
  
  if(np < 1)
    np = 1;
  
  int partition_nstates = config.nstates / np;
  if(config.nstates % np != 0)
    partition_nstates ++;
  
  int start_pe = (idx2d.index[0]/partition_nstates) * pes_per_state;
  int start_state = (idx2d.index[0]/partition_nstates) * partition_nstates;
  
  double cum_load = 0.0;
  double average_load = state_load * config.nstates / CkNumPes();

  for(int x = 0; x < numPlanes; x++) 
    for(int s = 0; s < partition_nstates; s++) {
      double curload = 0.0;
      double sload = 0.0;
      
      hackGSpacePlaneLoad(sim,x, &curload, &sload);

      cum_load += curload;

      if((idx2d.index[0] == s + start_state) && idx2d.index[1] == x) {
      
	//if(CkMyPe() == 0)
	//  CkPrintf("Load[%d] = %g\n", pe, load[pe]);
	
	double dpe = 0.0;
	dpe = cum_load / average_load;

	pe = (int)dpe;
	pe = pe % pes_per_state;
	pe += start_pe;

	return pe % CkNumPes();
      }
    }
  
  //  CkPrintf("Warning pe not found for index [%d, %d]\n",idx2d.index[0],idx2d.index[1]);
  return (idx2d.index[0]*1037+idx2d.index[1])%CkNumPes();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void hackGSpacePlaneLoad(CPcharmParaInfo *sim,int idx, double *line_load, 
                         double *pt_load){
//============================================================================

  int nplane_x         = sim->nplane_x;
  double *lines_per_plane = 
                 sim->lines_per_plane;
  double *pts_per_plane = 
                 sim->pts_per_plane;

  if( (idx < 0) || (idx >= nplane_x) ){
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkPrintf("GSpace plane index %d out of range: 0 < idx > %d\n",idx,nplane_x);
   CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
   CkExit();
  }//endif

  line_load[0] = lines_per_plane[idx];
  pt_load[0]   =   pts_per_plane[idx];

//============================================================================
  }//end routine
//============================================================================


//============================================================================
#include "cpaimd.def.h"
//============================================================================
