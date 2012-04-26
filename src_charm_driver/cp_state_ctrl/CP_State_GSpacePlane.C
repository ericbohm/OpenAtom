//======================================================
// Things to do : 
//    move resetiterstate
//======================================================
//#define _CP_DEBUG_WARN_SUSPEND_
//#define _CP_DEBUG_ORTHO_OFF_
//#define _CP_DEBUG_PSI_OFF_
//#define DEBUG_CP_GSPACE_PSIV
//#define BARRIER_CP_GSPACE_PSI
//#define BARRIER_CP_GSPACE_PSIV
//#define BARRIER_CP_GSPACE_NONLOCAL
//#define BARRIER_CP_GSPACE_IFFT
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_State_GSpacePlane.C
 * This is a description of the 'life' of a CP_State_GSpacePlane object.
 *
 *  At the beginning of the program, the constructor CP_State_GSpacePlane() is 
 *  called, to initialize the CP_State_GSpacePlane array. The GSpaceSlab within 
 *  the CP_State_GSpacePlane is initialized using the initGSpace(...) method.
 *
 *  To start off an iteration of program execution, the method doFFT() is 
 *  called. As part of this method, the CP_State_GSpacePlane object does a 
 *  forward 1-D fft on its slab of data, and sends off this data to the next
 *  stage in the computational loop. After this, the CP_State_GSpacePlane is 
 *  idle, waiting for a message to trigger some computation/communication.
 *
 *  The idle period of the CP_State_GSpacePlane is terminated by a called to 
 *  the method doIFFT(). In this method, the CP_State_GSpacePlane receives the 
 *  partially processed data from the CP_State_RealSpacePlanes and performs 
 *  1-D inverse FFT on its slab of data. Then the calculation of the "S" 
 *  matrix is started by  calling the sendPsi() method
 *
 *  After the forces are calculated using the inverse FFT, a check is done 
 *  to see if the forces from the particle calculations are ready. If so, the 
 *  forces are added up. The sendPsi() method is not called until the forces
 *  from the particle calculations and the forces from the quantum computation
 *  are ready.
 *
 *  The object is idle until the corrected g-space data from orthonormalization
 *  is received through the acceptNewPsi() method.
 */
//============================================================================


//============================================================================
#include "CP_State_GSpacePlane.h"
#include "CP_State_ParticlePlane.h"

#include "main/startupMessages.h"
#include "utility/util.h"
#include "main/AtomsCache.h"
#include "main/energyGroup.h"
#include "main/eesCache.h"
#include "main/TimeKeeper.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "structure_factor/StructFactorCache.h"
#include "main/CPcharmParaInfoGrp.h"
#include "main/cpaimd.h"
#include "main/InstanceController.h"
#ifdef PC_USE_RDMA
    #define ENABLE_RDMA_HANDSHAKES
#endif
#include "paircalc/RDMAMessages.h"

#include "charm++.h"

#include <iostream>
#include <fstream>
#include <cmath>

//---------------------------------------------------------------------------
#define CHARM_ON
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "../../src_piny_physics_v1.0/include/class_defs/allclass_gen.h"
#include "../../src_piny_physics_v1.0/include/class_defs/allclass_cp.h"
#include "../../src_piny_physics_v1.0/include/class_defs/PINY_INIT/PhysicsParamTrans.h"
//============================================================================


//============================================================================
extern Config config;
extern CProxy_main                    mainProxy;
extern CProxy_InstanceController      instControllerProxy;
extern CProxy_TimeKeeper              TimeKeeperProxy;
extern CkVec <CProxy_CP_State_RealSpacePlane> UrealSpacePlaneProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>    UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>            UgSpaceDriverProxy;
extern CkVec <CProxy_CP_State_ParticlePlane>  UparticlePlaneProxy;
extern CkVec <CProxy_AtomsCache>              UatomsCacheProxy;
extern CkVec <CProxy_StructureFactor>         UsfCompProxy;
extern CkVec <CProxy_EnergyGroup>             UegroupProxy;
extern CkVec <CProxy_FFTcache>                UfftCacheProxy;
extern CkVec <CProxy_StructFactCache>         UsfCacheProxy;
extern CkVec <CProxy_eesCache>                UeesCacheProxy;


extern CProxy_ComlibManager mgrProxy;
extern ComlibInstanceHandle gssInstance;
extern CkGroupID mCastGrpId;

extern int nstates;
extern int sizeX;
extern int nchareG;              // number of g-space chares <= sizeX and >=nplane_x

// Temporary global readonlys to hold the MeshStreamer group proxies
extern CProxy_ArrayMeshStreamer<streamedChunk, CkArrayIndex2D> fftStreamer;
extern CProxy_CompletionDetector completionDetector;

void testeke(int ,complex *,int *,int *,int *, int ,int);

//#define _CP_DEBUG_STATEG_VERBOSE_
//#define _CP_DEBUG_WARN_SUSPEND_



//============================================================================
// Entry method to resume execution after computing reduction over all planes
// and states to form psiCgOvlap (cg only) and magforPsi
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::psiCgOvlap(CkReductionMsg *msg){
//============================================================================
// Unpack
//  CkPrintf("{%d} GSP [%d,%d] psiCgOvlap\n",thisInstance.proxyOffset, thisIndex.x,thisIndex.y);
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  AtomsCache *ag         = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch(); // find me the local copy

  int cp_min_opt    = sim->cp_min_opt;
  double tol_cp_min = sim->tol_cp_min;
  double tol_cp_dyn = sim->tol_cp_dyn;
  int natm          = ag->natm;
  double rnatm      = ((double)natm)/96.0;  // relative to 32 waters

  double d0         = ((double *)msg->getData())[0];
  double d1         = ((double *)msg->getData())[1];
         d1         = sqrt(d1); // piny convention

  delete msg;  

//============================================================================
// Copy old/new : Set new values

  fovlap_old        = fovlap;          // CG ovlap (all chares need it)
  fmagPsi_total_old = fmagPsi_total;  

  fovlap            = d0;  
  fmagPsi_total     = d1;              // mag of psi force (all chares need it)

//============================================================================
// Output the mag force, send to the energy group, set the exit flag 

  if(thisIndex.x==0 && thisIndex.y==0){
      int iprintout   = iteration-1;

    fprintf(temperScreenFile, "Iter [%d] MagForPsi   =  %5.8lf | %5.8lf per entity\n", iprintout,d1,d1/rnatm);
    fprintf(temperScreenFile,"Iter [%d] Memory      =  %ld\n",iprintout,CmiMemoryUsage());
    computeEnergies(ENERGY_FMAG, d1);
  }//endif

  if(cp_min_opt==0 && fmagPsi_total>rnatm*tol_cp_dyn){
    exitFlag=1;
    if(thisIndex.x==0 && thisIndex.y==0){
      fprintf(temperScreenFile,"@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      fprintf(temperScreenFile, "Mag psi force %.10g > %.10g too large for CP dynamics. Ciao! \n",
	       fmagPsi_total/rnatm,tol_cp_dyn);
      fprintf(temperScreenFile,"@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
  }//endif
  int numBeads=config.UberImax;
  int numTempers=config.UberKmax;
#ifndef _CP_DEBUG_ORTHO_OFF_
  if (cp_min_opt==1)
    {
      if(numBeads==1 && numTempers==1)
	{
	  if(fmagPsi_total<=tol_cp_min){
#ifndef _CP_DEBUG_SCALC_ONLY_ 
	    exitFlag=1; outputFlag=1;
	    if(thisIndex.x==0 && thisIndex.y==0){
	      CkPrintf("----------------------------------------------\n");
	      CkPrintf("   CP wavefunction force tolerence reached!   \n");
	      CkPrintf("----------------------------------------------\n");
	    }//endif
	  }
#endif // _CP_DEBUG_SCALC_ONLY_ 
	}
      else
	{
#ifndef _CP_DEBUG_SCALC_ONLY_ 
	  if(fmagPsi_total<=tol_cp_min)
	    outputFlag=1;
	  if(exitFlagMin==1) // every bead and temper is minimized
	    {
	      exitFlag=1; outputFlag=1;
	      if(thisIndex.x==0 && thisIndex.y==0){
		CkPrintf("----------------------------------------------\n");
		CkPrintf("   CP wavefunction force tolerence reached!   \n");
		CkPrintf("----------------------------------------------\n");
	      }//endif
	    }
	  if(thisIndex.x==0 && thisIndex.y==0){
	    // can't let any bead stop until they all reach tolerance.
	    // but we only need one contributor from each replica.
	    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
	    int result=(fmagPsi_total <= tol_cp_min);
	    //	    CkPrintf("{%d} [%d,%d] tolcheck contrib %d %.5g %5g %5g\n",thisInstance.proxyOffset, thisIndex.x, thisIndex.y, result, fmagPsi_total, tol_cp_min,fmagPsi_total - tol_cp_min);
	    mcastGrp->contribute(sizeof(int), &result, CkReduction::logical_and, 
				 beadCookie);
	  }
#endif // _CP_DEBUG_SCALC_ONLY_ 
	}

    }
#endif

//============================================================================
// Do a little cputime management in GS class then resume

  if(thisIndex.x==0 && thisIndex.y==0){
#ifdef _CP_DEBUG_ORTHO_OFF_
      CkPrintf("==============================================\n");
      CkPrintf("        Running with PC/Ortho Off             \n");
#endif
     double cpuTimeOld = cpuTimeNow;
     cpuTimeNow        = CkWallTimer();
     if(iteration>1){
       fprintf(temperScreenFile, "Iter [%d] CpuTime(GSP)= %g\n",iteration-1,cpuTimeNow-cpuTimeOld);
       if(cp_min_opt==0){
         int heavyside = 1-(iteration-iterRotation >= 1 ? 1 : 0);
         fprintf(temperScreenFile, "Iter [%d] Step = %d : Step Last Rot = %d : Interval Rot = %d : Num Rot = %d : %d\n",iteration,
                   iteration,iterRotation,iteration-iterRotation,nrotation,heavyside);
       }//endif
     }//endif
  }//endif

  UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl();
//============================================================================
  }// end routine
//============================================================================
void CP_State_GSpacePlane::initBeadCookie(ICCookieMsg *m)
{
  //CkPrintf("{%d} [%d,%d] beadcookie initialized\n",thisInstance.proxyOffset, thisIndex.x, thisIndex.y);
  CkGetSectionInfo(beadCookie,m);
  //beadCookie=m->_cookie;
}

void CP_State_GSpacePlane::minimizeSync(ICCookieMsg *m)
{
  // CkPrintf("{%d} [%d,%d] minimizeSync %d\n",thisInstance.proxyOffset, thisIndex.x, thisIndex.y, m->junk);
  CkGetSectionInfo(beadCookie,m);
  if(m->junk==1)
    thisProxy.setExitFlag();
}

void CP_State_GSpacePlane::setExitFlag()
{
  exitFlagMin=1;
}
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void printForce(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d = ((double *)m->getData())[0];

  CkPrintf("printForces = %5.8f \n", d);
  delete m;
}
//============================================================================


//============================================================================
// int sizeX: # of planes in x per state 
// size2d size: size of plane (in y*z dimension)
// int GSpaceUnits: #planes per slab (Gspace: in x dimension)
// int realSpaceUnits: #planes per slab (Rspace: in y dimension)
// int s_grain: S matrix grain size
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_GSpacePlane::CP_State_GSpacePlane(int    sizeX, 
                                           int    gSpaceUnits, 
                                           int    realSpaceUnits, 
                                           int    s_grain,
					   int   _gforward,
					   int   _gbackward,
                       int   _fftFwd,
                       int   _fftBwd,
					   UberCollection _thisInstance
					   ) :
  forwardTimeKeep(_gforward),  backwardTimeKeep(_gbackward),
  fftFwdTimer(_fftFwd), fftBwdTimer(_fftBwd),
  thisInstance(_thisInstance)
//============================================================================
   {//begin routine
//============================================================================
//  ckout << "State G Space Constructor : "
//        << thisIndex.x << " " << thisIndex.y << " " <<CkMyPe() << endl;
//============================================================================

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt  = sim->cp_min_opt;
  int gen_wave    = sim->gen_wave;
  wallTimeArr=NULL;
  if(thisIndex.x==0 && thisIndex.y==0 && config.maxIter<30){
    wallTimeArr = new double[config.maxIter+2];
  }else{
    wallTimeArr = new double[30];
  }//endif
  wallTimeArr[0]=0.0;
  wallTimeArr[1]=0.0;

  myBeadIndex    = thisInstance.idxU.x;
  myKptIndex     = thisInstance.idxU.y;
  myTemperIndex  = thisInstance.idxU.z;
  mySpinIndex    = thisInstance.idxU.s;

//============================================================================

  istate_ind           = thisIndex.x;
  iplane_ind           = thisIndex.y;  
  ibead_ind            = thisInstance.idxU.x;
  kpoint_ind           = thisInstance.idxU.y;
  itemper_ind          = thisInstance.idxU.z;
  ispin_ind            = 0;                   //needs to be updated 
  initialized          = false;
  iteration            = 0;
  nrotation            = 0;
  iterRotation         = 0;
  gotHandles =0;

  total_energy      = 0.0;
  ehart_total       = 0.0;
  enl_total         = 0.0;
  eke_total         = 0.0;
  egga_total        = 0.0;
  eexc_total        = 0.0;
  eext_total        = 0.0;
  ewd_total         = 0.0;
  fovlap            = 0.0;
  fovlap_old        = 0.0;
  fmagPsi_total     = 0.0;
  fmagPsi_total0    = 0.0; // only chare(0,0) cares
  fmagPsi_total_old = 0.0;
  cpuTimeNow        = 0.0;
  fictEke_total     = 0.0;

  halfStepEvolve      = 1;
  ireset_cg           = 1;
  numReset_cg         = 0;
  exitFlag            = 0;
  exitFlagMin         = 0;
  outputFlag          = 0;
  iRecvRedPsi         = 1;  
  iSentRedPsi         = 1;
  iRecvRedPsiV        = 0;
  iSentRedPsiV        = 0;

  finishedCpIntegrate = 0;
  if(cp_min_opt==0){finishedCpIntegrate = 1;}// alternate entry point
  if(gen_wave==1){finishedCpIntegrate = 1;}// alternate entry point
  doneDoingIFFT       = false;
  isStreamerReady     = false;
  isForwardFftSendPending = false;
  acceptedPsi         = true;    // we start out with a psi
  acceptedVPsi        = true;    // we start out with a vpsi
  acceptedLambda      = false;   // no forces yet
  doneNewIter         = false;

  countFileOut    = 0;
  countRedPsi     = 0;
  countRedPsiV    = 0;
  countPsi        = 0;
  countLambda     = 0;
  countVPsi       = 0;
  countIFFT       = 0;
  ecount          = 0; //No energies have been received.
  cleanExitCalled = 0;

  countPsiO       = new int[config.numChunksSym];

  for(int i=0;i<config.numChunksSym;i++)
    countPsiO[i]=0;

  countVPsiO= new int[config.numChunksSym];
  for(int i=0;i<config.numChunksSym;i++)
    countVPsiO[i]=0;
  /* figure out how to record the off diag extra lambdas */
  /* some trickery like psi uses */
 
 
  countLambdaO= new int[config.numChunksAsym];
  for(int i=0;i<config.numChunksAsym;i++)
    countLambdaO[i]=0;

  numRecvRedPsi   = sim->RCommPkg[thisIndex.y].num_recv_tot;

 //---------------------------------------------
 //@todo: RV - try to understand the logic behind this calculation of All Psi/Lambda Expected
 // Symm PC accounting 
  int remainder=nstates%s_grain;
  int sizeoflastgrain=s_grain+remainder;
  int lastgrain=nstates-sizeoflastgrain;
  int ourgrain    = thisIndex.x/s_grain*s_grain; 
  if(nstates == s_grain){
     AllPsiExpected=1;
  }else{ 
    //    if(ourgrain<lastgrain){ // corner has no extras
    if(ourgrain<(nstates-s_grain)){ // corner has no extras
       AllPsiExpected=2;
    }else{
       AllPsiExpected=1;
    }//endif
  }//endif
  AllPsiExpected*=config.numChunksSym;

  int numGrains=nstates/s_grain;
  if(config.gSpaceSum){ // no reductions its all coming direct
      AllPsiExpected=numGrains*config.numChunksSym;
  }//endif

 //---------------------------------------------
 // Asymm PC accounting 
  if(cp_min_opt==0){    //expect non diagonal column results
    if(nstates == s_grain){
	AllLambdaExpected=1;
    }else {
	AllLambdaExpected=2;
    }//endif
  }else{  //no column reduction results in minimization
      AllLambdaExpected=1;
  }//endif

  if(config.gSpaceSum){ // no reductions its all coming direct
    if(cp_min_opt==0 && numGrains>1) 
      { AllLambdaExpected=(2*numGrains-1)*config.numChunksAsym;}
    else
      { AllLambdaExpected=numGrains*config.numChunksAsym*AllLambdaExpected;}
  }//endif
  else
    {
      AllLambdaExpected*=config.numChunksAsym;

    }

	/// Compute the number of RDMA links that I'll have with the symm/asymm PC chares
	#ifdef PC_USE_RDMA
		numRDMAlinksSymm  = numGrains * config.numChunksSym;
		numRDMAlinksAsymm = numGrains * config.numChunksAsym * 2;
	#else
		numRDMAlinksSymm  = 0;
		numRDMAlinksAsymm = 0;
	#endif
	
//============================================================================
// Just zero everything for now

  gSpaceNumPoints = 0;
  tpsi           = NULL;
  tvpsi          = NULL;
  tk_x           = NULL;
  tk_y           = NULL;
  tk_z           = NULL;
  int len_nhc_cp;
  int num_nhc_cp;
  int nck_nhc_cp;
  CPINTEGRATE::fetchNHCsize(&len_nhc_cp,&num_nhc_cp,&nck_nhc_cp);
  initGStateSlab(&gs,sizeX,sim->sizeY,sim->sizeZ,gSpaceUnits,realSpaceUnits,s_grain,
                 thisIndex.y,thisIndex.x,len_nhc_cp,num_nhc_cp,nck_nhc_cp);


//============================================================================
// Load Balancing etc

  usesAtSync = CmiTrue;
  if(config.lbgspace){
    setMigratable(true);
  }else{
    setMigratable(false);
  }//endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
   savedvksBf=NULL;
   savedforceBf=NULL;
#endif
#ifdef  _CP_GS_DEBUG_COMPARE_PSI_
   savedpsiBf=NULL;
   savedpsiBfp=NULL;
   savedpsiAf=NULL;
   savedlambdaBf=NULL;
   savedlambdaAf=NULL;
#endif

//============================================================================
// Pick a reduction plane

  redPlane = 0;
  if(nchareG>1){
     redPlane = 1;
#ifdef _FANCY_PLANES_
    if(config.numPes>2*nstates){
      CkVec <int> usedVec;
      CkVec <int> peUsedByNLZ;
      CkVec <int> planeUsedByNLZ;
      FILE *fp;
      if(thisIndex.x+1==nstates){fp=fopen();}
      for(int state=0; state<thisIndex.x;state++){
        redPlane=nchareG-1;
        while(redPlane>=0){
          bool used=false;
          int thisstateplaneproc=GSImaptable.get(state,redPlane)%CkNumPes();
          for(int i=0;i<usedVec.size();i++){
	    if(usedVec[i]==thisstateplaneproc){used=true;}
          }//endfor
	  if(!used || redPlane==0){
	      peUsedByNLZ.push_back(thisstateplaneproc);
	      planeUsedByNLZ.push_back(redPlane);
	      usedVec.push_back(thisstateplaneproc);
	      redPlane=-1;
	  }//endif
	  redPlane--;
        }//endwhile
      }//endfor
    }else{
      redPlane = (thisIndex.x % (nchareG-1));
    }//endif
#else 
      redPlane = (thisIndex.x % (nchareG-1));
#endif
    redPlane = (redPlane < 0 ? redPlane+nchareG : redPlane);
    redPlane = (redPlane > nchareG-1 ? redPlane-nchareG : redPlane);
  }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_GSpacePlane::CP_State_GSpacePlane(CkMigrateMessage *m) {}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_GSpacePlane::~CP_State_GSpacePlane(){
  if(initialized) {
    delete [] lambdaproxyother;
    delete [] lambdaproxy;
    delete [] psiproxyother;
    delete [] psiproxy;
  }//
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::pup(PUP::er &p) {
//============================================================================
  ArrayElement2D::pup(p);
  UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).pup(p);
  //control flags and functions reference by thread are public
  p|istate_ind;
  p|iplane_ind;
  p|ibead_ind; p|kpoint_ind; p|itemper_ind; p|ispin_ind;
  p|registrationFlag;
  p|initialized;
  p|istart_typ_cp;
  p|iteration;
  p|nrotation;
  p|exitFlag;
  p|exitFlagMin;
  p|outputFlag;
  p|cleanExitCalled;
  p|finishedCpIntegrate;
  p|iRecvRedPsi;
  p|iSentRedPsi;
  p|iRecvRedPsiV;
  p|iSentRedPsiV;
  p|ireset_cg;
  p|numReset_cg;
  p|countPsi;
  p|countVPsi;
  p|countLambda;
  p|countIFFT;
  p|countFileOut;
  p|ecount;
  p|countRedPsi;
  p|numRecvRedPsi;
  p|AllPsiExpected;
  p|AllLambdaExpected;
  p|numRDMAlinksSymm;
  p|numRDMAlinksAsymm;
  p|doneDoingIFFT;
  p|isStreamerReady;
  p|isForwardFftSendPending ;
  p|doneNewIter;
  p|acceptedPsi;
  p|acceptedVPsi;
  p|acceptedLambda;
  p|itemp; // 2 temporary variables for debugging in scope
  p|jtemp;
  p|myBeadIndex;
  p|myKptIndex;
  p|myTemperIndex;
  p|mySpinIndex;

  p|ehart_total;
  p|enl_total;
  p|eke_total;
  p|fictEke_total;
  p|fmagPsi_total0;
  p|fmagPsi_total;
  p|fmagPsi_total_old;
  p|fovlap;
  p|fovlap_old;
  p|egga_total;
  p|eexc_total;
  p|eext_total;
  p|ewd_total;
  p|total_energy;
  p|cpuTimeNow;

  p|real_proxy;   
  gs.pup(p);
  p|gSpaceNumPoints;
  if (p.isUnpacking()) {
    lambdaproxy=new CProxySection_PairCalculator[config.numChunksAsym];
    lambdaproxyother=new CProxySection_PairCalculator[config.numChunksAsym];
    psiproxy=new CProxySection_PairCalculator[config.numChunksSym];
    psiproxyother=new CProxySection_PairCalculator[config.numChunksSym];
    countPsiO= new int[config.numChunksSym];
    countVPsiO= new int[config.numChunksSym];
    countLambdaO= new int[config.numChunksAsym];
  }//endif  
  PUParray(p,lambdaproxy,config.numChunksAsym);
  PUParray(p,lambdaproxyother,config.numChunksAsym);
  PUParray(p,psiproxy,config.numChunksSym);
  PUParray(p,psiproxyother,config.numChunksSym);
  PUParray(p,countPsiO,config.numChunksSym);
  PUParray(p,countVPsiO,config.numChunksSym);
  PUParray(p,countLambdaO,config.numChunksAsym);
//-------------------------------------------------------
   }// end routine : pup
//============================================================================




void CP_State_GSpacePlane::acceptPairCalcAIDs(pcSetupMsg *msg)
{
    symmPCmgr  = PCCommManager(thisIndex, msg->symmCfg, msg->symmIDs);
    asymmPCmgr = PCCommManager(thisIndex, msg->asymmCfg, msg->asymmIDs);
    myOrtho    = CProxy_Ortho(msg->orthoAID);

//============================================================================
// Contribute to the reduction telling main we are done

  int constructed=1;
  contribute(sizeof(int), &constructed, CkReduction::sum_int, 
	     CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy), thisInstance.proxyOffset);
}




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// In this function data is read from files, and sent to the corresponding
// G-space planes. Data reading will be done in chunk 0 of each state
//============================================================================
void CP_State_GSpacePlane::readFile() {
//============================================================================
// Local pointers

  CP           *cp           = CP::get();
#include "../class_defs/allclass_strip_cp.h"
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int numData       = config.numData;
  int ibinary_opt   = sim->ibinary_opt;
  int istart_typ_cp = sim->istart_typ_cp;
  int gen_wave      = sim->gen_wave;
  int ncoef         = sim->ncoef;
  CkVec <RunDescriptor> *sortedRunDescriptors = sim->sortedRunDescriptors;
  int *npts_lgrp    = sim->npts_per_chareG;
  int *nline_lgrp   = sim->nlines_per_chareG;
  int *istrt_lgrp   = NULL;
  int *iend_lgrp    = NULL;

  int nx            = sizeX;
  int ny            = sim->sizeY;
  int nz            = sim->sizeZ;
  int ngridaNL      = sim->ngrid_nloc_a;
  int ngridbNL      = sim->ngrid_nloc_b;
  int ngridcNL      = sim->ngrid_nloc_c;
  int nkpoint       = sim->nkpoint;
  int cp_force_complex_psi = sim->cp_force_complex_psi;

//============================================================================
// Set the file name using the config path and state number

  char fname[1024];
  int ind_state=thisIndex.x;
  //------------------------------------------------------------------
  // Get the complex data, Psi(g) and the run descriptor (z-lines in g-space)

  complex *complexPoints  = (complex *)fftw_malloc(numData*sizeof(complex));
  complex *vcomplexPoints = NULL;
  if(istart_typ_cp>=3){vcomplexPoints = (complex *)fftw_malloc(numData*sizeof(complex));}

  int *kx=  (int *)fftw_malloc(numData*sizeof(int));
  int *ky=  (int *)fftw_malloc(numData*sizeof(int));
  int *kz=  (int *)fftw_malloc(numData*sizeof(int));
  int nlines_tot,nplane;

  if(istart_typ_cp>=3){
    sprintf(fname, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/vState%d.out",
            config.dataPath,mySpinIndex,myKptIndex,myBeadIndex,myTemperIndex,ind_state+1);
    readState(numData,vcomplexPoints,fname,ibinary_opt,&nlines_tot,&nplane, 
            kx,ky,kz,&nx,&ny,&nz,istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,0,1);
  }//endif

  if(gen_wave==0){
    sprintf(fname, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/state%d.out",
            config.dataPath,mySpinIndex,myKptIndex,myBeadIndex,myTemperIndex,ind_state+1);
    readState(numData,complexPoints,fname,ibinary_opt,&nlines_tot,&nplane, 
              kx,ky,kz,&nx,&ny,&nz,istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,0,0);
    if(cp_force_complex_psi==1){
      if(ind_state==0){
        CkPrintf("\n$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n");
        CkPrintf("Adding a phase to the states to debug kpt code!!\n");
        CkPrintf("in routine CP_State_GspacePlane.C\n");
        CkPrintf("$$$$$$$$$$$$$$$$$$$$_warning_$$$$$$$$$$$$$$$$$$$$\n\n");
      }//endif
      double phase = M_PI*((double)(ind_state+1))/((double)(nstates+1));
      for(int i=0;i<numData;i++){
        double re = cos(phase); double im = sin(phase);
        double ore = re*complexPoints[i].re - im*complexPoints[i].im; 
        double oim = re*complexPoints[i].im + im*complexPoints[i].re; 
        complexPoints[i].re = ore;
        complexPoints[i].im = oim;
      }//endfor
    }//endif
  }else{
    kx -= 1;  ky -= 1; kz -=1;
    PhysicsParamTransfer::fetch_state_kvecs(kx,ky,kz,ncoef,config.doublePack);
    kx += 1;  ky += 1; kz +=1;
    processState(numData,ncoef,complexPoints,fname,ibinary_opt,&nlines_tot,
                 &nplane,kx,ky,kz,istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,
                 0,0,0,ny); 
    double *xfull =UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->fastAtoms.x-1;
    double *yfull =UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->fastAtoms.y-1;
    double *zfull =UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->fastAtoms.z-1;
    cpgen_wave->create_coefs(kx,ky,kz,numData,ind_state,complexPoints,
                             xfull,yfull,zfull,kpoint_ind);
  }//*endif

  if(config.nGplane_x != nplane && config.doublePack){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Mismatch in planesize\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

// Test the input g-vectors and psi(g) using kinetic energy
//  testeke(numData,complexPoints,kx,ky,kz,1,ind_state);

//============================================================================
// Blast off the data to chares

  int ioff = 0;
  for(int x = 0; x < nchareG; x ++){

      int numPoints    = 0;
      for (int j = 0; j < sortedRunDescriptors[x].size(); j++){
	numPoints += sortedRunDescriptors[x][j].length;
      }//endfor

      complex *dataToBeSent  = (complex *)fftw_malloc(numPoints*sizeof(complex));
      complex *temp          = complexPoints+ioff;
      CmiMemcpy(dataToBeSent,temp,(sizeof(complex) * numPoints));

      int numPointsV;
      complex *vdataToBeSent;
      if(istart_typ_cp>=3){
         numPointsV          = numPoints;
         vdataToBeSent       = (complex *)fftw_malloc(numPointsV*sizeof(complex));
         complex *vtemp      = vcomplexPoints+ioff;
         CmiMemcpy(vdataToBeSent,vtemp,(sizeof(complex) * numPoints));
      }else{
         numPointsV          = 1;
         vdataToBeSent       = (complex *)fftw_malloc(numPointsV*sizeof(complex));
         vdataToBeSent[0].re = 0.0;
         vdataToBeSent[0].im = 0.0;
      }//endif

      if(ioff>numData){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Error reading\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif

      UgSpacePlaneProxy[thisInstance.proxyOffset](ind_state, x).initGSpace(
                                numPoints,dataToBeSent,numPointsV,vdataToBeSent,
				nx,ny,nz,ngridaNL,ngridbNL,ngridcNL,istart_typ_cp);
      fftw_free(dataToBeSent);
      fftw_free(vdataToBeSent);

      ioff += numPoints;
      CmiNetworkProgress();
  }//endfor : loop over all possible chares in g-space (pencils)

  CkAssert(numData==ioff);

//============================================================================
// Clean up

  fftw_free(complexPoints);
  if(istart_typ_cp>=3){fftw_free(vcomplexPoints);}
  fftw_free(kx);
  fftw_free(ky);
  fftw_free(kz);

//---------------------------------------------------------------------------
   }//read the file
//============================================================================


//============================================================================
/**
 * This method is used to accept the state data from some initializing
 * routine. Since a CP_State_GSpacePlane class can have more than one plane, 
 * the method used for initialization is as follows:
 * 
 * This call initializes the entire set of planes
 * runDescSize: number of run-descriptors
 * size: the total number of non-zero points
 * points: pointer to the total data.
n *
 * The data is copied into the planes.
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::initGSpace(int            size, 
                                      complex*       points,
                                      int            vsize, 
                                      complex*       vpoints,
                                      int nx, int ny, int nz, 
                                      int nxNL, int nyNL, int nzNL, 
                                      int istart_cp) 
//============================================================================
   { //begin routine
//============================================================================
#ifdef _CP_DEBUG_STATEG_VERBOSE_
    CkPrintf("GSpace[%d,%d] initGSpace %d\n",thisIndex.x,thisIndex.y,size);
#endif

    temperScreenFile = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->temperScreenFile;

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  registrationFlag  = 1;
  eesCache *eesData   = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  eesData->registerCacheGSP(thisIndex.x,thisIndex.y);

//============================================================================
// Section Reductions and SF proxy creation

  real_proxy = UrealSpacePlaneProxy[thisInstance.proxyOffset];

#ifdef USE_COMLIB
  if (config.useGssInsRealP){
     ComlibAssociateProxy(gssInstance,real_proxy);
  }//endif
#endif
  
//============================================================================
// Register with the cache : Eric's multiple reduction schemes ensure its done
//                           before we need it.

  int cp_min_opt    = sim->cp_min_opt;
  int cp_min_update = sim->cp_min_update;

  if (true == initialized) {
    ckerr << "Attempt to initialize a plane twice" << endl;
    return;
  }//endif
  initialized = true;

//============================================================================
// Setup gs

  istart_typ_cp  = istart_cp;

  gs.eke_ret     = 0.0;  
  gs.fictEke_ret = 0.0;  
  gs.ekeNhc_ret  = 0.0;  
  gs.cp_min_opt  = cp_min_opt;

  gs.numRuns     = eesData->GspData[iplane_ind]->numRuns;
  gs.numLines    = gs.numRuns/2;
  gs.numFull     = (gs.numLines)*nz;
  gs.numFullNL   = (gs.numLines)*nzNL;
  gs.istate_ind  = thisIndex.x;
  gs.iplane_ind  = thisIndex.y;
  gs.mysizeX     = sizeX;
  gs.fftReqd     = true;

  gs.xdim        = 1;
  gs.ydim        = ny;
  gs.zdim        = nz;

  gs.ngridaNL    = nxNL;
  gs.ngridbNL    = nyNL;
  gs.ngridcNL    = nzNL;

  gs.numPoints   = eesData->GspData[iplane_ind]->ncoef;;
  CkAssert(gs.numPoints == size);

  gs.ees_nonlocal        = sim->ees_nloc_on;

  gs.packedPlaneData     = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));
  gs.packedForceData     = (complex *)fftw_malloc(gs.numFull*sizeof(complex));
  gs.packedVelData       = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));
  CmiMemcpy(gs.packedPlaneData, points, sizeof(complex)*gs.numPoints);
  bzero(gs.packedForceData,sizeof(complex)*gs.numFull);

  if(cp_min_opt==0){
    gs.packedPlaneDataScr   = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));
    gs.packedPlaneDataTemp  = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));  
    memset(gs.packedPlaneDataScr, 0, sizeof(complex)*gs.numPoints);
    CmiMemcpy(gs.packedPlaneDataTemp, points, sizeof(complex)*gs.numPoints);
  }//endif

  if(cp_min_opt==1 && cp_min_update==0){
    gs.packedPlaneDataTemp = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));  
    CmiMemcpy(gs.packedPlaneDataTemp, points, sizeof(complex)*gs.numPoints);
  }//endif

#ifdef _CP_DEBUG_SCALC_ONLY_
  //  if(cp_min_opt==0){
    gs.packedPlaneDataTemp2 = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));  
    bzero(gs.packedPlaneDataTemp2,sizeof(complex)*gs.numPoints);
    //  }//endif
#endif

  // Under cp_min veldata is the conjugate gradient : always need it.
  if(istart_typ_cp>=3 && cp_min_opt==0){
    CkAssert(vsize == size);
    CmiMemcpy(gs.packedVelData, vpoints, sizeof(complex)*gs.numPoints);
  }else{
    memset(gs.packedVelData, 0, sizeof(complex)*gs.numPoints);
  }//endif

//============================================================================
// Setup gpspaceplane and particle plane

  CkAssert(UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal());
  CP_State_ParticlePlane *localParticlePlaneChare = UparticlePlaneProxy[thisInstance.proxyOffset](thisIndex.x, thisIndex.y).ckLocal();
  CkAssert(localParticlePlaneChare);
  localParticlePlaneChare->initKVectors();

//============================================================================
// Setup k-vector ranges, masses and zero the force overlap

  int *k_x          = eesData->GspData[iplane_ind]->ka;
  int *k_y          = eesData->GspData[iplane_ind]->kb;
  int *k_z          = eesData->GspData[iplane_ind]->kc;
  double *coef_mass = eesData->GspData[iplane_ind]->coef_mass;
  int mycoef        = eesData->GspData[iplane_ind]->ncoef;
  gSpaceNumPoints   = gs.numPoints;

  if(eesData->allowedGspChares[iplane_ind]==0 || mycoef != gSpaceNumPoints){
    CkPrintf("Plane %d of state %d toasy %d %d\n",iplane_ind,thisIndex.x,
              mycoef,gSpaceNumPoints);
    CkExit();
  }//endif

  gs.setKRange(gSpaceNumPoints,k_x,k_y,k_z);
  // if I have redundant coefs, I must receive them from someone (e.g. myself included)
  // numRecRedPsi is the number of chares that send master to coefs to me to overwrite my redundant guys
  if(numRecvRedPsi==0 && gs.nkx0_red>0){
    CkPrintf("Error in GSchare(%d %d) on proc %d : numRecvRedPsi=%d gs.nkx0_red=%d\n",thisIndex.x,thisIndex.y,
	     CkMyPe(),numRecvRedPsi,gs.nkx0_red);
    CkExit();
  }//endif

  fovlap      = 0.0; 

//============================================================================
// Init NHC, Sample velocities 

  int ncoef_use = gs.numPoints-gs.nkx0_red;
  int maxLenNHC = gs.len_nhc_cp; // double check
  int maxNumNHC = gs.num_nhc_cp; // double check
  int maxNckNHC = gs.nck_nhc_cp; // double check

  CPINTEGRATE::initCPNHC(ncoef_use,gs.nkx0_zero,maxLenNHC,maxNumNHC,maxNckNHC,
                         &gs.kTCP,&gs.tauNHCCP,&gs.degfree,&gs.degfreeNHC,
                         gs.mNHC,gs.v0NHC,gs.a2NHC,gs.a4NHC,
                         gs.degFreeSplt,gs.istrNHC,gs.iendNHC);

  if(cp_min_opt==0){
   CPINTEGRATE::CPSmplVel(gs.numPoints,coef_mass,gs.packedVelData,
                          gs.len_nhc_cp,gs.num_nhc_cp,gs.nck_nhc_cp,
                          gs.mNHC,gs.vNHC,gs.xNHC,gs.xNHCP,gs.a2NHC,
                          gs.kTCP,istart_typ_cp,gs.nkx0_red,gs.nkx0_uni,
                          gs.nkx0_zero,gs.degfree,gs.degfreeNHC);
  }//endif

#ifdef _CP_DEBUG_PSI_OFF_
  memset(gs.packedPlaneData, 0, sizeof(complex)*gs.numPoints);
#endif

#define _CP_DEBUG_DYNAMICS_ZVEL_OFF_
#ifdef _CP_DEBUG_DYNAMICS_ZVEL_
  bzero(gs.packedVelData,sizeof(complex)*gs.numPoints);
#endif

//============================================================================
// Send the k's to the structure factor 

  if(sim->ees_nloc_on == 0){
   for(int atm=0;atm<config.numSfGrps; atm++){ //each atm
    for(int dup=0;dup<config.numSfDups;dup++){ //each dup
      if(dup==thisIndex.x){
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("GSP [%d,%d] on PE %d sending KVectors to SF[%d,%d,%d]\n",
                    thisIndex.x, thisIndex.y, CkMyPe(), atm, thisIndex.y, dup);
#endif
        UsfCompProxy[thisInstance.proxyOffset](atm,thisIndex.y,dup).acceptKVectors(gSpaceNumPoints, k_x, k_y, k_z);
      }//endif
    }//endfor
   }//endfor
  }//endif

//============================================================================
//Some PC initialization that needs to happen here to avoid
//constructor race conditions

  makePCproxies();
#ifdef PC_USE_RDMA
  /** Now that we have paircalcids, proxies, and allocated data, setup RDMA with a handshake to each PC.
   * Give each PC a token as part of the handshake that identifies what we will be sending it.
   * The PC will fill the token with its ID information and return it to us.
   * This, along with the RDMA handle will allow us to complete the setup by pointing out the appropriate 
   * data to be sent to that PC.
   */
  
  /// Get the entry point index of the method that should be called to complete the RDMA handshake
  int rdmaConfirmEP = CkIndex_CP_State_GSpacePlane::completeRDMAhandshake(0);
  /// Create a callback object to this entry method
  CkCallback rdmaConfirmCB( rdmaConfirmEP, CkArrayIndex2D(thisIndex.x,thisIndex.y), thisProxy.ckGetArrayID() );
  /// Create a handshake token and fill it with my (sender) ID
  RDMApair_GSP_PC handshakeToken;
  handshakeToken.gspIndex = thisIndex;
  
  /// Send rdma requests to all asymmetric loop PCs who will get data from me
  handshakeToken.symmetric = false;
  handshakeToken.shouldSendLeft = true;
  sendLeftRDMARequest (handshakeToken,gs.numPoints,rdmaConfirmCB);
  handshakeToken.shouldSendLeft = false;
  sendRightRDMARequest(handshakeToken,gs.numPoints,rdmaConfirmCB);
  
  /// Send rdma requests to all symmetric loop PCs who will get data from me
  handshakeToken.symmetric = true;
  handshakeToken.shouldSendLeft = true;
  sendLeftRDMARequest (handshakeToken,gs.numPoints,rdmaConfirmCB);
  handshakeToken.shouldSendLeft = false;
  sendRightRDMARequest(handshakeToken,gs.numPoints,rdmaConfirmCB);
  
#else
//============================================================================
// This reduction is done to signal the end of initialization to main

  int i=1;
  contribute(sizeof(int), &i, CkReduction::sum_int, 
	     	       CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy), thisInstance.proxyOffset);

#endif
#ifdef _CP_SUBSTEP_TIMING_
#if USE_HPM
  (TimeKeeperProxy.ckLocalBranch())->initHPM();
#endif // HPM
#endif // _CP_SUBSTEP_TIMING_


//---------------------------------------------------------------------------
   
}// end routine
//============================================================================

void CP_State_GSpacePlane::acceptNewTemperature(double temp)
{
  // Hey GLENN do something with your new temperature here


  // when you're done
  int i=1;
  contribute(sizeof(int), &i, CkReduction::sum_int, 
	     	       CkCallback(CkIndex_InstanceController::gspDoneNewTemp(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy), thisInstance.proxyOffset);
}



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::makePCproxies(){

  lambdaproxy=new CProxySection_PairCalculator[config.numChunksAsym];
  lambdaproxyother=new CProxySection_PairCalculator[config.numChunksAsym];
  psiproxy=new CProxySection_PairCalculator[config.numChunksSym];
  psiproxyother=new CProxySection_PairCalculator[config.numChunksSym];
  //need one proxy per chunk
  if(!config.gSpaceSum){
      for(int chunk=0;chunk<config.numChunksAsym;chunk++){
	  lambdaproxy[chunk]=asymmPCmgr.makeOneResultSection_asym(chunk);
	  if(AllLambdaExpected/config.numChunksAsym == 2)//additional col. red. in dynamics
	    lambdaproxyother[chunk]=asymmPCmgr.makeOneResultSection_asym_column(chunk);
      }//endfor chunk
      for(int chunk=0; chunk < config.numChunksSym ;chunk++){
	  psiproxy[chunk]=symmPCmgr.makeOneResultSection_sym1(chunk);
	  if(AllPsiExpected / config.numChunksSym > 1)
	    psiproxyother[chunk]=symmPCmgr.makeOneResultSection_sym2(chunk);
      }//endfor chunk
  }//endif not gspacesum

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::startNewIter ()  {
//============================================================================
// Check for flow of control errors :

#ifdef _CP_DEBUG_SF_CACHE_
    CkPrintf("GSP [%d,%d] StartNewIter\n",thisIndex.x, thisIndex.y);
#endif
#if CMK_TRACE_ENABLED
    traceUserSuppliedData(iteration);
#endif 
  if(iteration>0){
   if(UegroupProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration_gsp != iteration || 
     UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration  != iteration){
      CkPrintf("{%d} @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n",thisInstance.proxyOffset);
      CkPrintf("{%d} Flow of Control Error : Starting new iter before\n",thisInstance.proxyOffset);
      CkPrintf("{%d} finishing atom integrate or iteration mismatch.\n",thisInstance.proxyOffset);
      CkPrintf("{%d} iter_gsp %d iter_energy %d iter_atm %d and %d and %d\n",
	       thisInstance.proxyOffset,
 	        iteration,UegroupProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration_gsp,
                UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration,
	       config.maxIter,cleanExitCalled);
      CkPrintf("{%d} chare %d %d\n",thisInstance.proxyOffset,thisIndex.x,thisIndex.y);
      CkPrintf("{%d} @@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n",thisInstance.proxyOffset);
      CkExit();
   }//endif
  } //endif

  if(!acceptedVPsi){
      CkPrintf("GSpace[%d,%d] Error: Flow of Control. Starting new iter (%d) before finishing Vpsi.\n",thisIndex.x,thisIndex.y,iteration+1);
      CkAbort("Error: GSpace cannot startNewIter() before finishing the PsiV loop\n");
  }//endif

  // Inform the meshStreamer to be ready for the forward FFT in the first step
  // Setup for subsequent steps happens in doIFFT()
  if (iteration == 0 && config.streamFFTs)
  {
      CkCallback startCb(CkIndex_CP_State_GSpacePlane::readyToStreamFFT(), thisProxy);
      CkCallback endCb(CkCallback::ignore);
      if (thisIndex.x == 0 && thisIndex.y == 0)
          fftStreamer.associateCallback(nstates*nchareG, startCb, endCb, completionDetector, config.rsfftpriority);
  }

  doneNewIter = true;
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  if(thisIndex.x==0 && thisIndex.y==0 ){
    if(!(sim->cp_min_opt==1)){
      fprintf(temperScreenFile,"-------------------------------------------------------------------------------\n");
	  int iii = iteration;
        if(!sim->gen_wave){iii+=1;}
	fprintf(temperScreenFile,"Iteration %d done\n",iii);
	fprintf(temperScreenFile,"===============================================================================\n");
	fprintf(temperScreenFile,"===============================================================================\n");
      }else{
	if(iteration>0){
	  fprintf(temperScreenFile,"===============================================================================\n");
	  fprintf(temperScreenFile,"===============================================================================\n");
	}//endif
        if(iteration<config.maxIter){
  	  fprintf(temperScreenFile,"Beginning Iteration %d \n", iteration);
	}else{
  	  fprintf(temperScreenFile,"Completing Iteration %d \n", iteration-1);
	}//endif
	fprintf(temperScreenFile,"-------------------------------------------------------------------------------\n");
    }//endif
  }//endif
//============================================================================
// Reset all the counters that need to be reset (not more not less)
// otherwise race conditions can leak in.  Rely on the constructor
// for initialization.  Reset set your flags as soon as you are done
// with the tests that require them.

  int cp_min_opt = sim->cp_min_opt;

  // Finished integrate and red psi are safe.
  // You can't get to these until you get passed through this routine
  finishedCpIntegrate = 0;
  iRecvRedPsi      = 1;   if(numRecvRedPsi>0){iRecvRedPsi  = 0;}
  iRecvRedPsiV     = 1;   if(cp_min_opt==0 && numRecvRedPsi>0){iRecvRedPsiV = 0;}
  iSentRedPsi      = 0;
  iSentRedPsiV     = 1;   if(cp_min_opt==0){iSentRedPsiV = 0;}

  iteration++;   // my iteration # : not exactly in sync with other chares
                 //                  but should agree when chares meet.
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("GSP [%d,%d] StartNewIter : %d\n",thisIndex.x, thisIndex.y,iteration);
#endif

//============================================================================
// Check Load Balancing, Increment counter, set done flags equal to false.
/*  if(iteration=3)
    {
      CmiMemoryMark();
    }
  if(iteration=4)
    {
      CmiMemorySweep("GSP");
    }
*/
#if CMK_TRACE_ENABLED
  if(iteration==TRACE_ON_STEP ){(TimeKeeperProxy.ckLocalBranch())->startTrace();}
  if(iteration==TRACE_OFF_STEP){(TimeKeeperProxy.ckLocalBranch())->stopTrace();}
#endif

    if(config.lbgspace || config.lbpaircalc ||config.lbdensity){
	if((iteration % (FIRST_BALANCE_STEP - PRE_BALANCE_STEP) == 0)  || 
           (iteration % (LOAD_BALANCE_STEP - PRE_BALANCE_STEP) == 0)){
	    LBTurnInstrumentOn();
	}//endif
    }//endif

//============================================================================
// Output psi at start of minimization for debugging
    
    if(iteration==1 && cp_min_opt==1){screenOutputPsi(0);}
#ifdef _CP_SUBSTEP_TIMING_
  if(forwardTimeKeep>0){
      double gstart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&gstart,CkReduction::min_double, cb , forwardTimeKeep);
  }//endif
#if USE_HPM
  if(iteration==HPM_ON_STEP ){(TimeKeeperProxy.ckLocalBranch())->startHPM("OneStep");}
  if(iteration==HPM_OFF_STEP ){(TimeKeeperProxy.ckLocalBranch())->stopHPM("OneStep");}
#endif // HPM
#endif // TIMING

//---------------------------------------------------------------------------
    }//end routine

void CP_State_GSpacePlane::screenPrintWallTimes()
{
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  // do new iteration output once globally from the 0th instance
  if(thisIndex.x==0 && thisIndex.y==0 
     && thisInstance.idxU.x==0 && thisInstance.idxU.y==0 
     && thisInstance.idxU.z==0 && thisInstance.idxU.s==0 )
    {
      int iprintout   = config.maxIter;
      if(!(sim->cp_min_opt==1) && !sim->gen_wave){iprintout-=1;}
      int itime       = iteration;
      if(config.maxIter>=30){itime=1; wallTimeArr[0]=wallTimeArr[1];}
      wallTimeArr[itime] = CkWallTimer();
      if (iteration == iprintout && config.maxIter<30) {
	CkPrintf("-------------------------------------------------------------------------------\n");
	CkPrintf("Wall Times from within GSP\n\n");
	for (int t = 1; t < iprintout; t++){
	  CkPrintf("%g\n",wallTimeArr[t] - wallTimeArr[t-1]);
	}//endfor
	if(itime>0)
	  {
	    CkPrintf("%g\n", wallTimeArr[itime] - wallTimeArr[itime-1]);
	    CkPrintf("-------------------------------------------------------------------------------\n");
	  }
      }else{
	if(iteration>0){
	  CkPrintf("Iteration time (GSP) : %g\n", 
		   wallTimeArr[itime] - wallTimeArr[itime-1]);
	}//endif
      }//endif
    }//endif
}

//============================================================================
/**
 * This method is used to start the forward ffts in the CP_State_GSpacePlanes.
 * The work done here could have been done in the acceptData method, but for
 * the fact that we need to synchronize all the CP_State_GSpacePlanes before 
 * doing the forward-ffts. This is a feature, not a bug!
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::doFFT() {

#ifdef _CP_DEBUG_STATEG_VERBOSE_
    CkPrintf("dofft %d.%d \n",thisIndex.x,thisIndex.y);
#endif

  // If there is no data to send, return immediately
  if (gs.numNonZeroPlanes == 0 || gs.fftReqd == false){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude, no data to send : Why launch the stategpsaceplane %d %d %d\n",
               thisIndex.x,thisIndex.y, gs.numNonZeroPlanes);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================

#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif

// Do fft in forward direction, 1-D, in z direction
// A local function not a message : get pointer to memory for fft group

  eesCache *eesData   = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  RunDescriptor *runs = eesData->GspData[iplane_ind]->runs;

  UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->doStpFFTGtoR_Gchare(gs.packedPlaneData,gs.packedForceData, 
	           gs.numFull,gs.numPoints,gs.numLines,gs.numRuns,runs,gs.zdim,iplane_ind);

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(GspaceFwFFT_, StartTime, CmiWallTimer());
#endif   
}
//============================================================================


void CP_State_GSpacePlane::readyToStreamFFT()
{
    isStreamerReady = true;
    if (isForwardFftSendPending)
    {
        sendFFTData();
        isForwardFftSendPending = false;
    }
}

//============================================================================
// Send result to realSpacePlane : perform the transpose
// Force data cannot be overwritten due to all to all nature of comm.
// Until everyone sends, no one gets anything back
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CP_State_GSpacePlane::sendFFTData ()
{

    #ifdef _CP_DEBUG_STATEG_VERBOSE_
        CkPrintf("sendfft %d.%d \n",thisIndex.x,thisIndex.y);
    #endif
    complex *data_out = gs.packedForceData;
    int numLines = gs.numLines; // same amount of data to each realspace chare puppy
    int sizeZ    = gs.planeSize[1];

    //============================================================================
    // Do a Comlib Dance
    #ifdef USE_COMLIB
    #ifdef OLD_COMMLIB
        if (config.useGssInsRealP){gssInstance.beginIteration();}
    #else
        // if (config.useGssInsRealP){ComlibBegin(real_proxy, 0);}
    #endif
    #endif

    #ifdef _CP_SUBSTEP_TIMING_
        if(fftFwdTimer > 0)
        {
            double fftStart = CmiWallTimer();
            CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
            contribute(sizeof(double), &fftStart, CkReduction::min_double, cb, fftFwdTimer);
        }
    #endif
    //============================================================================
    // Send your (x,y,z) to processors z.
    if (config.streamFFTs)
    {
        for(int z=0; z < sizeZ; z++)
        {
            // Hand over the fft data to the MeshStreamer chunk by chunk
            for (int i=0, seq=0; i < numLines; seq++)
            {
                streamedChunk fftchunk(thisIndex.y, thisIndex.x, z, numLines, seq);
                do
                {
                    fftchunk.data[i%streamedChunk::sz] = data_out[z+i*sizeZ];
                } while( (++i < numLines) && (i % streamedChunk::sz != 0) );

                CkArrayIndex2D destIdx(thisIndex.x, z);
                fftStreamer.ckLocalBranch()->insertData(fftchunk, destIdx);
            }
            //CkPrintf("GSpace[%d,%d] sending %d lines in %f chunks to RealSpace[%d,%d]\n",
            //thisIndex.x, thisIndex.y, numLines, std::ceil((double)numLines/streamedChunk::sz),thisIndex.x,z);
        }
        // Now that the streamer is grappling with this load of data,
        fftStreamer.ckLocalBranch()->done();
        // it will become ready only when we set it up for the next iteration
        // This happens in doIFFT(). Until then, its not ready for the next step
        isStreamerReady = false;
    }
    else
    {
        for(int z=0; z < sizeZ; z++)
        {
            // Malloc and prio the message
            RSFFTMsg *msg    = new (numLines,8*sizeof(int)) RSFFTMsg;
            msg->size        = numLines;
            msg->senderIndex = thisIndex.y;  // planenumber
            msg->senderJndex = thisIndex.x;  // statenumber
            msg->senderKndex = z;            // planenumber of rstate
            msg->numPlanes   = gs.numNonZeroPlanes; // unity baby

            if(config.prioFFTMsg)
            {
                CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                *(int*)CkPriorityPtr(msg) = config.rsfftpriority + thisIndex.x*gs.planeSize[0] + thisIndex.y;
            }
            // beam out all points with same z to chare array index z
            complex *data = msg->data;
            for (int i=0,j=z; i<numLines; i++,j+=sizeZ)
                data[i] = data_out[j];
            // send to same state,realspace index [z]
            real_proxy(thisIndex.x, z).acceptFFT(msg);

            // progress engine baby
            CmiNetworkProgress();
        }
    }

    //============================================================================
    // Finish up
    #ifdef USE_COMLIB
    #ifdef OLD_COMMLIB
        if (config.useGssInsRealP){gssInstance.endIteration();}
    #else
        //if (config.useGssInsRealP){ComlibEnd(real_proxy, 0);}
    #endif
    #endif

    #ifdef _CP_SUBSTEP_TIMING_
        if(forwardTimeKeep>0)
        {
            double gend = CmiWallTimer();
            CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
            contribute(sizeof(double),&gend,CkReduction::max_double, cb , forwardTimeKeep);
        }
    #endif
}



/**
 * Force cannot be overwitten because we can't receive until all chares send. The beauty of all to all comm
 * Forces are initialized HERE. No need to zero them etc. elsewhere.
 */
void CP_State_GSpacePlane::acceptIFFT(GSIFFTMsg *msg) 
{
#ifdef _CP_SUBSTEP_TIMING_
  if(backwardTimeKeep>0)
    {
      double gstart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&gstart,CkReduction::min_double, cb , backwardTimeKeep);
    }
#endif

#ifdef _CP_DEBUG_STATEG_VERBOSE_
	CkPrintf("GSpace[%d,%d] acceptIFFT\n",thisIndex.x,thisIndex.y);
#endif
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
      CkAssert(finite(msg->data[i].re));
      CkAssert(finite(msg->data[i].im));
    }
#endif

  int size             = msg->size;
  int offset           = msg->offset;
  complex *partlyIFFTd = msg->data;

  complex *data_in     = gs.packedForceData;
  int numLines         = gs.numLines;
  int sizeZ            = gs.planeSize[1];
  int expandedDataSize = numLines*sizeZ;

  CkAssert(numLines == size);

//============================================================================
// Recv the message

  countIFFT++;

  // This is not a reduction. Don't zero me please. Every elements is set.
  // z=offset is inner index : collections of z-lines of constant (gx,gy)
  for(int i=0,j=offset; i< numLines; i++,j+=sizeZ){data_in[j] = partlyIFFTd[i];}

  delete msg;

//============================================================================
// If you have recved from every z plane, go on

  if (countIFFT == gs.planeSize[1]) {
    #ifdef _CP_SUBSTEP_TIMING_
        if(fftBwdTimer > 0)
        {
            double fftEnd = CmiWallTimer();
            CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
            contribute(sizeof(double), &fftEnd, CkReduction::max_double, cb , fftBwdTimer);
        }
    #endif
    countIFFT = 0;
        UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl();
  }//endif : has everyone arrived?

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::doIFFT() 
{
    // Now do the IFFT in place
    #if CMK_TRACE_ENABLED
        double StartTime=CmiWallTimer();
    #endif

    if (config.streamFFTs && thisIndex.x == 0 && thisIndex.y == 0)
    {
        // Inform the meshStreamer to be ready for the forward FFT in the next step
        CkCallback startCb(CkIndex_CP_State_GSpacePlane::readyToStreamFFT(), thisProxy);
        CkCallback endCb(CkCallback::ignore);
        fftStreamer.associateCallback(nstates*nchareG, startCb, endCb, completionDetector, config.rsfftpriority);
    }

    eesCache *eesData   = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    RunDescriptor *runs = eesData->GspData[iplane_ind]->runs;
    FFTcache *fftcache        = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    
    fftcache->getCacheMem("CP_State_GSpacePlane::doIFFT");
    complex *forcTmp = fftcache->tmpData;
    fftcache->doStpFFTRtoG_Gchare(gs.packedForceData,forcTmp,
    gs.numFull,gs.numPoints,gs.numLines,gs.numRuns,runs,gs.zdim,iplane_ind);
    
    #if CMK_TRACE_ENABLED
        traceUserBracketEvent(GspaceBwFFT_, StartTime, CmiWallTimer());
    #endif
    
    // Finish up by multiplying the the FFT scaling factor
    CPcharmParaInfo *sim = CPcharmParaInfo::get();
    double scaleFactor = -2.0/double(sim->sizeX * 
                                     sim->sizeY * 
                                     sim->sizeZ);
    complex *forces = gs.packedForceData;
    for(int index=0; index<gs.numPoints; index++)
        forces[index] = forcTmp[index]*scaleFactor;
    fftcache->freeCacheMem("CP_State_GSpacePlane::doIFFT");
    CmiNetworkProgress();
    
    // Report that you are done
    #ifdef _CP_DEBUG_STATEG_VERBOSE_
        if(thisIndex.x==0)
            CkPrintf("Done doing Ifft gsp : %d %d\n",thisIndex.x,thisIndex.y);
    #endif

    /// If there is a barrier after the IFFTs, contribute to the reduction barrier that will sync all GSpace chares
    #ifdef BARRIER_CP_GSPACE_IFFT
        //put contribute here to reduction with a broadcast client
        int wehaveours=1;
        contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
        CkCallback(CkIndex_GSpaceDriver::allDoneIFFT(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
    #else
        doneDoingIFFT = true;
    #endif
}//end routine
//=============================================================================================================



//=============================================================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================================================
void CP_State_GSpacePlane::combineForcesGetEke(){

#ifdef _CP_DEBUG_VKS_OFF_
  if(thisIndex.x==0 && thisIndex.y==0){
    CkPrintf("EHART       = OFF FOR DEBUGGING\n");
    CkPrintf("EExt        = OFF FOR DEBUGGING\n");
    CkPrintf("EWALD_recip = OFF FOR DEBUGGING\n");
    CkPrintf("EEXC        = OFF FOR DEBUGGING\n");
    CkPrintf("EGGA        = OFF FOR DEBUGGING\n");
    CkPrintf("EEXC+EGGA   = OFF FOR DEBUGGING\n");
  }//endif
#endif

//================================================================================
// Check stuff from the gsplane data cache and particle plane

  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int ncoef         = gs.numPoints;
  int *k_x          = eesData->GspData[iplane_ind]->ka;
  int *k_y          = eesData->GspData[iplane_ind]->kb;
  int *k_z          = eesData->GspData[iplane_ind]->kc;
  double **g2       = eesData->GspData[iplane_ind]->g2;

//================================================================================
// Add forces from particle plane to forces from IFFT then zero them

#ifdef _CP_DEBUG_VKS_FORCES_
  if(thisIndex.x==0 && thisIndex.y==0){
    FILE *fp = fopen("vks_forces_s0_p0.out","w");
    int ncoef       = gs.numPoints;
    complex *forces = gs.packedForceData;
    for(int i=0;i<ncoef;i++)
      fprintf(fp,"%d %d %d : %g %g\n",k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
    fclose(fp);
  }//endif
#endif

  complex *forces = gs.packedForceData;
  complex *ppForces = UparticlePlaneProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).ckLocal()->myForces; 
  
#ifdef _CP_DEBUG_SFNL_OFF_
  bzero(ppForces,ncoef*sizeof(complex));
#endif
#ifdef _CP_DEBUG_VKS_OFF_ 
  bzero(forces,ncoef*sizeof(complex));
  bzero(ppForces,ncoef*sizeof(complex));
#endif


#ifdef _CP_GS_DUMP_VKS_
    dumpMatrix("vksBf",(double *)ppForces, 1, 
                gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    dumpMatrix("forceBf",(double *)forces, 1, 
                gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  if(savedvksBf==NULL){ // load it
      savedvksBf= new complex[gs.numPoints];
      loadMatrix("vksBf",(double *)savedvksBf, 1, 
                 gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  if(savedforceBf==NULL){ // load it
      savedforceBf= new complex[gs.numPoints];
      loadMatrix("forceBf",(double *)savedforceBf, 1, 
                gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif

  for(int i=0;i<gs.numPoints;i++){
      if(fabs(ppForces[i].re-savedvksBf[i].re)>0.0001){
	 fprintf(stderr, "GSP [%d,%d] %d element vks  %.10g not %.10g\n",
              thisIndex.x, thisIndex.y,i, ppForces[i].re, savedvksBf[i].re);
      }//endif
      CkAssert(fabs(ppForces[i].re-savedvksBf[i].re)<0.0001);
      CkAssert(fabs(ppForces[i].im-savedvksBf[i].im)<0.0001);
  }//endfor
  for(int i=0;i<gs.numPoints;i++){
      if(fabs(forces[i].re-savedforceBf[i].re)>0.0001){
	  fprintf(stderr, "GSP [%d,%d] %d element force  %.10g not %.10g\n",
              thisIndex.x, thisIndex.y,i, forces[i].re, savedforceBf[i].re);
      }//endif
      CkAssert(fabs(forces[i].re-savedforceBf[i].re)<0.0001);
      CkAssert(fabs(forces[i].im-savedforceBf[i].im)<0.0001);
  }//endfor
#endif

  gs.addForces(ppForces,k_x);
  bzero(ppForces,gs.numPoints*sizeof(complex));
  CmiNetworkProgress();

//================================================================================
// Compute force due to quantum kinetic energy and add it in.
// Reduce quantum kinetic energy or eke

  int istate        = gs.istate_ind;
  int nkx0          = gs.nkx0;
  complex *psi_g    = gs.packedPlaneData;
  double *eke_ret   = &(gs.eke_ret);

  CPNONLOCAL::CP_eke_calc(ncoef,istate,forces,psi_g,k_x,k_y,k_z,g2,eke_ret,config.doublePack,nkx0,
                          kpoint_ind,config.nfreq_cpnonlocal_eke);
  contribute(sizeof(double), &gs.eke_ret, CkReduction::sum_double, 
	     CkCallback(CkIndex_InstanceController::printEnergyEke(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
  //isEnergyReductionDone = false; ///@note: This doesnt seem necessary here and commenting out has not affected simple tests. This flag is reset at the start of the iter itself.

#ifdef _CP_DEBUG_SCALC_ONLY_ 
  bzero(forces,ncoef*sizeof(complex));
#endif

//========================================================================
// Debug output

#ifdef _CP_DEBUG_OLDFORCE_
  if(ncoef >0){
    FILE *fp = fopen("force_old.out", "a+");
    for(int i = 0; i < ncoef; i++){
      if(k_x[i]==0 && k_y[i]==1 && k_z[i]==4){
	fprintf(fp,
		"old force H+Ext+Exc+Eke+Enl : is=%d %d %d %d : %g %g\n",
		istate,k_x[i],k_y[i],k_z[i],forces[i].re,forces[i].im);
      }
    }
    fclose(fp);
  }//endif
#endif

//-----------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
// Atoms are launched by allDoneCpForces when all after ALL planes and states 
// have reported.
//==============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::launchAtoms() {
  //  CkPrintf("{%d} GSP [%d,%d] launchAtoms\n",thisInstance.proxyOffset, thisIndex.x,thisIndex.y);

//==============================================================================
// The usual stuff
#ifdef _CP_DEBUG_PSI_OFF_
//  iteration++;
  if(iteration==config.maxIter+1){
    int i=0;
#ifdef _CP_SUBSTEP_TIMING_
#if USE_HPM
    (TimeKeeperProxy.ckLocalBranch())->printHPM();
#endif	
#endif
    cleanExitCalled = 1;

    contribute(sizeof(int),&cleanExitCalled,CkReduction::sum_int,  CkCallback(CkIndex_InstanceController::cleanExit(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
  }else{
#endif
   int i=0;
   contribute(sizeof(int),&i,CkReduction::sum_int,  CkCallback(CkIndex_InstanceController::allDoneCPForces(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));
#ifdef _CP_DEBUG_PSI_OFF_
  }//endif
#endif
}//end routine



//==============================================================================
// After MY Cp forces have arrived : sendLambda
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendLambda() {
//==============================================================================
// Reset set lambda (not done) and force counters (not done for NEXT step):

   acceptedLambda        = false;
   doneDoingIFFT         = false;
   doneNewIter           = false;

//==============================================================================
// Get nice local variables and if dynamics make a copy

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  complex *psi   = gs.packedPlaneData;
  complex *force = gs.packedForceData;
  int cp_min_opt = sim->cp_min_opt;

  if(cp_min_opt==0){ // dynamics in on : minimization is off
    int ncoef           = gs.numPoints;
    complex *psi_g      = gs.packedPlaneData;
    complex *psi_g_tmp  = gs.packedPlaneDataTemp;
    CmiMemcpy(psi_g_tmp,psi_g,sizeof(complex)*ncoef);
  }//endif

//==============================================================================
// Debug output schmoo : contrib to diagonal element of lambda from this chare
//                     : Most effective for 1 gspace chare when you
//                     : get the diangonal elemens of lambda (see below)


#define _CP_GSPACE_DUMP_LMAT_DIAGONAL_VALS_OFF_

#ifdef _CP_GSPACE_DUMP_LMAT_DIAGONAL_VALS_
  /* The lambda matrix is the reduced result of the forward path GEMMs in the
   * asymmetric paircalcs that arrives at Ortho::aceptSectionLambda.  There is
   * an LMAT dump macro over there. Here, we manually compute the product of
   * psi and the forces to obtain the contribution of this chare to the
   * diagonal element of the L matrix.  To verify that the forward path is
   * producing correct results, we can compare these hand-computed diagonal
   * elements with the lambda matrix dumped in Ortho.
   *
   * Running the code with just one plane (reduce gExpandFact to get just one
   * plane) should ensure that there is just one GSpace chare for each state
   * and that the diagonal elements are computed in this chare itself.
   * Otherwise, the diagonal elements have to be computed by summing the
   * corresponding values dumped by all the GSpace chares in a state (ie across
   * all planes).
   */

    complex mylambda_diag = 0;
    for(int i=0; i<gs.numPoints; i++){
        mylambda_diag += force[i] * psi[i].conj();
    }/*endfor*/
    if(config.doublePack==0){
      CkPrintf("lambda[%d, %d] = %.12g %.12g at plane %d\n",
          thisIndex.x, thisIndex.x, mylambda_diag.re, mylambda_diag.im, thisIndex.y);
    }else{
      CkPrintf("lambda[%d, %d] = %.12g %.12g at plane %d\n",
          thisIndex.x, thisIndex.x, mylambda_diag.re,0.0, thisIndex.y);
    }/*endif*/
#endif

//==============================================================================
// Debug output schmoo : output forces before lambda-ization!
//                     : 


#define _CP_GSPACE_PSI_FORCE_OUTPUT_BEFORE_LAMBDA_OFF_
#ifdef  _CP_GSPACE_PSI_FORCE_OUTPUT_BEFORE_LAMBDA_
    eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    int *ka           = eesData->GspData[iplane_ind]->ka;
    int *kb           = eesData->GspData[iplane_ind]->kb;
    int *kc           = eesData->GspData[iplane_ind]->kc;

    char fname[1024];
    sprintf(fname,"psi_forces_before_lambda_state%d_plane%d.out",thisIndex.x,thisIndex.y);
    FILE *fp = cfopen(fname,"w");
    for(int i=0; i<gs.numPoints; i++){
        int igo = 0; double wght = 0.5;
        if (ka[i] > 0) igo=1;
        if (ka[i] == 0 && kb[i]>0){ igo=1; wght=1.0; }
        if (ka[i] == 0 && kb[i]==0 && kc[i]>=0){ igo=1; wght=1.0; }
        if(config.doublePack==0) wght=1.0;
        if(igo==1)
          fprintf(fp,"%d %d %d : %.10g %.10g\n",ka[i],kb[i],kc[i],force[i].re*wght,force[i].im*wght);
    }/*endfor*/
    fclose(fp);
#endif

//==============================================================================
// Scale the variables for dynamics and double packing : Single pack
// will require modification of factors of 2 in PC or scaling all the
// variables which stinks

#ifndef PAIRCALC_TEST_DUMP
#ifndef _CP_DEBUG_ORTHO_OFF_
  if(cp_min_opt==0){
    double rad2i = 1.0/sqrt(2.0);
    double rad2  = sqrt(2.0);
    if(gs.ihave_kx0==1 && config.doublePack==1){
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i]   *= rad2i;}
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){force[i] *= rad2;}
    }else{
      for(int i=0; i<gs.numPoints; i++){
        psi[i]   *= rad2i; 
        force[i] *= rad2;
      }//endfor
    }//endif
  }//endif
#endif
#endif

//==============================================================================
// Enormous debug schmoo

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  for(int i=0;i<config.numChunksAsym;i++)
    CkAssert(countLambdaO[i]==0);
#endif
#ifdef _CP_GS_DUMP_LAMBDA_
    dumpMatrix("lambdaBf",(double *)force, 1, 
                     gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    dumpMatrix("psiBf",(double *)psi, 1, 
                     gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  if(savedlambdaBf==NULL){ // load it
      savedlambdaBf= new complex[gs.numPoints];
      loadMatrix("lambdaBf",(double *)savedlambdaBf, 1, 
                       gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  if(savedpsiBf==NULL){ // load it
      savedpsiBf= new complex[gs.numPoints];
      loadMatrix("psiBf",(double *)savedpsiBf, 1, 
                        gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  double testvalue=0.00000001;
  for(int i=0;i<gs.numPoints;i++){
      if(fabs(force[i].re-savedlambdaBf[i].re)>testvalue){
	  fprintf(stderr, "GSP [%d,%d] %d element lambda  %.10g not %.10g\n",
                  thisIndex.x, thisIndex.y,i, force[i].re, savedlambdaBf[i].re);
      }//endif
      CkAssert(fabs(force[i].re-savedlambdaBf[i].re)<testvalue);
      CkAssert(fabs(force[i].im-savedlambdaBf[i].im)<testvalue);
  }//endfor
  for(int i=0;i<gs.numPoints;i++){
      CkAssert(fabs(psi[i].re-savedpsiBf[i].re)<testvalue);
      CkAssert(fabs(psi[i].im-savedpsiBf[i].im)<testvalue);
  }//endfor
#endif

//==============================================================================
// Send to lambda : Finally

  int numPoints   = gs.numPoints;
#ifndef _CP_DEBUG_ORTHO_OFF_
  asymmPCmgr.sendLeftData(numPoints,psi,false);
  CmiNetworkProgress();
  asymmPCmgr.sendRightData(numPoints,force,false);
#else
  acceptedLambda=true;
  bzero(force,sizeof(complex)*numPoints);
#endif
#ifdef _CP_DEBUG_STATEG_VERBOSE_ 
  if(thisIndex.x==0){CkPrintf("Sent Lambda %d %d\n",thisIndex.y,cleanExitCalled);}
#endif

//==============================================================================
// Eric's micro-timing 

#ifdef _CP_SUBSTEP_TIMING_
  if(backwardTimeKeep>0){
      double gend=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&gend,CkReduction::max_double, cb , backwardTimeKeep);
  }//endif
#endif

//-----------------------------------------------------------------------------
   }// end routine 
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptLambda(CkReductionMsg *msg) {
//==============================================================================
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt    = sim->cp_min_opt;
  int cp_lsda       = sim->nspin - 1;  // 1 for lsda and 0 for lda
  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind]->ka;

  int       N    = msg->getSize()/sizeof(complex);
  complex *data  = (complex *)msg->getData();
  int offset     = msg->getUserFlag();  if(offset<0){offset=0;}
  //  CkPrintf("[%d %d] accepts lambda %d \n", thisIndex.x, thisIndex.y,offset);
  complex *force = gs.packedForceData;
  int chunksize  = gs.numPoints/config.numChunksAsym;
  int chunkoffset=offset*chunksize; // how far into the points this contribution lies

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->getSize()/sizeof(double) ;i++){
      CkAssert(finite(((double*) msg->getData())[i]));
  }//endfor
#endif

//==============================================================================
// Count then debug

  countLambda++;//lambda arrives in as many as 2 * numChunks reductions

/*
  if(thisIndex.y==0){
    dumpMatrix("lambdab4",(double *)force, 1, gs.numPoints*2,
                      thisIndex.y,thisIndex.x,thisIndex.x,0,false);
  }//endif

  if(thisIndex.y==0){
      CkPrintf("LAMBDA [%d %d], offset %d chunkoffset %d N %d countLambdao %d\n", 
                thisIndex.x, thisIndex.y, offset, chunkoffset, N, countLambdaO[offset]);
  }//endif
*/

//==============================================================================
// Add the forces 

  //---------------------------------------------------
  // A) BGL STuff
#ifdef CMK_BLUEGENEL
#pragma disjoint(*force, *data)
  //      __alignx(16,force);
  //      __alignx(16,data);
#endif

  //---------------------------------------------------
  // B) Double Pack 
  if(config.doublePack==1){
   if(cp_min_opt==1){
     double overlap = (cp_lsda==0 ? 2.0 : 1.0);
     double ws = 1.0/overlap; double wd = 2.0/overlap;
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
#ifndef PAIRCALC_TEST_DUMP
     for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       double wght  = (k_x[idest]==0 ? ws : wd);
       force[idest].re -= wght*data[i].re;
       force[idest].im -= wght*data[i].im;
     }//endfor
#endif
   }else{
     if(countLambdaO[offset]<1){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
     }else{
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
     }// coutlambda
   }//endif : cpmin
 
  }//endif : double pack

  //---------------------------------------------------
  // C) Double Pack is off
  if(config.doublePack==0){
   if(cp_min_opt==1){
     double overlap = (cp_lsda==0 ? 2.0 : 1.0);
     double ws = 1.0/overlap;
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
    for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       force[idest].re -= ws*data[i].re;
       force[idest].im -= ws*data[i].im;
    }//endfor
   }else{
     if(countLambdaO[offset]<1){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
     }else{
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
     }// endif : 1st guy
   }//endif : minimization
  }//endif : singlePack = kpts

  //---------------------------------------------------
  // D) Cleanup
  delete msg;  

//==============================================================================
// Do we have everything?

  countLambdaO[offset]++;
  if(countLambda==AllLambdaExpected){ 
    thisProxy(thisIndex.x,thisIndex.y).doLambda();
#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("doLambda %d %d\n",thisIndex.y,cleanExitCalled);
#endif
  }//endif

//==============================================================================
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptLambda(partialResultMsg *msg) {
//==============================================================================
// 0) unpack the message : pop out variables from groups

  complex *data     = msg->result;
  int       N       = msg->N;
  int offset        = msg->myoffset;
  if(offset<0){offset=0;}

#ifdef _NAN_CHECK_
  for(int i=0; i < N; i++)
  {
      CkAssert(finite(data[i].re));
      CkAssert(finite(data[i].im));
  }
#endif

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt    = sim->cp_min_opt;
  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind]->ka;

  complex *force    = gs.packedForceData;
  int chunksize     = gs.numPoints/config.numChunksAsym;
  int chunkoffset   = offset*chunksize; // how far into the points this contribution lies

//==============================================================================
// I) Increment the counter and do some checking

  countLambda++;  //lambda arrives in as many as 2 * numChunks reductions

  //  CkPrintf("GSP [%d,%d] acceptLambda N %d offset %d countLambda %d expecting %d\n",thisIndex.x,thisIndex.y,N,offset,countLambda,AllLambdaExpected);

//=============================================================================
// (II) Add it in to our forces : Careful about offsets, doublepack and cpmin/cp

 //----------------------------------------------------------
 //A) BlueGene nonsense

#ifdef CMK_BLUEGENEL
#pragma disjoint(*force, *data)
  //      __alignx(16,force);
  //      __alignx(16,data);
#endif
 //----------------------------------------------------------
 //B) Double Pack

  if(config.doublePack==1){
   if(cp_min_opt==1){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
#ifndef PAIRCALC_TEST_DUMP
     for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       double wght  = (k_x[idest]==0 ? 0.5 : 1);
       force[idest].re -= wght*data[i].re;
       force[idest].im -= wght*data[i].im;
     }//endfor
#endif
   }else{
     if(countLambdaO[offset]<1){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
        for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
     }else{
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
        for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
     }//endif : off set thingy
   }//endif : cp_min_on
  }//endif : double pack

 //----------------------------------------------------------
 //C) Single pack

  if(config.doublePack==0){
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
#endif
    for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       force[idest].re -= 0.5*data[i].re;
       force[idest].im -= 0.5*data[i].im;
    }//endfor

  }//endif : single pack

 //----------------------------------------------------------
 //D) Clean up

  delete msg;  

//=============================================================================
// are we done?

  countLambdaO[offset]++;
  if(countLambda==AllLambdaExpected){ 
    thisProxy(thisIndex.x,thisIndex.y).doLambda();
#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("doLambda %d %d\n",thisIndex.y,cleanExitCalled);
#endif
  }//endif

//==============================================================================
    }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doLambda() {
//==============================================================================
// (I) If you have got it all : Rescale it and resume
  // CkPrintf("{%d} GSP [%d,%d] doLambda\n",thisInstance.proxyOffset, thisIndex.x,thisIndex.y);
  CkAssert(countLambda==AllLambdaExpected);
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt = sim->cp_min_opt;
  complex *force = gs.packedForceData;
#ifdef _NAN_CHECK_
  for(int i=0;i<gs.numPoints ;i++){
      CkAssert(finite(force[i].re));
      CkAssert(finite(force[i].im));
    }
#endif

  if(cp_min_opt==0){
    // dynamics scale it out
    if(gs.ihave_kx0==1){
      double rad2i = 1.0/sqrt(2.0);
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){force[i] *= rad2i;}
    }//endif
  }//endif

//==============================================================================
// Retrieve Non-orthog psi

  if(cp_min_opt==0){
    int ncoef          = gs.numPoints;
    complex *psi_g_scr = gs.packedPlaneDataScr;
    complex *psi_g     = gs.packedPlaneData;
    CmiMemcpy(psi_g,psi_g_scr,sizeof(complex)*ncoef); // overwrite ortho with non-ortho
  }//endif
    
//==============================================================================
// Resume 

  acceptedLambda=true;
  countLambda=0;
  bzero(countLambdaO,config.numChunksAsym*sizeof(int));

//==============================================================================
// Debug Schmoo

#ifdef _CP_GS_DUMP_LAMBDA_
    dumpMatrix("lambdaAf",(double *)force, 1, gs.numPoints*2,
                      thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  double testvalue=0.00000001;
  if(savedlambdaAf==NULL){ // load it
      savedlambdaAf= new complex[gs.numPoints];
      loadMatrix("lambdaAf",(double *)savedlambdaAf, 1, 
                        gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  for(int i=0;i<gs.numPoints;i++){
      CkAssert(fabs(force[i].re-savedlambdaAf[i].re)<testvalue);
      CkAssert(fabs(force[i].im-savedlambdaAf[i].im)<testvalue);
  }//endfor
#endif

#ifdef _CP_DEBUG_STATEG_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("doLambda %d %d\n",thisIndex.y,cleanExitCalled);
#endif

//==============================================================================
// Debug : Write out forces after Lambda-ization


#define _CP_GSPACE_PSI_FORCE_OUTPUT_AFTER_LAMBDA_OFF_
#ifdef  _CP_GSPACE_PSI_FORCE_OUTPUT_AFTER_LAMBDA_
    eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    int *ka           = eesData->GspData[iplane_ind]->ka;
    int *kb           = eesData->GspData[iplane_ind]->kb;
    int *kc           = eesData->GspData[iplane_ind]->kc;

    char fname[1024];
    sprintf(fname,"psi_forces_after_lambda_state%d_plane%d.out",thisIndex.x,thisIndex.y);
    FILE *fp = cfopen(fname,"w");
    for(int i=0; i<gs.numPoints; i++){
        int igo = 0; double wght = 0.5;
        if (ka[i] > 0) igo=1;
        if (ka[i] == 0 && kb[i]>0){ igo=1; wght=1.0; }
        if (ka[i] == 0 && kb[i]==0 && kc[i]>=0){ igo=1; wght=1.0; }
        if(config.doublePack==0) wght=1.0;
        if(igo==1)
          fprintf(fp,"%d %d %d : %.10g %.10g\n",ka[i],kb[i],kc[i],force[i].re*wght,force[i].im*wght);
    }/*endfor*/
    fclose(fp);
#endif

//==============================================================================
// Back to the threaded loop

  UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl();

//==============================================================================
  }//end routine
//==============================================================================

//=========================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=========================================================================
// If minimization :
//      compute lambda matrix (Lagrange multipliers)
//      modify forces using the lambda matrix
//      compute the |force|^2 = fovlap
//=========================================================================
void CP_State_GSpacePlane::computeCgOverlap() {
//=========================================================================
// Flow of control check and local variables
//  CkPrintf("{%d} GSP [%d,%d] computeCgOverlap\n",thisInstance.proxyOffset, thisIndex.x,thisIndex.y);
  if(!acceptedLambda){
     CkPrintf("GSpace[%d,%d] Flow of Control Error : Attempting to Cg ovlap without lambda correction\n",thisIndex.x,thisIndex.y);
     CkAbort("Error: Attempting to compute cg overlap without lambda correction\n");
   }//endif

   CPcharmParaInfo *sim = CPcharmParaInfo::get();
   int istate      = gs.istate_ind;
   int ncoef       = gs.numPoints;
   complex *forces = gs.packedForceData;
   int cp_min_opt  = sim->cp_min_opt;
   int cp_min_cg   = sim->cp_min_cg;

//=========================================================================

   double fovlap_loc = 0.0;
   if(cp_min_opt==1 && cp_min_cg==1){
     CPINTEGRATE::CP_fovlap_calc(ncoef,istate,forces,&fovlap_loc);
   }//endif
   double force_sq_sum_loc=0.0;
   for(int i=0; i<gs.numPoints; i++){force_sq_sum_loc+= forces[i].getMagSqr();}

   double redforc[2];
   redforc[0] = fovlap_loc;
   redforc[1] = force_sq_sum_loc;
   contribute(2*sizeof(double),redforc,CkReduction::sum_double,
	     CkCallback(CkIndex_CP_State_GSpacePlane::psiCgOvlap(NULL),thisProxy));

#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("ovlpa %d %d\n",thisIndex.y,cleanExitCalled);
#endif

//----------------------------------------------------------------------------
  }// end routine : computeCgOverlap
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// In this function data is written to files the simpliest way possible
//============================================================================
void CP_State_GSpacePlane::writeStateDumpFile()
{
  // Local pointers, variables and error checking
  if(!acceptedPsi || !acceptedLambda || !acceptedVPsi)
    CkAbort("{%d} Flow of Control Error : Attempting to write states without completing psi, vpsi and Lambda\n");


  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt = (sim->cp_min_opt);

  int ind_state  = (thisIndex.x+1);
  int ind_chare  = (thisIndex.y+1);
  int ncoef      = gSpaceNumPoints;

  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind]->ka;
  int *k_y          = eesData->GspData[iplane_ind]->kb;
  int *k_z          = eesData->GspData[iplane_ind]->kc;
  double *coef_mass = eesData->GspData[iplane_ind]->coef_mass;

  complex *psi      = gs.packedPlaneData;
  complex *vpsi     = gs.packedVelData;       // to evolve psi to time, t.
  complex *forces   = gs.packedForceData;
  double ***xNHC    = gs.xNHC;
  double ***xNHCP   = gs.xNHCP;
  double ***vNHC    = gs.vNHC;
  double ***fNHC    = gs.fNHC;          
  double *mNHC      = gs.mNHC;          
  int len_nhc_cp    = gs.len_nhc_cp;
  int num_nhc_cp    = gs.num_nhc_cp;
  int nck_nhc_cp    = gs.nck_nhc_cp;
  int nkx0_red      = gs.nkx0_red;
  int nkx0_uni      = gs.nkx0_uni;
  int nkx0_zero     = gs.nkx0_zero;
  double kTCP       = gs.kTCP;

  int myiteration = iteration;
  if(cp_min_opt==0){myiteration=iteration-1;}

  //============================================================================
  // Set the file names and write the files
  if(ind_state==1 && ind_chare==1){
    CkPrintf("-----------------------------------\n");
    CkPrintf("Writing states to disk at time %d\n",myiteration);
    CkPrintf("-----------------------------------\n");
  }//endif
  //------------------------------------------------------------------
  // Update the velocities into scratch as we are between steps
  if(cp_min_opt==0 && halfStepEvolve ==1){
    halfStepEvolve = 0;
    CPINTEGRATE::cp_evolve_vel(ncoef,forces,vpsi,coef_mass,
			       len_nhc_cp,num_nhc_cp,nck_nhc_cp,fNHC,vNHC,xNHC,xNHCP,mNHC,
			       gs.v0NHC,gs.a2NHC,gs.a4NHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,
			       2,iteration,gs.degfree,gs.degfreeNHC,gs.degFreeSplt,
			       gs.istrNHC,gs.iendNHC,1);
  }//endif
  //------------------------------------------------------------------
  // Pack the message and send it to your plane 0
  GStateOutMsg *msg  = new (ncoef,ncoef,ncoef,ncoef,ncoef,
			    8*sizeof(int)) GStateOutMsg;
  if(config.prioFFTMsg){
    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
    *(int*)CkPriorityPtr(msg) = config.rsfftpriority + 
      thisIndex.x*gs.planeSize[0]+thisIndex.y;
  }//endif
  msg->size        = ncoef;
  msg->senderIndex = thisIndex.y;  // planenumber
  complex *data    = msg->data;
  complex *vdata   = msg->vdata; 
  int *mk_x        = msg->k_x;
  int *mk_y        = msg->k_y;
  int *mk_z        = msg->k_z;
  if(cp_min_opt==0){
    for(int i=0;i<ncoef;i++){vdata[i] = vpsi[i];}
  }else{
    for(int i=0;i<ncoef;i++){vdata[i] = 0.0;}
  }//endif
  for (int i=0;i<ncoef; i++){
    data[i]  = psi[i];  
    mk_x[i]  = k_x[i];  mk_y[i]  = k_y[i];  mk_z[i]  = k_z[i];
  }//endfor
  UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x,redPlane).collectFileOutput(msg);
  //------------------------------------------------------------------
  // If you are not plane redPlane, you are done. Invoke the correct reduction.
  if(thisIndex.y!=redPlane){
    if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==1)
      UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).readyToExit();
    else
      {
	int i = 0;
	contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(CkIndex_GSpaceDriver::allDoneWritingPsi(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
      }
  }//endif
}//write the file
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::collectFileOutput(GStateOutMsg *msg){
//============================================================================

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int sizeX        = (sim->sizeX);
  int sizeY        = (sim->sizeY);
  int sizeZ        = (sim->sizeZ);
  int npts_tot     = (sim->npts_tot);
  int *ipacked_off = (sim->index_output_off);
  int cp_min_opt   = (sim->cp_min_opt);
  int ibinary_write_opt = sim->ibinary_write_opt;

  int myiteration = iteration;
  if(cp_min_opt==0){myiteration=iteration-1;}

//============================================================================
// Receive the message

  if(countFileOut==0){
    tpsi  = (complex *)fftw_malloc(npts_tot*sizeof(complex));
    tvpsi = (complex *)fftw_malloc(npts_tot*sizeof(complex));
    tk_x  = (int *)fftw_malloc(npts_tot*sizeof(int));
    tk_y  = (int *)fftw_malloc(npts_tot*sizeof(int));
    tk_z  = (int *)fftw_malloc(npts_tot*sizeof(int));
  }//endif
  countFileOut++;

  int ncoef      = msg->size;
  int myplane    = msg->senderIndex;
  complex *data  = msg->data;
  complex *vdata = msg->vdata; 
  int *mk_x      = msg->k_x;
  int *mk_y      = msg->k_y;
  int *mk_z      = msg->k_z;
  int ioff       = ipacked_off[myplane];
  for(int i=0,j=ioff;i<ncoef;i++,j++){
    tpsi[j]  = data[i];
    tvpsi[j] = vdata[i];
    tk_x[j]  = mk_x[i];
    tk_y[j]  = mk_y[i];
    tk_z[j]  = mk_z[i];
  }//endfor

  delete msg;

//============================================================================
// If you've got the whole state, write it out and then invoke the reduction.

  if(countFileOut==nchareG){
     countFileOut = 0;
     int ind_state = thisIndex.x+1;
     int ibead     = myBeadIndex;
     int ikpt      = myKptIndex;
     int itemper   = myTemperIndex;
     int ispin     = mySpinIndex;

     char psiName[400]; char vpsiName[400];

     sprintf(psiName,  "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/state%d.out",
                       config.dataPathOut,ispin,ikpt,ibead,itemper,ind_state);
     sprintf(vpsiName, "%s/Spin.%d_Kpt.%d_Bead.%d_Temper.%d/vState%d.out",
                       config.dataPathOut,ispin,ikpt,ibead,itemper,ind_state);
     writeStateFile(npts_tot,tpsi,tvpsi,tk_x,tk_y,tk_z,cp_min_opt,
                    sizeX,sizeY,sizeZ,psiName,vpsiName,ibinary_write_opt,
                    myiteration,ind_state,ispin,ikpt,ibead,itemper);
     fftw_free(tpsi); tpsi  = NULL;
     fftw_free(tvpsi);tvpsi = NULL;
     fftw_free(tk_x); tk_x  = NULL;
     fftw_free(tk_y); tk_y  = NULL;
     fftw_free(tk_z); tk_z  = NULL;
	if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==1)
	  UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).readyToExit();
	else{
   	  int i=0;
	  contribute(sizeof(int),&i,CkReduction::sum_int,
          CkCallback(CkIndex_GSpaceDriver::allDoneWritingPsi(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
	}//endif
    }//endif

//============================================================================
  }//end routine
//============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::integrateModForce() {
//==============================================================================
// I. Error Conditions : you haven't accepted present steps lambda
//                       you haven't accepted previous steps psi or vpsi

  if(!acceptedLambda || !acceptedVPsi || !acceptedPsi){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Flow of Control Error : Attempting to cp integrate\n");
    CkPrintf("without lambda correction to forces or psi or vpsi\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("integrate %d %d\n",thisIndex.y,cleanExitCalled);
#endif

//==============================================================================
// II. Local pointers    

  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind]->ka;
  int *k_y          = eesData->GspData[iplane_ind]->kb;
  int *k_z          = eesData->GspData[iplane_ind]->kc;
  double *coef_mass = eesData->GspData[iplane_ind]->coef_mass;

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt     = sim->cp_min_opt;
  int cp_min_cg      = sim->cp_min_cg;

  int istate         = gs.istate_ind;
  int ncoef          = gs.numPoints;
  int len_nhc        = gs.len_nhc_cp; 
  int num_nhc        = gs.num_nhc_cp;
  int nck_nhc        = gs.nck_nhc_cp;
  complex *psi_g     = gs.packedPlaneData; 
  complex *forces    = gs.packedForceData; 
  complex *forcesold = gs.packedVelData; // for minimization not cp
  complex *vpsi_g    = gs.packedVelData; // for cp not minimization
  double ***fNHC     = gs.fNHC;
  double ***vNHC     = gs.vNHC;
  double ***xNHC     = gs.xNHC;
  double ***xNHCP    = gs.xNHCP;
  double *mNHC       = gs.mNHC;
  double *v0NHC      = gs.v0NHC;
  double *a2NHC      = gs.a2NHC;
  double *a4NHC      = gs.a4NHC;
  double kTCP        = gs.kTCP;
  int nkx0_red       = gs.nkx0_red;
  int nkx0_uni       = gs.nkx0_uni;
  int nkx0_zero      = gs.nkx0_zero;
  double degfree     = gs.degfree;
  double degfreeNHC  = gs.degfreeNHC;
  double *degFreeSplt= gs.degFreeSplt;
  int *istrNHC       = gs.istrNHC;
  int *iendNHC       = gs.iendNHC;

  double fictEke = 0.0;
  double ekeNhc  = 0.0;
  double potNHC  = 0.0;

//==========================================================================
// III. Set conjugate gradient parameter : Don't reset too often

  ireset_cg = 0;
  if(iteration==1){ireset_cg=1;numReset_cg=0;}
  if(iteration>1 && cp_min_opt==1 && cp_min_cg==1){
    if( (fmagPsi_total>1.05*fmagPsi_total_old && numReset_cg>=5) || 
        (fmagPsi_total>2.0*fmagPsi_total_old)){
      ireset_cg   = 1;
      numReset_cg = 0;
      if(thisIndex.x==0 && thisIndex.y==0){
         CkPrintf("----------------------------------------------\n");
         CkPrintf("   Resetting Conjugate Gradient!              \n");
         CkPrintf("----------------------------------------------\n");
      }//endif
    }//endif
  }//endif

  double gamma_conj_grad = 0.0; // start CG up again
  if( (cp_min_opt==1) && (cp_min_cg == 1) && (ireset_cg==0)){
    gamma_conj_grad = fovlap/fovlap_old; //continue evolving CG
  }//endif
  ireset_cg = 0;
  numReset_cg++;

//==========================================================================
// IV. Evolve the states using the forces/conjugate direction

//---------------------------------------------------------------
// (A) Debug output before integration

#define _CP_DEBUG_DYNAMICS_OFF
#ifdef _CP_DEBUG_DYNAMICS_
  if(cp_min_opt!=1){
    if(istate<3){
      char forcefile[80];
      snprintf(forcefile,80,"Bpsi_force_%d_%d_%d.out",thisIndex.x,thisIndex.y,iteration);
      FILE *fp=fopen(forcefile,"w");
      for(int i=0;i <ncoef;i++){
	  fprintf(fp,"%d %d %d : %.10g %.10g : %g %g : %g %g : %g\n",
                  k_x[i], k_y[i], k_z[i],psi_g[i].re,psi_g[i].im,
                  vpsi_g[i].re,vpsi_g[i].im,
                  forces[i].re,forces[i].im, coef_mass[i]);
      }//endfor
      fclose(fp);
    }//endif
  }//endif
#endif

//---------------------------------------------------------------
// (B) Numerical integration

#if CMK_TRACE_ENABLED
      double StartTime=CmiWallTimer();
#endif

#define _GLENN_CHECK_DYNAMICS_OFF_
#ifdef _GLENN_CHECK_DYNAMICS_
  if(iteration==1){bzero(vpsi_g,ncoef*sizeof(complex));}
  if(thisIndex.x==0&&thisIndex.y==0){CkPrintf("Before Integrate : iteration %d\n",iteration);}
#endif

#ifdef _GLENN_CHECK_INTEGRATE_
  for(int i=0;i<ncoef;i++){
    if( (k_x[i]==0 &&k_y[i]==1 && k_z[i]==4) ||
        (k_x[i]==2 &&k_y[i]==1 && k_z[i]==3)){
     if(thisIndex.x==0 || thisIndex.x==127){
       CkPrintf(" Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g %.15g %.15g\n",
		 gs.istate_ind+1,k_x[i],k_y[i],k_z[i],
  	         psi_g[i].re,psi_g[i].im,
                 forces[i].re,forces[i].im);
     }//endif
   }//endif
  }//endfor
#endif

#ifdef _CP_DEBUG_SCALC_ONLY_ 
  bzero(forces,ncoef*sizeof(complex));
  if(cp_min_opt==0){
    psi_g = gs.packedPlaneDataTemp2;
  }
#endif

#ifdef _NAN_CHECK_
  for(int i=0; i < ncoef; i++)
  {
      CkAssert(finite(vpsi_g[i].re));
      CkAssert(finite(vpsi_g[i].im));
      CkAssert(finite(psi_g[i].re));
      CkAssert(finite(psi_g[i].im));
      CkAssert(finite(forces[i].re));
      CkAssert(finite(forces[i].im));
  }
#endif

  fictEke = 0.0; ekeNhc=0.0; potNHC=0.0;
  CPINTEGRATE::CP_integrate(ncoef,istate,iteration,forces,forcesold,psi_g,
               coef_mass,k_x,k_y,k_z,len_nhc,num_nhc,nck_nhc,fNHC,vNHC,xNHC,xNHCP,
   	       mNHC,v0NHC,a2NHC,a4NHC,kTCP,gamma_conj_grad,&fictEke,
               nkx0_red,nkx0_uni,nkx0_zero,&ekeNhc,&potNHC,degfree,degfreeNHC,
	       degFreeSplt,istrNHC,iendNHC,halfStepEvolve,config.nfreq_cpintegrate);
  halfStepEvolve = 1;


#ifdef _GLENN_CHECK_INTEGRATE_
  for(int i=0;i<ncoef;i++){
    if( (k_x[i]==0 &&k_y[i]==1 && k_z[i]==4) ||
        (k_x[i]==2 &&k_y[i]==1 && k_z[i]==3)){
     if(thisIndex.x==0 || thisIndex.x==127){
       CkPrintf(" Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g %.15g %.15g\n",
		 gs.istate_ind+1,k_x[i],k_y[i],k_z[i],
  	         psi_g[i].re,psi_g[i].im,
                 forces[i].re,forces[i].im);
     }//endif
   }//endif
  }//endfor
#endif

#if CMK_TRACE_ENABLED
      traceUserBracketEvent(IntegrateModForces_, StartTime, CmiWallTimer());
#endif

//---------------------------------------------------------------
// (C) Debug output after integration

#ifdef _NAN_CHECK_
  for(int i=0; i < ncoef; i++)
  {
      CkAssert(finite(vpsi_g[i].re));
      CkAssert(finite(vpsi_g[i].im));
      CkAssert(finite(psi_g[i].re));
      CkAssert(finite(psi_g[i].im));
      CkAssert(finite(forces[i].re));
      CkAssert(finite(forces[i].im));
  }
#endif

#ifdef _CP_DEBUG_DYNAMICS_
  if(cp_min_opt!=1){
    if(thisIndex.x<3){
      char forcefile[80];
      snprintf(forcefile,80,"Apsi_force_%d_%d_%d.out",thisIndex.x,thisIndex.y,iteration);
      FILE *fp=fopen(forcefile,"w");
      for(int i=0;i <ncoef;i++){
	  fprintf(fp,"%d %d %d : %.10g %.10g : %g %g : %g %g : %g\n",
                  k_x[i], k_y[i], k_z[i],psi_g[i].re,psi_g[i].im,
                  vpsi_g[i].re,vpsi_g[i].im,
                  forces[i].re,forces[i].im, coef_mass[i]);
      }//endfor
      fclose(fp);
    }//endif
  }//endif
  if(iteration==2){CkPrintf("later debugging dyn\n");CkExit();}
#endif


//==========================================================================
// V. Contribute FictKe : Output and store in energy group

  double redSize = ((double) (nchareG*nstates));
  double sendme[5];
  sendme[0] = 2.0*BOLTZ*(fictEke/degfree)/redSize;
  sendme[1] = 2.0*BOLTZ*(ekeNhc/degfreeNHC)/redSize;
  sendme[2] = fictEke;
  sendme[3] = ekeNhc;
  sendme[4] = potNHC;
  contribute(5*sizeof(double),sendme,CkReduction::sum_double, 
               CkCallback(CkIndex_InstanceController::printFictEke(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy));

//==========================================================================
// VI. Pack up and set the flag that indicating you've finished integrating.

  gs.fictEke_ret      = fictEke;
  gs.ekeNhc_ret       = ekeNhc;
  gs.potNHC_ret       = potNHC;
  finishedCpIntegrate = 1;

//------------------------------------------------------------------------------
   }// end CP_State_GSpacePlane::integrateModForce
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::sendRedPsi() {
//==============================================================================


  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  complex *sendData = gs.packedPlaneData;
  int isend         = thisIndex.y;       // my g-space chare index
  int  *num_send    = RCommPkg[isend].num_send;
  int **lst_send    = RCommPkg[isend].lst_send;
  int num_send_tot  = RCommPkg[isend].num_send_tot;

  if(finishedCpIntegrate==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Flow of Control Error : Attempting to sendredpsi\n");
    CkPrintf("without integration\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  
//==============================================================================

  int iii=0; int jjj=0;
  if(num_send_tot>0){ 
    for(int irecv = 0; irecv < nchareG; irecv ++){
    
      int ncoef       = num_send[irecv];
      jjj += ncoef;
      if(ncoef>0){
         GSRedPsiMsg *msg  = new (ncoef,8*sizeof(int)) GSRedPsiMsg;
         if(config.prioFFTMsg){
            CkSetQueueing(msg, CK_QUEUEING_IFIFO);
            *(int*)CkPriorityPtr(msg) = config.rsfftpriority + 
                                     thisIndex.x*gs.planeSize[0]+thisIndex.y;
         }//endif
         if(ncoef>0){iii++;}
         msg->size        = ncoef;
         msg->senderIndex = isend;  // my gspace chare index
         complex *msgData = msg->data;
         for (int i=0;i<ncoef; i++){msgData[i] = sendData[lst_send[irecv][i]];}
         UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x,irecv).acceptRedPsi(msg);
      }//endif

    }//endfor

  }//endif : no one to which I have to send

//==============================================================================
// Check for errors 

  if(iii!=num_send_tot){
    CkPrintf("Error in GSchare %d %d : %d %d : sendredPsi.1\n",thisIndex.x,thisIndex.y,
 	                                           num_send_tot,iii);
    CkExit();
  }//endif
  if(jjj != gs.nkx0_uni-gs.nkx0_zero){
    CkPrintf("Error in GSchare %d %d : %d %d\n sendredPsi.3",thisIndex.x,thisIndex.y,
	                                        jjj,gs.nkx0_uni);
    CkExit();
  }//endif

//==============================================================================
// I sent the stuff

  iSentRedPsi = 1;

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptRedPsi(GSRedPsiMsg *msg) {
//==============================================================================

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  int ncoef         = msg->size;
  int isend         = msg->senderIndex;
  complex *msgData  = msg->data;
  complex *recvData = gs.packedRedPsi;
  int irecv         = thisIndex.y;       // my g-space chare index
  int  *num_recv    = RCommPkg[irecv].num_recv;
  int **lst_recv    = RCommPkg[irecv].lst_recv;

//==============================================================================
// unpack

  if(num_recv[isend]!=ncoef){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Number sent not equal to number reciever expected \n");
    CkPrintf("Sender %d size %d Reciever %d size %d \n",isend,ncoef,irecv,num_recv[isend]);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(countRedPsi==0){jtemp=0;}
  jtemp+= ncoef;

  for(int i=0;i<ncoef;i++){
    recvData[lst_recv[isend][i]].re= msgData[i].re;
    recvData[lst_recv[isend][i]].im=-msgData[i].im;
  }//endfor

  delete msg;

//==============================================================================
// Done

  countRedPsi++;
  if(countRedPsi==numRecvRedPsi){
    countRedPsi=0;
    iRecvRedPsi=1;
    if(jtemp!=gs.nkx0_red){
      CkPrintf("Error in GSchare recv cnt %d %d : %d %d\n",thisIndex.x,thisIndex.y,
 	                                                   gs.nkx0_red,jtemp);
      CkExit();
    }//endif
    // If sent before I received then I resume
    if(iSentRedPsi==1){ 
      UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl(); 
    }//endif
  }//endif

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doneRedPsiIntegrate() {

    int ncoef    = gs.nkx0_red;
    if(ncoef>0){
      complex *psi     = gs.packedPlaneData;
      complex *psi_red = gs.packedRedPsi;
      for(int i=0;i<ncoef;i++){psi[i]=psi_red[i];}
    }//endif

}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::sendPsi() {
//==============================================================================
// Error checking

#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("sendpsi %d %d\n",thisIndex.y,cleanExitCalled);
#endif

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt       = sim->cp_min_opt;

  if(finishedCpIntegrate==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't sendPsi without completing integrate\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(iteration>0){
    if(iRecvRedPsi!=1 || iSentRedPsi!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you can't sendPsi without receiving Redpsi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//==============================================================================
// Prepare the data : If cp dynamics is going, save the non-orthogonal puppies.

  acceptedPsi    = false;

  complex *psi   = gs.packedPlaneData;
  int numPoints  = gs.numPoints;

  if(cp_min_opt==0){
     int ncoef     = gs.numPoints;
     complex *scr  = gs.packedPlaneDataScr; //save non-orthog psi
     CmiMemcpy(scr,psi,sizeof(complex)*ncoef);
  }//endif

#ifndef PAIRCALC_TEST_DUMP
#ifndef _CP_DEBUG_ORTHO_OFF_
  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i] *= rad2i;}
  }//endif
#endif
#endif

//==============================================================================
// Debugging 

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  for(int i=0;i<config.numChunksSym;i++){CkAssert(countPsiO[i]==0);}
#endif
  
#ifdef _CP_GS_DUMP_PSI_
    dumpMatrix("psiBfp",(double *)psi, 1, gs.numPoints*2,
                     thisIndex.y,thisIndex.x,thisIndex.x,0,false);     
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  double testvalue=0.00000001;
  if(savedpsiBfp==NULL){ // load it
      savedpsiBfp= new complex[gs.numPoints];
      loadMatrix("psiBfp",(double *)savedpsiBfp, 1, gs.numPoints*2,
                        thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  for(int i=0;i<gs.numPoints;i++){
      if(fabs(psi[i].re-savedpsiBfp[i].re)>testvalue){
	  fprintf(stderr, "GSP [%d,%d] %d element psi  %.10g not %.10g\n",
                  thisIndex.x, thisIndex.y,i, psi[i].re, savedpsiBfp[i].re);
      }//endif
      CkAssert(fabs(psi[i].re-savedpsiBfp[i].re)<testvalue);
      CkAssert(fabs(psi[i].im-savedpsiBfp[i].im)<testvalue);
  }//endfor
#endif

//==============================================================================
// Start the calculator

#ifndef _CP_DEBUG_ORTHO_OFF_
  symmPCmgr.sendLeftData(numPoints, psi, false);
  /// Symm loop PC chares in the top left [*,0,0,*] will not receive any right matrix data. Hence, if you're in such a PC's block, dont send right
  if(thisIndex.x >= symmPCmgr.pcCfg.grainSize)
      symmPCmgr.sendRightData(numPoints, psi, false);
#else
  acceptedPsi=true;
  if((iteration==config.maxIter || exitFlag==1) && cp_min_opt==1 && config.stateOutput==0)
  	UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).readyToExit();
#endif

//----------------------------------------------------------------------------
    }// end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNewPsi(CkReductionMsg *msg){
//=============================================================================
// (0) Fuss with the redundant psis

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt       = sim->cp_min_opt;
  if(iteration>0){
    if(iRecvRedPsi!=1 || iSentRedPsi!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you can't acceptPsi without receiving Redpsi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//=============================================================================
// (0) Nan Check and output

  int N           = msg->getSize()/sizeof(complex);
  complex *data   = (complex *)msg->getData();
  int offset      = msg->getUserFlag();  if(offset<0){offset=0;}
  complex *psi    = gs.packedPlaneData;
  int chunksize   = gs.numPoints/config.numChunksSym;
  int chunkoffset = offset*chunksize; // how far into the points this contribution lies

#ifdef _NAN_CHECK_
  for(int i=0;i<N ;i++){
      CkAssert(finite(((complex *) msg->getData())[i].re));
      CkAssert(finite(((complex *) msg->getData())[i].im));
  }//endfor
#endif

/*
  if(thisIndex.y==0){
      CkPrintf("PSI [%d %d], offset %d chunkoffset %d N %d countPsiO %d\n", 
           thisIndex.x, thisIndex.y, offset, chunkoffset, N, countPsiO[offset]);
  }//endif
*/

//=============================================================================
// (I) Unpack the contribution to newpsi (orthonormal psi)

  int idest=chunkoffset;

  if(countPsiO[offset]<1){
    CmiMemcpy(&(psi[idest]), &(data[0]), N*sizeof(complex)); //slower?
    //for(int i=0; i<N; i++,idest++){psi[idest] = data[i];}
  }else{
    //    for(int i=0; i<N; i++,idest++){psi[idest] += data[i];}
    fastAdd((double *) &psi[idest],(double *)data,N*2);
  }//endif

  delete msg;

//=============================================================================
// (II) If you have got it all : Rescale it, produce some output

  countPsi++;//psi arrives in as many as 2 *numblock reductions
  countPsiO[offset]++;//psi arrives in as many as 2 
  if(countPsi==AllPsiExpected){ 
    thisProxy(thisIndex.x,thisIndex.y).doNewPsi();
#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("acceptpsi %d %d\n",thisIndex.y,cleanExitCalled);
#endif
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/**
 * When the psi calculator is sending to us directly
 */
//==============================================================================
void CP_State_GSpacePlane::acceptNewPsi(partialResultMsg *msg){
//=============================================================================
// (0) Fuss with the redundant psis

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt       = sim->cp_min_opt;
  if(iteration>0){
    if(iRecvRedPsi!=1 || iSentRedPsi!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you can't acceptPsi without receiving Redpsi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//=============================================================================
// (0) Check for Nans



  int N           = msg->N;
  complex *data   = msg->result;
  int offset      = msg->myoffset;  if(offset<0){offset=0;}
  complex *psi    = gs.packedPlaneData;
  int chunksize   = gs.numPoints/config.numChunksSym;
  int chunkoffset = offset*chunksize; // how far into the points this contribution lies

#ifdef _NAN_CHECK_
  for(int i=0;i<N ;i++){
      if((!finite(data[i].re)) || (!finite(data[i].im))){
	 CkPrintf("GSP [%d,%d] acceptNewPsi offset %d %d of %d is nan\n",
               thisIndex.x,thisIndex.y, msg->myoffset, i,N);
      }
      CkAssert(finite(data[i].re));
      CkAssert(finite(data[i].im));
  }//endif
#endif

/*
  if(thisIndex.y==0){
      CkPrintf("PSI [%d %d], offset %d chunkoffset %d N %d countPsiO %d\n", 
               thisIndex.x, thisIndex.y, offset, chunkoffset, N, countPsiO[offset]);
  }//endif
*/

//=============================================================================
// (I) Unpack the contribution to newpsi (orthonormal psi)

  int idest = chunkoffset;
  if(countPsiO[offset]<1){
    CmiMemcpy(&(psi[idest]), &(data[0]), N*sizeof(complex)); //slower?
    //    for(int i=0; i<N; i++,idest++){psi[idest] = data[i];}
  }else{
    //for(int i=0; i<N; i++,idest++){psi[idest] += data[i];}
    fastAdd((double *) &psi[idest],(double *)data,N*2);
  }//endif

  delete msg;

//==============================================================================
// When you are done, continue

  countPsi++;         //psi arrives in as many as 2 *numblock * numgrain reductions
  countPsiO[offset]++;//psi arrives in as many as 2 * numgrain
  //
  if(countPsi==AllPsiExpected){ 
    thisProxy(thisIndex.x,thisIndex.y).doNewPsi();
#ifdef _CP_DEBUG_STATEG_VERBOSE_
    if(thisIndex.x==0){CkPrintf("aceeptpsi %d %d\n",thisIndex.y,cleanExitCalled);}
#endif
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
/**
 * All Psi have arrived, finish the new Psi process
 */
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doNewPsi(){
//=============================================================================
// (0) Fuss with the redundant psis

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt       = sim->cp_min_opt;
  int cp_min_update    = sim->cp_min_update;

  if(iteration>0){
    if(iRecvRedPsi!=1 || iSentRedPsi!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you can't acceptPsi without receiving Redpsi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//=============================================================================
// (I) If you have got it all : Rescale it, produce some output

#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("donewpsi %d %d\n",thisIndex.y,cleanExitCalled);
#endif

  CkAssert(countPsi==AllPsiExpected); 

  complex *psi  = gs.packedPlaneData;
#ifdef _NAN_CHECK_
  for(int i=0;i<gs.numPoints ;i++){
      CkAssert(finite(psi[i].re));
      CkAssert(finite(psi[i].im));
  }//endfor
#endif

//=============================================================================
// (A) Reset counters and rescale the kx=0 stuff

  acceptedPsi = true;
  countPsi    = 0;
  bzero(countPsiO,config.numChunksSym*sizeof(int));

#ifndef PAIRCALC_TEST_DUMP
  if(gs.ihave_kx0==1){
    double rad2 = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i] *= rad2;}
  }//endif
#endif

//=============================================================================
// (B) Generate some screen output of orthogonal psi

  if(iteration>0){screenOutputPsi(iteration);}

//=============================================================================
// (D) Go back to the top or exit

  if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==0)
	UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).readyToExit();

  if((iteration==config.maxIter || exitFlag==1) && cp_min_opt==1 && config.stateOutput==0)
	UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).readyToExit();

//=============================================================================
// (E) Debug psi

#ifdef _CP_GS_DUMP_PSI_
  dumpMatrix("psiAf",(double *)psi, 1, gs.numPoints*2,thisIndex.y,thisIndex.x,
                    thisIndex.x,0,false);
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  double testvalue=0.00000001;
  if(savedpsiAf==NULL){
      savedpsiAf= new complex[gs.numPoints];
      loadMatrix("psiAf",(double *)savedpsiAf, 1, gs.numPoints*2,
                        thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  for(int i=0;i<gs.numPoints;i++){
      if(fabs(psi[i].re-savedpsiAf[i].re)>testvalue){
        fprintf(stderr, "GSP [%d,%d] %d element psi  %.10g not %.10g\n",
        thisIndex.x, thisIndex.y,i, psi[i].re, savedpsiAf[i].re);
      }//endif
      CkAssert(fabs(psi[i].re-savedpsiAf[i].re)<testvalue);
      CkAssert(fabs(psi[i].im-savedpsiAf[i].im)<testvalue);
  }//endfor
#endif

//=============================================================================
// (E) Reset psi 

  if(cp_min_opt==1 && cp_min_update==0 && iteration>0){
    CmiMemcpy(gs.packedPlaneData,gs.packedPlaneDataTemp,
	      sizeof(complex)*gs.numPoints);
    memset(gs.packedVelData,0,sizeof(complex)*gs.numPoints);
  }//endif

  if(cp_min_opt==1 && cp_min_update==0 && iteration==0){
    CmiMemcpy(gs.packedPlaneDataTemp,gs.packedPlaneData,
	      sizeof(complex)*gs.numPoints);
    memset(gs.packedVelData,0,sizeof(complex)*gs.numPoints);
  }//endif

//==============================================================================
// Back to the threaded loop.
#ifdef BARRIER_CP_GSPACE_PSI
  int wehaveours=1;
  contribute(sizeof(int),&wehaveours,CkReduction::sum_int,CkCallback(CkIndex_GSpaceDriver::allDonePsi(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
#else
  UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl();
#endif

//----------------------------------------------------------------------------
  }//end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
//! dynamics triggers send of orthoT to asymm calc when psi is done
void CP_State_GSpacePlane::launchOrthoT(){
//=============================================================================
  CkPrintf("[%d,%d] launchOrthoT \n",thisIndex.x, thisIndex.y);
  if(thisIndex.x==0 && thisIndex.y==0)
    myOrtho.sendOrthoTtoAsymm();
//----------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::sendRedPsiV(){

#ifdef DEBUG_CP_GSPACE_PSIV
	CkPrintf("GSpace[%d,%d] sendRedPsiV: Going to send redundant PsiV data\n",thisIndex.x,thisIndex.y);
#endif
//==============================================================================
// I) Local Pointers

  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  double *coef_mass = eesData->GspData[iplane_ind]->coef_mass;

  int ncoef         = gSpaceNumPoints;
  complex *psi      = gs.packedPlaneData;
  complex *vpsi     = gs.packedVelData;       // to evolve psi to time, t.
  complex *forces   = gs.packedForceData;
  double ***xNHC    = gs.xNHC;
  double ***xNHCP   = gs.xNHCP;
  double ***vNHC    = gs.vNHC;
  double ***fNHC    = gs.fNHC;          
  double *mNHC      = gs.mNHC;          
  int len_nhc_cp    = gs.len_nhc_cp;
  int num_nhc_cp    = gs.num_nhc_cp;
  int nck_nhc_cp    = gs.nck_nhc_cp;
  int nkx0_red      = gs.nkx0_red;
  int nkx0_uni      = gs.nkx0_uni;
  int nkx0_zero     = gs.nkx0_zero;
  double kTCP       = gs.kTCP;

//=============================================================================
// II) Sync yourself with psi by integrating to time t if output has not done it
//     you have the wrong force but thats OK until you put in a better rotation

  halfStepEvolve = 1; 

#ifdef JUNK
  halfStepEvolve = 0; // do the 1/2 step update now
  CPINTEGRATE::cp_evolve_vel(ncoef,forces,vpsi,coef_mass,
                             len_nhc_cp,num_nhc_cp,nck_nhc_cp,fNHC,vNHC,xNHC,xNHCP,mNHC,
                             gs.v0NHC,gs.a2NHC,gs.a4NHC,kTCP,nkx0_red,nkx0_uni,nkx0_zero,
                             2,iteration,gs.degfree,gs.degfreeNHC,gs.degFreeSplt,
                             gs.istrNHC,gs.iendNHC,1);
#endif

//=============================================================================
// III) We still have these funky g=0 plane guys that may be on other procs and
//      we have sync those guys, also, or PC won't work correctly

  CPcharmParaInfo *sim       = CPcharmParaInfo::get();
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  complex *sendData = gs.packedVelData;
  int isend         = thisIndex.y;       // my g-space chare index
  int  *num_send    = RCommPkg[isend].num_send;
  int **lst_send    = RCommPkg[isend].lst_send;
  int num_send_tot  = RCommPkg[isend].num_send_tot;

  int iii=0; int jjj=0;
  if(num_send_tot>0){ 
    for(int irecv = 0; irecv < nchareG; irecv ++){
    
      int ncoef       = num_send[irecv];
      jjj += ncoef;
      if(ncoef>0){
         GSRedPsiMsg *msg  = new (ncoef,8*sizeof(int)) GSRedPsiMsg;
         if(config.prioFFTMsg){
            CkSetQueueing(msg, CK_QUEUEING_IFIFO);
            *(int*)CkPriorityPtr(msg) = config.rsfftpriority + 
                                     thisIndex.x*gs.planeSize[0]+thisIndex.y;
         }//endif
         if(ncoef>0){iii++;}
         msg->size        = ncoef;
         msg->senderIndex = isend;  // my gspace chare index
         complex *msgData = msg->data;
         for (int i=0;i<ncoef; i++){msgData[i] = sendData[lst_send[irecv][i]];}
         UgSpacePlaneProxy[thisInstance.proxyOffset](thisIndex.x,irecv).acceptRedPsiV(msg);
      }//endif

    }//endfor

  }//endif : no one to which I have to send

//==============================================================================
// Check for errors 

  if(iii!=num_send_tot){
    CkPrintf("Error in GSchare %d %d : %d %d : sendRedPsiV.1\n",thisIndex.x,thisIndex.y,
 	                                           num_send_tot,iii);
    CkExit();
  }//endif
  if(numRecvRedPsi==0 && gs.nkx0_red>0){
    CkPrintf("Error in GSchare %d %d : %d %d : sendRedPsiV.2\n",thisIndex.x,thisIndex.y,
	                                        numRecvRedPsi,gs.nkx0_red);
    CkExit();
  }//endif
  if(jjj != gs.nkx0_uni-gs.nkx0_zero){
    CkPrintf("Error in GSchare %d %d : %d %d : sendRedPsiV.3\n",thisIndex.x,thisIndex.y,
	                                        jjj,gs.nkx0_uni);
    CkExit();
  }//endif

//==============================================================================
// I send the stuff and I need a new velocity

 iSentRedPsiV = 1;
 acceptedVPsi = false;

//==============================================================================
  }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptRedPsiV(GSRedPsiMsg *msg) {
//==============================================================================

  CPcharmParaInfo *sim       = CPcharmParaInfo::get();
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  int ncoef         = msg->size;
  int isend         = msg->senderIndex;
  complex *msgData  = msg->data;
  complex *recvData = gs.packedRedPsiV;

  int irecv         = thisIndex.y;       // my g-space chare index
  int  *num_recv    = RCommPkg[irecv].num_recv;
  int **lst_recv    = RCommPkg[irecv].lst_recv;

#ifdef DEBUG_CP_GSPACE_PSIV
		CkPrintf("GSpace[%d,%d] acceptRedPsiV Received redundant PsiV values from sender %d\n",thisIndex.x,thisIndex.y,isend);
#endif
//==============================================================================
// unpack

  if(num_recv[isend]!=ncoef){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Number sent not equal to number reciever expected \n");
    CkPrintf("Sender %d size %d Reciever %d size %d \n",isend,ncoef,irecv,num_recv[isend]);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(countRedPsiV==0){jtemp=0;}
  jtemp+= ncoef;

  for(int i=0;i<ncoef;i++){
    recvData[lst_recv[isend][i]].re= msgData[i].re;
    recvData[lst_recv[isend][i]].im=-msgData[i].im;
  }//endfor

  delete msg;

//==============================================================================
// Done

  countRedPsiV++;
  if(countRedPsiV==numRecvRedPsi){
#ifdef DEBUG_CP_GSPACE_PSIV
		CkPrintf("GSpace[%d,%d] acceptRedPsiV received all %d GSRedPsi messages carrying redundant PsiV data\n",thisIndex.x,thisIndex.y,countRedPsiV);
#endif
    countRedPsiV = 0;
    iRecvRedPsiV  = 1;
    if(jtemp!=gs.nkx0_red){
      CkPrintf("Error in GSchare recv cnt %d %d : %d %d\n",thisIndex.x,thisIndex.y,
 	                                                   gs.nkx0_red,jtemp);
      CkExit();
    }//endif
    // If sent before I received then I resume
    if(iSentRedPsiV == 1){
        UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl();
    }//endif
  }//endif

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doneRedPsiVIntegrate() {

  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  double *coef_mass = eesData->GspData[iplane_ind]->coef_mass;
  int ncoef         = gs.numPoints;
  int ncoef_red     = gs.nkx0_red;
  complex *vpsi     = gs.packedVelData;
  int istrt0        = gs.nkx0_red;
  int istrt         = gs.nkx0_red+gs.nkx0_zero;
  int iend          = gs.nkx0_red+gs.nkx0_uni;

  if(ncoef_red>0){
    complex *vpsi_red = gs.packedRedPsiV;
    for(int i=0;i<ncoef_red;i++){vpsi[i]=vpsi_red[i];}
  }//endif

  ake_old = 0.0;
  for(int i=istrt0;i<istrt;i++){ake_old += vpsi[i].getMagSqr()*coef_mass[i];}       // g=0
  for(int i=istrt;i<iend;i++)  {ake_old += vpsi[i].getMagSqr()*(2.0*coef_mass[i]);} // gx=0
  for(int i=iend;i<ncoef;i++)  {ake_old += vpsi[i].getMagSqr()*coef_mass[i];}       // gx!=0

}//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendPsiV() {
//==============================================================================
// Error Check

  CPcharmParaInfo *sim       = CPcharmParaInfo::get();
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("GSpace[%d,%d] Error: You can't sendPsiV() without sending/receiving the redundant psiV values around: finished %d %d : %d %d\n",thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkAbort("Error: GSpace cannot sendPsiV() without sending/receiving the redundant psi values around\n");
  }//endif

#ifdef DEBUG_CP_GSPACE_PSIV
	CkPrintf("GSpace[%d,%d] sendPsiV\n",thisIndex.x,thisIndex.y);
#endif

  acceptedVPsi = false;

//==============================================================================
//
  nrotation++;
  iterRotation = iteration+1;

  int ncoef     = gs.numPoints;
  complex *data = gs.packedVelData;
  complex *scr  = gs.packedPlaneDataScr;  // replace no-ortho psi 
  complex *psi  = gs.packedPlaneData;     // by orthonormal psi when norb rotating
  CmiMemcpy(scr,psi,sizeof(complex)*ncoef);

#ifndef PAIRCALC_TEST_DUMP
  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){data[i] *= rad2i;}
  }//endif
#endif

  int numPoints = gs.numPoints;
  symmPCmgr.sendLeftData(numPoints,data,true);
  /// Symm loop PC chares in the top left [*,0,0,*] will not receive any right matrix data. Hence, if you're in such a PC's block, dont send right
  if(thisIndex.x >= symmPCmgr.pcCfg.grainSize)
      symmPCmgr.sendRightData(numPoints,data,true);

//----------------------------------------------------------------------------
}// end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNewPsiV(CkReductionMsg *msg){
//=============================================================================
// (I) Local pointers

  int N           = msg->getSize()/sizeof(complex);
  complex *data   = (complex *)msg->getData();
  int offset      = msg->getUserFlag();  if(offset<0){offset=0;}

  complex *vpsi   = gs.packedVelData;
  int chunksize   = gs.numPoints/config.numChunksSym;
  int chunkoffset = offset*chunksize;; // how far into the points this contribution lies

#ifdef DEBUG_CP_GSPACE_PSIV
		CkPrintf("GSpace[%d,%d] acceptNewPsiV(reductionMsg) PCs have sent new PsiV data\n",thisIndex.x,thisIndex.y);
#endif
  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("GSpace[%d,%d] Error: You can't acceptNewPsiV() without sending/receiving the redundant PsiV values around: finished %d %d : %d %d\n",thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkAbort("Error: GSpace cannot acceptNewPsiV() without sending/receiving the redundant PsiV values around\n");
  }//endif

//=============================================================================
// (II) Unpack the contribution to newpsi (orthonormal psi)

  int idest=chunkoffset;
  if(countVPsiO[offset]<1){
     CmiMemcpy(&(vpsi[idest]), &(data[0]), N*sizeof(complex));//slower?
     //for(int i=0; i<N; i++,idest++){vpsi[idest] = data[i];}
  }else{
    //    for(int i=0; i<N; i++,idest++){vpsi[idest] += data[i];}
    fastAdd((double *) &vpsi[idest],(double *)data,N*2);    
  }//endif  

  delete msg;

//=============================================================================
// (III) When all has arrive, onward to victory

  countVPsi++;         //psi arrives in as many as 2 reductions
  countVPsiO[offset]++;//psi arrives in as many as 2 reductions

  if(countVPsi==AllPsiExpected){ 
#ifdef DEBUG_CP_GSPACE_PSIV
	CkPrintf("GSpace[%d,%d] Received all PsiV data from PCs (%d reductions).\n",thisIndex.x,thisIndex.y,AllPsiExpected);
#endif
    thisProxy(thisIndex.x,thisIndex.y).doNewPsiV();
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNewPsiV(partialResultMsg *msg){
//=============================================================================
// (I) Local pointers

  int N           = msg->N;
  complex *data   = msg->result;
  int offset      = msg->myoffset;  if(offset<0){offset=0;}

#ifdef _NAN_CHECK_
  for(int i=0; i < N; i++)
  {
      CkAssert(finite(data[i].re));
      CkAssert(finite(data[i].im));
  }
#endif

  complex *vpsi   = gs.packedVelData;
  int chunksize   = gs.numPoints/config.numChunksSym;
  int chunkoffset = offset*chunksize;; // how far into the points this contribution lies

#ifdef DEBUG_CP_GSPACE_PSIV
		CkPrintf("GSpace[%d,%d] acceptNewPsiV(partialResultMsg) Received new PsiV data (msg %d of %d) from PC [%d,%d,%d,%d] (offset %d)\n",thisIndex.x,thisIndex.y,countVPsi+1,AllPsiExpected,
		msg->sndr.w,msg->sndr.x,msg->sndr.y,msg->sndr.z,offset);
#endif
  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("GSpace[%d,%d] Error: You can't acceptNewPsiV() without sending/receiving the redundant PsiV values around: finished %d %d : %d %d\n",thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkAbort("Error: GSpace cannot acceptNewPsiV() without sending/receiving the redundant PsiV values around\n");
  }//endif

//=============================================================================
// (I) Unpack the contribution to newpsi (orthonormal psi)

  int idest=chunkoffset;
  if(countVPsiO[offset]<1){
    CmiMemcpy(&(vpsi[idest]), &(data[0]), N*sizeof(complex));//slower?
     //for(int i=0; i<N; i++,idest++){vpsi[idest] = data[i];}
  }else{
    //    for(int i=0; i<N; i++,idest++){vpsi[idest] += data[i];}
    fastAdd((double *) &vpsi[idest],(double *)data,N*2);
  }//endif
  
  delete msg;

//=============================================================================
// (II) Continue

  countVPsi++;//psi arrives in as many as 2 reductions
  countVPsiO[offset]++;//psi arrives in as many as 2 reductions

  if(countVPsi==AllPsiExpected){ 
#ifdef DEBUG_CP_GSPACE_PSIV
		CkPrintf("GSpace[%d,%d] Received all PsiV data from PCs (%d messages).\n",thisIndex.x,thisIndex.y,AllPsiExpected);
#endif
    thisProxy(thisIndex.x,thisIndex.y).doNewPsiV();
  }//endif

//----------------------------------------------------------------------------
  }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doNewPsiV(){
//=============================================================================
// (0) Error check
//  CkPrintf("[%d %d] GSP doNewPsiV \n",thisIndex.x, thisIndex.y);
#ifdef DEBUG_CP_GSPACE_PSIV
		CkPrintf("GSpace[%d,%d] doNewPsiV\n",thisIndex.x,thisIndex.y);
#endif

  CkAssert(countVPsi==AllPsiExpected); 

  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("GSpace[%d,%d] Error: You can't doNewPsiV() without sending/receiving the redundant PsiV values around: finished %d %d : %d %d\n",thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkAbort("Error: GSpace cannot doNewPsiV() without sending the redundant/receiving PsiV values around\n");
  }//endif

//=============================================================================
// (I) Reset counters and rescale the kx=0 stuff

  acceptedVPsi    = true;
  countVPsi       = 0;

  complex *vpsi = gs.packedVelData;

#ifndef PAIRCALC_TEST_DUMP
  for(int i=0;i<config.numChunksSym;i++){countVPsiO[i]=0;}
  if(gs.ihave_kx0==1){
    double rad2 = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){vpsi[i] *= rad2;}
  }//endif
#endif

//=============================================================================
// III) Replace by finite difference until update is better

 CPcharmParaInfo *sim = CPcharmParaInfo::get();
 eesCache *eesData    = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();

 double *coef_mass    = eesData->GspData[iplane_ind]->coef_mass;
 double dt           = sim->dt;
 double dt2          = 2.0*dt;
 int ncoef           = gs.numPoints;
 complex *psi_g      = gs.packedPlaneData;
 complex *psi_g_tmp  = gs.packedPlaneDataTemp;
 int istrt0          = gs.nkx0_red;
 int istrt           = gs.nkx0_red+gs.nkx0_zero;
 int iend            = gs.nkx0_red+gs.nkx0_uni;

 for(int i=0;i<ncoef;i++){
   double vre = (psi_g[i].re-psi_g_tmp[i].re)/dt;
   double vim = (psi_g[i].im-psi_g_tmp[i].im)/dt;
   vpsi[i].re = vre;
   vpsi[i].im = vim;
 }//endif

 double ake_new = 0.0;
 for(int i=istrt0;i<istrt;i++){ake_new += vpsi[i].getMagSqr()*coef_mass[i];}       // g=0
 for(int i=istrt;i<iend;i++)  {ake_new += vpsi[i].getMagSqr()*(2.0*coef_mass[i]);} // gx=0
 for(int i=iend;i<ncoef;i++)  {ake_new += vpsi[i].getMagSqr()*coef_mass[i];}       // gx!=0

 if(sim->cp_norb_rot_kescal==1){
   double scale = sqrt(ake_old/ake_new);
   for(int i=0;i<ncoef;i++){
     vpsi[i].re *= scale;
     vpsi[i].im *= scale;
   }//endfor
 }//endif

//=============================================================================
/// II) A Barrier for debugging

#ifdef BARRIER_CP_GSPACE_PSIV
	int wehaveours=1;
	contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
	CkCallback(CkIndex_GSpaceDriver::allDonePsiV(NULL),UgSpaceDriverProxy[thisInstance.proxyOffset]));
#else
	UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).resumeControl();
#endif

//------------------------------------------------------------------------------
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::screenOutputPsi(int iprintout){
//==============================================================================
#ifdef _CP_DEBUG_STATEG_VERBOSE_
  if(thisIndex.x==0){CkPrintf("output %d %d\n",thisIndex.y,cleanExitCalled);}
#endif

  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind]->ka;
  int *k_y          = eesData->GspData[iplane_ind]->kb;
  int *k_z          = eesData->GspData[iplane_ind]->kc;

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt    = sim->cp_min_opt;
  int nstates       = sim->nstates;

  complex *vpsi     = gs.packedVelData;
  complex *psi      = gs.packedPlaneData;       //orthogonal psi
  if(cp_min_opt==0){psi=gs.packedPlaneDataScr;} //non-orthogonal psi

  int ntime = config.maxIter;
  if(cp_min_opt==0){ntime-=1;}

//==============================================================================
// Screen Output

#ifdef _CP_DEBUG_COEF_SCREEN_
  if(iteration<=ntime){
    if(gs.istate_ind==0 || gs.istate_ind==nstates-1){
      for(int i = 0; i < gs.numPoints; i++){
	if(k_x[i]==0 && k_y[i]==1 && k_z[i]==4 ){
	  fprintf(temperScreenFile,"-------------------------------------------------------------------------------\n");
	  fprintf(temperScreenFile,"Iter [%d] Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		  iprintout,
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],psi[i].re,psi[i].im);
          if(cp_min_opt==0){
            double vre=vpsi[i].re;
            double vim=vpsi[i].im;
 	    fprintf(temperScreenFile,"Iter [%d] VPsi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		     iprintout,
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],vre,vim);
	  }//endif
	  fprintf(temperScreenFile,"-------------------------------------------------------------------------------\n");
	}//endif
	if(k_x[i]==2 && k_y[i]==1 && k_z[i]==3){
          double vre=vpsi[i].re;
          double vim=vpsi[i].im;
	  fprintf(temperScreenFile,"-------------------------------------------------------------------------------\n");
	  fprintf(temperScreenFile,"Iter [%d] Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   iprintout,
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],psi[i].re,psi[i].im);
          if(cp_min_opt==0){
 	    fprintf(temperScreenFile,"Iter [%d] VPsi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		     iprintout,
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],vre,vim);
	  }//endif
	  fprintf(temperScreenFile,"-------------------------------------------------------------------------------\n");
	}//endif
      }//endfor
    }//endif
  }//endif
#endif

//==============================================================================
// II) Tell the world you are done with the output

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("GSpace[%d,%d] screenwrite: %d\n",thisIndex.x, thisIndex.y,iteration);
#endif

//----------------------------------------------------------------------------
  }//end routine




//==============================================================================
// Once All chares have completed energy computation and reduction, this 
// storage routine is invoked. Since eke is the last energy, its reduction client
// invokes this guy on chare (0,0). The routine then bcasts its results to the 
// energy group (one per processor) thereby making information available on all 
// procs for tolerence testing via a cklocal type deal.
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::computeEnergies(int param, double d){
//==============================================================================

  switch(param){
      
    case ENERGY_EHART : 
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received Hart\n");
#endif
      ehart_total = d;
      total_energy += d;
      ecount++;
      break;
      
    case ENERGY_ENL :
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received Enl\n");
#endif
      enl_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EKE : 
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received Eke\n");
#endif
      eke_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EGGA :
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received GGA\n");
#endif
      egga_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EEXC :
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received EEXC\n");
#endif
      eexc_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EEXT :
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received EEXT\n");
#endif
      eext_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EWD :
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received EWD\n");
#endif
      ewd_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_FICTEKE : 
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received FICTEKE\n");
#endif
      fictEke_total = d;
      ecount++;
      break;
      
    case ENERGY_FMAG : 
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
      CkPrintf("Received FMAG\n");
#endif
      fmagPsi_total0 = d;
      ecount++;
      break;
      
    default :
      CkAbort("unknown energy");
      break;
  }//end switch

//==============================================================================
// if your debugging or natm_nl==0 you get fewer energies

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int isub =0;
  int natm_nl = sim->natm_nl;
  if(natm_nl==0){isub++;}
#ifdef _CP_DEBUG_SFNL_OFF_
  if(natm_nl!=0){isub++;}
#endif
#ifdef  _CP_DEBUG_RHO_OFF_
  isub+=5;
#endif
#ifdef _CP_DEBUG_HARTEEXT_OFF_
#ifndef  _CP_DEBUG_RHO_OFF_
  isub+=3;
#endif
#endif

  if(thisInstance.idxU.y>0)
    { // you get no rho
      isub+=5;
    }

  int myid = CkMyPe();
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("ecount %d %d %d\n",ecount,NUM_ENERGIES-isub,myid);
#endif
  if(ecount == NUM_ENERGIES-isub){
    EnergyStruct estruct;
    estruct.enl             = enl_total;
    estruct.eke             = eke_total;
    estruct.eext            = eext_total;
    estruct.ehart           = ehart_total;
    estruct.eewald_recip    = ewd_total;
    estruct.egga            = egga_total;
    estruct.eexc            = eexc_total;
    estruct.fictEke         = fictEke_total;
    estruct.totalElecEnergy = total_energy;
    estruct.fmagPsi         = fmagPsi_total0;
    estruct.iteration_gsp   = iteration;
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
    CkPrintf("Bcasting to energygrp %d\n",myid);
#endif
    UegroupProxy[thisInstance.proxyOffset].updateEnergiesFromGS(estruct,thisInstance); // broadcast the electronic energies
                                               //  so that all procs have them
    total_energy        = 0.0;
    ecount              = 0;

  }// got all the energies

//-----------------------------------------------------------------------------
   }// end routine : computenergies
//==============================================================================




//==============================================================================
// Test program : Does ReadFile() parse psi(g) and g-vectors properly
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void testeke(int ncoef,complex *psi_g,int *k_x,int *k_y,int *k_z, int iflag,int index)
//==============================================================================
   {//begin routine
//==============================================================================
// Local pointers

  GENERAL_DATA *general_data = GENERAL_DATA::get();
  CP           *cp           = CP::get();

#include "../class_defs/allclass_strip_gen.h"
#include "../class_defs/allclass_strip_cp.h"

  double gx, gy, gz, g2;
  double *hmati    = gencell->hmati;
  double ecut      = cpcoeffs_info->ecut_psi; // KS-state cutoff in Ryd
  double tpi       = 2.0*M_PI;
  double wght      = 2.0;
  double norm      = 0.0;
  double norm2     = 0.0;
  double eke       = 0.0;
  double eke2      = 0.0;

//==============================================================================
// Compute some eke

  double eke_i[1000];
  int    kxmax, kxmin;
  for(int i=0; i<1000; i++){eke_i[i]=0.0;}

  kxmax=-10000;
  kxmin=100000;
  for(int i = 0; i < ncoef; i++){
    kxmax=(kxmax < k_x[i] ? k_x[i] : kxmax);
    kxmin=(kxmin > k_x[i] ? k_x[i] : kxmin);
  }//endfor

  for(int i = 0; i < ncoef; i++){

    gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
    gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
    gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
    g2 = gx*gx + gy*gy + gz*gz;

    if(g2<=ecut){
       double wght_now = 1.0;
       if(config.doublePack){
         wght_now=2.0;
         if(k_x[i]==0 && k_y[i]<0){wght_now=0.0;}
         if(k_x[i]==0 && k_y[i]==0 && k_z[i]<0){wght_now=0.0;}
         if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){wght_now=1.0;}
       }//endif
       eke_i[k_x[i]-kxmin]+=(wght_now*g2)*psi_g[i].getMagSqr();
       eke       += (wght_now*g2)*psi_g[i].getMagSqr();
       norm      += (wght_now)*psi_g[i].getMagSqr();
       wght_now = (k_x[i]==0 ? 1.0 : wght);
       eke2      += (wght_now*g2)*psi_g[i].getMagSqr();
       norm2     += (wght_now)*psi_g[i].getMagSqr();
    }else{
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkPrintf("Why the cutoff\n");
       CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
       CkExit();
    }//endif

  }/* endfor */

   eke/=2.0;
   eke2/=2.0;

   //   if(index==0){
   //  CkPrintf("hmati :");
   //  for(int i=1;i<=9;i++){CkPrintf(" %g",hmati[i]);}
   //  CkPrintf("\n");
   //}/*endif*/

   CkPrintf("%.12g %.12g %.12g %.12g: %d : eke\n",eke,eke2,norm,norm2,index);

    /****************************************
     FILE* fp;
     char junk[1000];
     sprintf(junk,"eke.%d.out",index);
     fp=fopen(junk,"w");
     fprintf(fp,"%.12g %.12g: %d : eke\n",eke,norm,index);
     double sum=0;
     for(int i=0; i<kxmax-kxmin+1; i++){
       fprintf(fp,"%d %.12g\n",i+kxmin,eke_i[i]);
       sum += eke_i[i];
     }//endfor
     fprintf(fp,"%.12g %.12g: %d : eke\n",eke,sum/2,index);
     fclose(fp);
    ***********************************/

//-----------------------------------------------------------------------------
   }// end routine : testeke
//==============================================================================



//==============================================================================
// RDMA routines for ibverbs CmiDirect optimization
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
#include <sstream>

void CP_State_GSpacePlane::completeRDMAhandshake(RDMASetupConfirmationMsg<RDMApair_GSP_PC> *msg)
{
	#ifndef PC_USE_RDMA
		CkPrintf("GSpace[%d,%d] aborting because someone called an RDMA method when it has been turned off\n",thisIndex.x,thisIndex.y);
		CkAbort("GSpace aborting RDMA setup completion");
	#else
	
	/// Retrieve the handshake token and the rdma handle from the message
	RDMApair_GSP_PC token = msg->token();
	rdmaHandleType ourHandle = msg->handle();
#ifdef DEBUG_CP_PAIRCALC_RDMA
        std::stringstream dbgStr; 
        dbgStr<<token;
		CkPrintf("%s : Received RDMA setup confirmation from paircalc. Now have %d handles of %d (%d symm + %d asymm)\n", dbgStr.str().c_str(),
			thisIndex.x,thisIndex.y, gotHandles+1, numRDMAlinksSymm+numRDMAlinksAsymm, numRDMAlinksSymm, numRDMAlinksAsymm );
#endif
	/// Determine which loop (symm/Asymm) this PC that has sent setup confirmation, belongs to
    cp::gspace::PCCommManager *pcMgr;
	if (token.symmetric)
        pcMgr  = &symmPCmgr;
	else
        pcMgr  = &asymmPCmgr;
		
	/// Compute the location and amount of data to be sent
	int chunkSize= gs.numPoints / pcMgr->pcCfg.numChunks;
	int offset   = token.pcIndex.z * chunkSize;
	int dataSize = chunkSize;
	/// The last chunk of data should get whatever is remaining
	if( (pcMgr->pcCfg.numChunks > 1) && (token.pcIndex.z == pcMgr->pcCfg.numChunks-1) )
		dataSize += gs.numPoints % pcMgr->pcCfg.numChunks;
		
	/// If the PC should be sent left matrix data ...
	if (token.shouldSendLeft)
	{
		/// Stuff the location of the data to be sent into the handle
		CmiDirect_assocLocalBuffer(&ourHandle,&(gs.packedPlaneData[offset]),dataSize*sizeof(complex));
		/// Store the rdma handle in the appropriate array
		pcMgr->leftDestinationHandles.push_back(ourHandle);
	}
	/// ... else if the PC should be sent right matrix data
	else
	{
		/// Stuff the location of the data to be sent into the handle
		if (token.symmetric)
			CmiDirect_assocLocalBuffer(&ourHandle,&(gs.packedPlaneData[offset]),dataSize*sizeof(complex));
		else
			CmiDirect_assocLocalBuffer(&ourHandle,&(gs.packedForceData[offset]),dataSize*sizeof(complex));
		/// Store the rdma handle in the appropriate array
		pcMgr->rightDestinationHandles.push_back(ourHandle);
	}
	
#ifdef DEBUG_CP_PAIRCALC_RDMA
        CkPrintf("%s : Will RDMA-put %d units of data at an offset of %d units from %p on proc %d to %p on proc %d\n",
			dbgStr.str().c_str(), dataSize, offset, 
            (!token.shouldSendLeft && !token.symmetric)? &(gs.packedForceData) : &(gs.packedPlaneData), ourHandle.senderNode,
            ourHandle.recverBuf,ourHandle.recverNode); 
#endif

	/// Call a reduction that signals the end of the initialization phase to main
	if(++gotHandles == numRDMAlinksSymm + numRDMAlinksAsymm)
	{
#ifdef DEBUG_CP_PAIRCALC_RDMA
			CkPrintf("GSpace[%d,%d] received RDMA setup confirmation from all %d PCs (%d symm + %d asymm) I was expecting. Triggering reduction to indicate end of init phase.\n",
				thisIndex.x,thisIndex.y,gotHandles,numRDMAlinksSymm,numRDMAlinksAsymm);
#endif
		int i=1;
		CkCallback cbDoneInit = CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy);
		contribute(sizeof(int), &i, CkReduction::sum_int, cbDoneInit, thisInstance.proxyOffset);
	}
	delete msg;
	#endif // PC_USE_RDMA
}

#include "gStatePlane.def.h"

