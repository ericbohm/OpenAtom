//======================================================
// Things to do : 
//    move resetiterstate
//======================================================

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
#include "charm++.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>
//---------------------------------------------------------------------------
#include "../../include/debug_flags.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"
#include "StructFactorCache.h"
//---------------------------------------------------------------------------
#define CHARM_ON
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpintegrate.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "../../src_piny_physics_v1.0/include/class_defs/allclass_gen.h"
#include "../../src_piny_physics_v1.0/include/class_defs/allclass_cp.h"
//============================================================================


//============================================================================
extern Config config;
extern PairCalcID pairCalcID1;
extern PairCalcID pairCalcID2;

extern CProxy_main mainProxy;
extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern CProxy_Ortho orthoProxy;
extern CProxy_CP_State_ParticlePlane particlePlaneProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_ComlibManager mgrProxy;
extern ComlibInstanceHandle mssInstance;
extern CProxy_StructureFactor sfCompProxy;

extern CProxy_EnergyGroup egroupProxy; //energy group proxy
extern int nstates;
extern int sizeX;
extern int nchareG;              // number of g-space chares <= sizeX and >=nplane_x
extern int atom_integrate_done;  // not a readonly global : a group of one element
extern ComlibInstanceHandle mcastInstancePP;
extern CProxy_FFTcache fftCacheProxy;
extern CProxy_StructFactCache sfCacheProxy;

void cleanExit(void *, void *);
void printMagForcePsi(void *, void *);
void printFictEke(void *, void *);
void printEnergyEke(void *, void *);
void testeke(int ,complex *,int *,int *,int *, int ,int);

//============================================================================


//============================================================================
// This routine defines the control flow of GSpacePlane
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
RTH_Routine_locals(CP_State_GSpacePlane,run)
  RTH_Routine_code(CP_State_GSpacePlane,run) {
//============================================================================
  while(1) {
 //===========================================================================
 // (I) Compute the forces and coef evolution code block

    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 1 || 
       c->first_step != 1){
    //------------------------------------------------------------------------
    // (A) Start new iteration : Reset counters
       c->startNewIter();
    //------------------------------------------------------------------------
    // (B) Start SF/computeZ, FFT psi(gx,gy,gz)->psi(gx,gy,z), Send psi to real
       /*       if(c->weneedPsiV()) 
	 {
	   RTH_Suspend(); // wait for all psiv to be done
	 }
       */
       c->releaseSFComputeZ();
       c->doFFT(); 
       c->sendFFTData();
       RTH_Suspend(); // wait for (psi*vks)=F[gx,gy,z] to arive from RealSpace
    //------------------------------------------------------------------------
    // (C) Complete IFFT of F(gx,gy,z) to F(gx,gy,gz)
       c->doIFFT();  // Message from realspace arrives : doifft(msg) resumes
    //------------------------------------------------------------------------
    // (D) Combine non-local and vks forces then compute eke forces
    //     If NL-pseudo forces done, completedExtExcNlForces calls combineForcesGetEke
       if (!(c->completedExtExcNlForces())){
         RTH_Suspend(); // If NL-pseudo forces are not finished then `suspend'.
                        // PP calls combineForcesGetEke() which invokes resume
       }//endif
#ifdef GIFFT_BARRIER
       if(!(c->allDoneIFFT())){
	  RTH_Suspend(); // wait for broadcast that all gspace is done  
       }//endif
#endif
    //------------------------------------------------------------------------
    // (E) Add contraint forces (rotate forces to non-orthogonal frame)
       c->sendLambda();        
       RTH_Suspend(); // wait for forces to be fixed up : acceptLambda resumes
    //------------------------------------------------------------------------
    // (F) If CG minimization : construct conjugate gradient
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_cg == 1){
	c->computeCgOverlap();
	RTH_Suspend(); // wait for cg reduction : psiCgOvlap resumes
      }// endif : CP-CG minimization
    //------------------------------------------------------------------------
    // (G) Output the states for cp dynamics
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0){
         c->writeStateDumpFile();// wait for output : psiwritecomplete resumes
  	 if(c->iwrite_now==1){RTH_Suspend();}
      }//endif
    //------------------------------------------------------------------------
    // (H) Evolve the electrons to the next step (Atom integration is hidden)
#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif
      c->integrateModForce();
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(IntegrateModForces_, StartTime, CmiWallTimer());    
#endif
    }// endif determine entry point

 //==========================================================================
 // (III) Orthogonalization and output code block
   //------------------------------------------------------------------------
   // (A) Orthogonalize
    c->sendPsi();      // send to Pair Calculator
    RTH_Suspend();     // Wait for new Psi : resume is called in acceptNewPsi
   //------------------------------------------------------------------------
   // (B) Output the states for minimization
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 1){
       c->writeStateDumpFile();// wait for output : psiwritecomplete resumes
       if(c->iwrite_now==1){RTH_Suspend();}
    }//endif
    if(c->weneedPsiV()) //ortho will have told us this
      {
	c->sendPsiV();
	RTH_Suspend();     // Wait for new PsiV : resume is called in acceptNewPsiV
      }
    c->first_step = 0; // its not the first step anymore!

  } //end while: Go back to top of loop (no suspending : no pausing)
//============================================================================

//--------------------------------------------------------------------------
  } RTH_Routine_end(CP_State_GSpacePlane,run)
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::gdoneIFFT(CkReductionMsg *msg){
      delete msg;
      //everybody lambda!
      allgdoneifft=true;
      RTH_Runtime_resume(run_thread);
  }
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::gdonePsiV(CkReductionMsg *msg){
      delete msg;
      //let my nonlocals go!
      needPsiV=false;
      RTH_Runtime_resume(run_thread);
  }
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::psiWriteComplete(CkReductionMsg *msg){
  delete msg;
  RTH_Runtime_resume(run_thread);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void printEnergyEke(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d = ((double *)m->getData())[0];
  delete m;

  CkPrintf("EKE         = %5.8lf\n", d);
  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_EKE, d);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void printFictEke(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d   = ((double *)m->getData())[0];
  double dd  = ((double *)m->getData())[1];
  delete m;

  int iopt   = (int) dd;
  if(iopt==0){CkPrintf("Fict Eke    =  %.10g\n", d);}
  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_FICTEKE, d);  

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void printMagForcePsi(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d = ((double *)m->getData())[0];
  delete m;

  CkPrintf("MagForPsi   =  %5.8lf\n", d);
  CkPrintf("Memory      =  %d\n",CmiMemoryUsage());
  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_FMAG, d);  

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void cleanExit(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  delete m;
  CkPrintf("======================================================\n\n\n");

  CkPrintf("======================================================\n");
  CkPrintf("         Open Atom Simulation Complete                \n");
  CkPrintf("======================================================\n");
  CkExit();
}
//============================================================================


//============================================================================
// start the thread that controls execution of GSpacePlane object
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::run () {
  run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_State_GSpacePlane,run),this);
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
// entry method to resume execution
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::resumeThread (PPDummyMsg *dmsg) {
  delete dmsg;
  RTH_Runtime_resume(run_thread);
}
//============================================================================

//============================================================================
// entry method to resume execution
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::resumePsiV (CkReductionMsg *msg) {
  delete msg;
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::psiCgOvlap(CkReductionMsg *msg){

  double d = ((double *)msg->getData())[0];
  fovlap = d;

  delete msg;  
  RTH_Runtime_resume(run_thread);

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::makePCproxies(){
  lambdaproxy=makeOneResultSection_asym(&gpairCalcID2, thisIndex.x, thisIndex.y);
  psiproxy=makeOneResultSection_sym1(&gpairCalcID1, thisIndex.x, thisIndex.y);
  if(AllExpected>1)
    psiproxyother=makeOneResultSection_sym2(&gpairCalcID1, thisIndex.x, thisIndex.y);
}
//============================================================================

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
                                           size2d size, 
                                           int    gSpaceUnits, 
                                           int    realSpaceUnits, 
                                           int    s_grain) 
//============================================================================
   {//begin routine
//============================================================================

//  ckout << "State G Space Constructor : "
//	<< thisIndex.x << " " << thisIndex.y << " " <<CkMyPe() << endl;

//============================================================================

  iwrite_now = 0;
  first_step = 1;
  iteration=0;
  partialCount=0;
  ireset=1;
  count = 0;
  gSpaceNumPoints = 0;

  ecount              = 0; //No energies have been received.
  allEnergiesReceived = 0;
  sendFFTDataSize= 0;

  total_energy   = 0.0;
  ehart_total    = 0.0;
  enl_total      = 0.0;
  eke_total      = 0.0;
  egga_total     = 0.0;
  eexc_total     = 0.0;
  eext_total     = 0.0;
  ewd_total      = 0.0;

  localState     = 0;
  displace_count = 0;
  fovlap         = 0.0;
  fovlap_old     = 0.0;

  ffttempdata    = NULL;
  k_x            = NULL;
  k_y            = NULL;
  k_z            = NULL;
  coef_mass      = NULL;
  int ourgrain=thisIndex.x/config.sGrainSize*config.sGrainSize; 
  if(nstates == config.sGrainSize)
    AllExpected=1;
  else if(ourgrain<(nstates-config.sGrainSize)) // corner has no extras
    AllExpected=2;
  else
    AllExpected=1;
  initialized   = false;
  doneDoingIFFT = false;
  allgdoneifft=false;
  acceptedPsi=true; // we start out with a psi
  needPsiV=false; // don't need tolerance check in first step
//============================================================================

  initGStateSlab(&gs, sizeX, size, gSpaceUnits, realSpaceUnits, 
		 s_grain, thisIndex.y,thisIndex.x);

//============================================================================

  usesAtSync = CmiTrue;
  flagsSent = false;
  if(config.lbgspace){
    setMigratable(true);
  }else{
    setMigratable(false);
  }//endif

  real_proxy = realSpacePlaneProxy;
  if (config.useCommlib)
      ComlibAssociateProxy(&mssInstance,real_proxy);

  // create structure factor proxy
  if(thisIndex.x==0){
      // dups must be less than the number of states because thats
      // the maximum number of time you can duplicate a plane
      int dups=config.numSfDups;
      if(config.numSfDups>nstates){dups=nstates;}
      CkVec <CkArrayIndex3D> sfelems;
      for(int dup=0;dup<dups; dup++){ //each dup
	  for(int atm=0;atm<config.numSfGrps; atm++){ //each atm
	      sfelems.push_back(CkArrayIndex3D(atm, thisIndex.y,dup));
	  }//endfor : atm groups
      }//endfor : dup groups
      sfCompSectionProxy = 
	CProxySection_StructureFactor::ckNew(sfCompProxy.ckGetArrayID(),
					     (CkArrayIndexMax *) 
					     sfelems.getVec(), sfelems.size());
  }//endif : state=0 

//---------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_GSpacePlane::CP_State_GSpacePlane(CkMigrateMessage *m) {
  run_thread = NULL;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_GSpacePlane::~CP_State_GSpacePlane(){
  if(initialized) {
    delete [] k_x;
    delete [] k_y;
    delete [] k_z;
    delete [] coef_mass;
    k_x = NULL;
    k_y = NULL;
    k_z = NULL;
    coef_mass = NULL;
  }//
}
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// In this function data is read from files, and sent to the corresponding
// G-space planes. Data reading will be done in chunk 0 of each state
//============================================================================
void CP_State_GSpacePlane::readFile() {
//============================================================================
// Local pointers

  int numData     = config.numData;
  int ibinary_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->ibinary_opt;
  CkVec <RunDescriptor> *sortedRunDescriptors=
                scProxy.ckLocalBranch()->cpcharmParaInfo->sortedRunDescriptors;
  int *npts_lgrp  = scProxy.ckLocalBranch()->cpcharmParaInfo->npts_per_chareG;
  int *nline_lgrp = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_per_chareG;
  int *istrt_lgrp = NULL;
  int *iend_lgrp  = NULL;

//============================================================================
// Set the file name using the config path and state number

  char fname[1024];
  int ind_state=thisIndex.x;
  sprintf(fname, "%s/state%d.out", config.dataPath, ind_state + 1);
  //------------------------------------------------------------------
  // Get the complex data, Psi(g) and the run descriptor (z-lines in g-space)

  complex *complexPoints = new complex[numData];
  int *kx=  new int[numData];
  int *ky=  new int[numData];
  int *kz=  new int[numData];
  int nlines_tot,nplane,nx,ny,nz;

  readState(numData,complexPoints,fname,ibinary_opt,&nlines_tot,&nplane, 
            kx,ky,kz,&nx,&ny,&nz,istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,0);

  if(config.low_x_size != nplane && config.doublePack){
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

      int runsToBeSent = sortedRunDescriptors[x].size();
      int numPoints    = 0;
      for (int j = 0; j < sortedRunDescriptors[x].size(); j++){
	  numPoints += sortedRunDescriptors[x][j].length;
      }//endfor

      complex *dataToBeSent  = new complex[numPoints];
      RunDescriptor *runDesc = new RunDescriptor[runsToBeSent];
      complex *temp          = complexPoints+ioff;
      CmiMemcpy(dataToBeSent,temp,(sizeof(complex) * numPoints));

      for (int j = 0; j < sortedRunDescriptors[x].size(); j++) {
	  runDesc[j] = sortedRunDescriptors[x][j];
      }//endfor
      if(ioff>numData){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Error reading\n");
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif

      gSpacePlaneProxy(ind_state, x).initGSpace(runsToBeSent,runDesc,
                                                numPoints,dataToBeSent,nx,ny,nz);
      delete [] dataToBeSent;
      delete [] runDesc;

      ioff += numPoints;
  }//endfor : loop over all possible chares in g-space (pencils)

  CkAssert(numData==ioff);

//============================================================================
// Clean up

  delete [] complexPoints;
  delete [] kx;
  delete [] ky;
  delete [] kz;

//---------------------------------------------------------------------------
   }//read the file
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// In this function data is written to files the simpliest way possible
//============================================================================
void CP_State_GSpacePlane::writeStateDumpFile() {
//============================================================================
// Local pointers and variables

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int cp_min_opt = (sim->cp_min_opt);
  int ndump_frq  = (sim->ndump_frq);
  int sizeX      = (sim->sizeX);
  int sizeY      = (sim->sizeY);
  int sizeZ      = (sim->sizeZ);
  int ind_state  = (thisIndex.x+1);
  int ind_chare  = (thisIndex.y+1);
  int ncoef      = gSpaceNumPoints;
  complex *psi   = gs.packedPlaneData;

//============================================================================
// Set the file names and write the files

  iwrite_now = 0;
  if(config.stateOutputOn==1){

    if( ((iteration % ndump_frq)==0) || (iteration==config.maxIter) ){
      complex *vpsi     = gs.packedPlaneDataScr;
      if(cp_min_opt==0){
        complex *forces   = gs.packedForceData;
        complex *vpsi_old = gs.packedVelData;
        memcpy(vpsi,vpsi_old,sizeof(complex)*ncoef);
        CPINTEGRATE::CP_integrate_half_vel(ncoef,iteration,forces,vpsi,psi,
                                           k_x,k_y,k_z,coef_mass);
      }//endif
      iwrite_now = 1;
      if(ind_state==1 && ind_chare==1){
        CkPrintf("-----------------------------------\n");
        CkPrintf("Writing states to disk on step %d\n",iteration);
        CkPrintf("-----------------------------------\n");
      }//endif
      char psiName[200]; char vpsiName[200];
      sprintf(psiName, "%s/newState%d_%d.out", config.dataPath,ind_state,ind_chare);
      sprintf(vpsiName,"%s/newVstate%d_%d.out",config.dataPath,ind_state,ind_chare);
      writePartState(ncoef,psi,vpsi,k_x,k_y,k_z,cp_min_opt,sizeX,sizeY,sizeZ,
                     psiName,vpsiName);
      int i=0;
      contribute(sizeof(int),&i,CkReduction::sum_int,
       CkCallback(CkIndex_CP_State_GSpacePlane::psiWriteComplete(NULL),gSpacePlaneProxy));
    }//endif : its time to write

  }//endif : it is useful to write
  
//---------------------------------------------------------------------------
   }//write the file
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
 *
 * The data is copied into the planes.
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::initGSpace(int            runDescSize, 
                                      RunDescriptor* runs, 
                                      int            size, 
                                      complex*       points,
                                      int nx, int ny, int nz) 
//============================================================================
   { //begin routine
//============================================================================

//#ifdef _CP_DEBUG_STATEG_VERBOSE_
//  CkPrintf("initGSpace %d.%d %d\n",thisIndex.x,thisIndex.y,size);
//#endif

  if (true == initialized) {
    ckerr << "Attempt to initialize a plane twice" << endl;
    return;
  }//endif
  initialized = true;

//============================================================================
// Setup gs

  gs.eke_ret  = 0.0;
    
  gs.numRuns     = runDescSize;
  gs.numLines    = runDescSize/2;
  gs.numFull     = (gs.numLines)*nz;
  gs.runs        = new RunDescriptor[gs.numRuns];
  gs.istate_ind  = thisIndex.x;
  gs.iplane_ind  = thisIndex.y;
  gs.mysizeX     = sizeX;
  gs.fftReqd     = true;
  gs.xdim        = 1;
  gs.ydim        = ny;
  gs.zdim        = nz;

  gs.numPoints=0;
  for (int r = 0; r < gs.numRuns; r++) {
    gs.numPoints += runs[r].length;
    gs.runs[r] = runs[r];
  }//endfor
  CkAssert(gs.numPoints == size);

  gs.packedPlaneData     = new complex[gs.numPoints];
  gs.packedPlaneDataScr  = new complex[gs.numPoints];
  gs.packedForceData     = new complex[gs.numPoints];
  gs.packedPlaneDataTemp = new complex[gs.numPoints];
  gs.packedVelData       = new complex[gs.numPoints];

  CmiMemcpy(gs.packedPlaneData, points, sizeof(complex)*gs.numPoints);
  CmiMemcpy(gs.packedPlaneDataTemp, points, sizeof(complex)*gs.numPoints);
  memset(gs.packedForceData, 0, sizeof(complex)*gs.numPoints);
  memset(gs.packedPlaneDataScr, 0, sizeof(complex)*gs.numPoints);
  memset(gs.packedVelData, 0, sizeof(complex)*gs.numPoints);

//============================================================================
// Setup gpspaceplane and particle plane

  CmiAssert(gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal());
  CmiAssert(particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal());
  particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal()->initKVectors(&gs);

//============================================================================
// Setup k-vectors and masses and zero the force overlap

  gs.setKVectors(&gSpaceNumPoints, &k_x, &k_y, &k_z);
  CkAssert(gSpaceNumPoints == size);

  coef_mass        = new double[gSpaceNumPoints];
  int mydoublePack = config.doublePack;
  CPINTEGRATE::CP_create_mass(gSpaceNumPoints,k_x,k_y,k_z,coef_mass,mydoublePack);

  fovlap      = 0.0; 

  complex *psi_g  = gs.packedPlaneData;
  complex *forces = gs.packedForceData;
  double *eke_ret = &(gs.eke_ret);
  int ncoef       = gSpaceNumPoints;

//  Use the kinetic energy to test the input
//  CPNONLOCAL::CP_eke_calc(ncoef,gs.istate_ind,forces,psi_g,k_x,k_y,k_z,eke_ret,
//			    config.doublePack);
//  CkPrintf("eke %g : %d %d \n",gs.eke_ret,gs.istate_ind,gs.iplane_ind);

  gs.eke_ret = 0;

//============================================================================
// Send the k's to the structure factor 

  for(int atm=0;atm<config.numSfGrps; atm++){ //each atm
    for(int dup=0;dup<config.numSfDups;dup++){ //each dup
      if(dup==thisIndex.x){
        sfCompProxy(atm,thisIndex.y,dup).acceptKVectors(gSpaceNumPoints, k_x, k_y, k_z);
      }//endif
    }//endfor
  }//endfor

//============================================================================
/* This reduction is done to signal the end of initialization */

  int i=1;
  contribute(sizeof(int), &i, CkReduction::sum_int, 
	     CkCallback(CkIndex_main::doneInit(0),mainProxy));

//============================================================================
//Some PC initialization that needs to happen here to avoid
//constructor race conditions

  gpairCalcID1=pairCalcID1;
  gpairCalcID2=pairCalcID2;
  makePCproxies();

//---------------------------------------------------------------------------
   }// end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CP_State_GSpacePlane::pup(PUP::er &p) {

//============================================================================

  ArrayElement2D::pup(p);
  p|needPsiV;
  p|doneDoingIFFT;
  p|allgdoneifft;
  p|initialized;
  p|iteration;
  p|ireset;
  p|count;
  gs.pup(p);
  p|flagsSent;
  p|partialCount;
  p|gSpaceNumPoints;
  p|real_proxy;
  // k_x, k_y, k_z need to be puped
  if (p.isUnpacking()) {
    k_x       = new int[gSpaceNumPoints];
    k_y       = new int[gSpaceNumPoints];
    k_z       = new int[gSpaceNumPoints];
    coef_mass = new double[gSpaceNumPoints];
  }//endif
    
  p(k_x, gSpaceNumPoints);
  p(k_y, gSpaceNumPoints);
  p(k_z, gSpaceNumPoints);
  p(coef_mass,gSpaceNumPoints);
  p|sendFFTDataSize;
  p|ehart_total;
  p|enl_total;
  p|eke_total;
  p|fmagPsi_total;
  p|fovlap;
  p|fovlap_old;
  p|egga_total;
  p|eexc_total;
  p|eext_total;
  p|ewd_total;
  p|ecount;
  p|displace_count;
  p|total_energy;
  p|allEnergiesReceived;
  p|localState;
  p|acceptedPsi;
  p | AllExpected;
  p | sfCompSectionProxy;
  p | gpairCalcID1;
  p | gpairCalcID2;
  p | iwrite_now;
    
//-------------------------------------------------------
   }// end routine : pup
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::startNewIter ()  {
//============================================================================
// Check for flow of control errors

    if(atom_integrate_done==0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Starting new iter before\n");
      CkPrintf("finishing atom integrate\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    } //endif

//============================================================================
// Reset all the counters etc.

  ecount              = 0;
  allEnergiesReceived = 0;
  total_energy        = 0.0;
  displace_count      = 0;
  doneDoingIFFT       = false;
  count               = 0;   // 'count' is used to check if all IFFT'd data 
                             // has arrived from RealSpacePlane
  doneDoingIFFT       = false;
  allgdoneifft        = false;

  CP_State_ParticlePlane *pp = 
                   particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();        
  pp->doneGettingForces = false;
  pp->doneEnl           = 0;
  pp->enl               = 0.0;
  pp->doneForces        = 0;


//============================================================================
// Check Load Balancing, Increment counter, set done flags equal to false.

#ifndef CMK_OPTIMIZE
  if(iteration==TRACE_ON_STEP ){traceBegin();}
  if(iteration==TRACE_OFF_STEP){traceEnd();}
#endif

    iteration++;
    if(config.lbgspace || config.lbpaircalc){
	if((iteration % (FIRST_BALANCE_STEP - PRE_BALANCE_STEP) == 0)  || 
           (iteration % (LOAD_BALANCE_STEP - PRE_BALANCE_STEP) == 0)){
	    LBTurnInstrumentOn();
	}//endif
    }//endif

//============================================================================
// Output psi at start of minimization for debugging

  int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  if(iteration==1 && cp_min_opt==1){screenOutputPsi();}

//---------------------------------------------------------------------------
    }//end routine
//============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::releaseSFComputeZ() {
//==============================================================================

  if(thisIndex.x==0){
       //multicast to all states of our plane and dups using the section proxy
	SFDummyMsg *msg = new(8*sizeof(int)) SFDummyMsg;
	CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	*(int*)CkPriorityPtr(msg) = config.sfpriority;
	sfCompSectionProxy.computeSF(msg);
    }//endif

//==============================================================================
//check all SFs
  CP_State_ParticlePlane *pp = particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  for(int i=0;i<config.numSfGrps;i++){
      if(pp->haveSFAtmGrp[i]>=0){
	    PPDummyMsg *pmsg = new (8*sizeof(int)) PPDummyMsg;
	    pmsg->atmGrp=i;
	    pmsg->sfindex=pp->haveSFAtmGrp[i];
	    CkSetQueueing(pmsg, CK_QUEUEING_IFIFO);
	    *(int*)CkPriorityPtr(pmsg) = config.sfpriority+i+config.numSfGrps; 
            //lower than sf and sfcache
	    particlePlaneProxy(thisIndex.x, thisIndex.y).computeZ(pmsg);
      }//endif
  }//endfor
  
//----------------------------------------------------------------------------
  }//end routine
//==============================================================================


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

  // If there is no data to send, return immediately
  if (gs.numNonZeroPlanes == 0 || gs.fftReqd == false){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude, no data to send : Why launch the stategpsaceplane\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// destroys forces
// creates new ones
// copies test states to the working states

#ifndef CMK_OPTIMIZE    
  double StartTime=CmiWallTimer();
#endif

#ifdef  _CP_DEBUG_UPDATE_OFF_
  CmiMemcpy(gs.packedPlaneData,gs.packedPlaneDataTemp,
            sizeof(complex)*gs.numPoints);
  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
    memset(gs.packedVelData,0,sizeof(complex)*gs.numPoints);
  }//endif
#endif      

// Do fft in forward direction, 1-D, in z direction
// A local function not a message

  ffttempdata = fftCacheProxy.ckLocalBranch()->doGSRealFwFFT(gs.packedPlaneData, 
   	        gs.runs, gs.numRuns, gs.numLines, gs.numFull, gs.numPoints,
                gs.zdim, gs.fftReqd);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(GspaceFwFFT_, StartTime, CmiWallTimer());
#endif   
 
}
//============================================================================


//============================================================================
// Send result to realSpacePlane : perform the transpose
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CP_State_GSpacePlane::sendFFTData () {

//============================================================================
// Do a Comlib Dance

  if (config.useCommlib){mssInstance.beginIteration();}

//============================================================================
// Send your (x,y,z) to processors z.

  int numLines = gs.numLines; // same amount of data to each realspace chare puppy
  int sizeZ    = gs.planeSize[1];

  for(int z=0; z < sizeZ; z++) {

    RSFFTMsg *msg    = new (numLines,8*sizeof(int)) RSFFTMsg;
    msg->size        = numLines;
    msg->senderIndex = thisIndex.y;  // planenumber
    msg->numPlanes   = gs.numNonZeroPlanes; // unity baby
    
    if(config.prioFFTMsg){
       CkSetQueueing(msg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(msg) = config.rsfftpriority + 
                                   thisIndex.x*gs.planeSize[0]+thisIndex.y;
    }//endif
    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = ffttempdata[j];}
    real_proxy(thisIndex.x, z).doFFT(msg);  // same state, realspace char[z]

  }//endfor
    
  if (config.useCommlib){mssInstance.endIteration();}
    
  ffttempdata = NULL; // its memory from the group : don't touch it

//----------------------------------------------------------------------
  }//end routine 
//============================================================================


//============================================================================
/**
 * This is used to send data to the CP_State_GSpacePlanes, which  do the
 * inverse ffts (upon receiving data from all the corresponding
 * realSpacePlanes)
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::doIFFT(GSIFFTMsg *msg) {
//============================================================================

  int size             = msg->size;
  int offset           = msg->offset;
  complex *partlyIFFTd = msg->data;

  int numLines         = gs.numLines;
  int sizeZ            = gs.planeSize[1];
  int expandedDataSize = numLines*sizeZ;

  CkAssert(numLines == size);

  if(ffttempdata==NULL) {
    ffttempdata = (complex *)fftw_malloc(expandedDataSize *sizeof(complex));
    memset(ffttempdata, 0, sizeof(complex)*expandedDataSize);
  }//endif

  // z=offset is inner index : collections of z-lines of constant (gx,gy)
  for(int i=0,j=offset; i< numLines; i++,j+=sizeZ){ffttempdata[j] = partlyIFFTd[i];}

  delete msg;
  // receive 1 message from each z-chare    
  count++;
  if (count == gs.planeSize[1]) {
    count = 0;
#ifdef GIFFT_BARRIER
    //put contribute here to reduction with a broadcast client
    int wehaveours=1;
    contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
      CkCallback(CkIndex_CP_State_GSpacePlane::gdoneIFFT(NULL),gSpacePlaneProxy));
#endif
    RTH_Runtime_resume(run_thread);
  }//endif : has everyone arrived?

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::doIFFT () {
//============================================================================
//now do the IFFT  does the ifft in place (this is done on the forceArr)

#ifndef CMK_OPTIMIZE    
  double StartTime=CmiWallTimer();
#endif

  gs.doBwFFT(ffttempdata);
  fftw_free(ffttempdata);
  ffttempdata = NULL;
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(GspaceBwFFT_, StartTime, CmiWallTimer());
#endif    

  ffttempdata = NULL;
  /* Scale the forces by 1/(128^3) */
  double scaleFactor = 1/double(scProxy.ckLocalBranch()->cpcharmParaInfo->sizeX * 
                                scProxy.ckLocalBranch()->cpcharmParaInfo->sizeY * 
                                scProxy.ckLocalBranch()->cpcharmParaInfo->sizeZ);
//============================================================================
// set the forces from vks

  complex *forces = gs.packedForceData;
  for(int index=0; index<gs.numPoints; index++){
    forces[index].re = -2.0 * forces[index].re * scaleFactor;
    forces[index].im = -2.0 * forces[index].im * scaleFactor;
  }/*endfor*/

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
/**
 * This function calls the force integration function if non-local pseudo
 * computation is complete
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
bool CP_State_GSpacePlane::completedExtExcNlForces () {
//============================================================================
// Test if the Non-local pseudo computation part is done

  CP_State_ParticlePlane *pp = 
                  particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();        

  if (pp->doneGettingForces) { 
    CkAssert(!doneDoingIFFT); //if doneDoingIFFT were already true,
			      // we shouldn't be here so exit
    combineForcesGetEke();
    doneDoingIFFT = true;
    return true;
  }else {
    doneDoingIFFT = true;
    return false;
  }//endif

//-----------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::combineForcesGetEke(){
//================================================================================
// add forces from particle plane to forces from IFFT then zero them

  CP_State_ParticlePlane *pp = 
    particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();

#ifdef _CP_DEBUG_VKS_FORCES_
  if(thisIndex.x==0 && thisIndex.y==0){
    FILE *fp = fopen("vks_forces_s0_p0.out","w");
    int ncoef       = gs.numPoints;
    complex *forces = gs.packedForceData;
    for(int i=0;i<ncoef;i++){
      fprintf(fp,"%d %d %d : %g %g\n",pp->k_x[i],pp->k_y[i],pp->k_z[i],
	      forces[i].re,forces[i].im);
    }//endfor
    fclose(fp);
  }//endif
#endif

  gs.addForces(pp->myForces,pp->k_x);
  bzero(pp->myForces,gs.numPoints*sizeof(complex));

//================================================================================
// Compute force due to quantum kinetic energy and add it in.
// Reduce quantum kinetic energy or eke

  int istate      = gs.istate_ind;
  int ncoef       = gs.numPoints;
  complex *psi_g  = gs.packedPlaneData;
  complex *forces = gs.packedForceData;
  double *eke_ret = &(gs.eke_ret);
  CPNONLOCAL::CP_eke_calc(ncoef,istate,forces,psi_g,k_x,k_y,k_z,eke_ret,
			  config.doublePack);
  contribute(sizeof(double), &gs.eke_ret, CkReduction::sum_double, 
	     CkCallback(printEnergyEke, NULL));

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

//========================================================================
// If called from ParticlePlane, then resume the GSpacePlane computation
// If called from here, you aren't suspended so don't resume.

  if (doneDoingIFFT){
    RTH_Runtime_resume(run_thread);
  }//endif

//-----------------------------------------------------------------------------
  }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendLambda() {
//==============================================================================

  complex *psi   = gs.packedPlaneData;
  complex *force = gs.packedForceData;
  int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  if(gs.ihave_kx0==1 && cp_min_opt==0){
    double rad2i = 1.0/sqrt(2.0);
    double rad2  = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i]   *= rad2i;}
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){force[i] *= rad2;}
  }//endif

  if(config.gSpaceNumChunks!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, while gSpaceNumChunk!=1 is cool, I'm \n");
    CkPrintf("afraid you'll have to do the implementation yourself!\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  int numPoints   = gs.numPoints / config.gSpaceNumChunks;
  int dataCovered = gs.numPoints;
  int c = 0;
  int toSend = (c == config.gSpaceNumChunks - 1) ? dataCovered : numPoints;

  startPairCalcLeft(&gpairCalcID2, toSend, psi + c * numPoints, 
		    thisIndex.x, thisIndex.y, false);
  startPairCalcRight(&gpairCalcID2, toSend, force + c * numPoints, 
		     thisIndex.x, thisIndex.y);

//-----------------------------------------------------------------------------
   }// end routine 
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptLambda(CkReductionMsg *msg) {
//==============================================================================

  complex *data = (complex *)msg->getData();
  complex *force = gs.packedForceData;
  int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;

//==============================================================================
// Get the modified forces

  if(config.doublePack==1){
   if(cp_min_opt==1){
    for(int i=0; i<gs.numPoints; i++){
       double wght  = (k_x[i]==0 ? 0.5 : 1);
       force[i].re -= wght*data[i].re;
       force[i].im -= wght*data[i].im;
     }//endfor
   }else{
     for(int i=0; i<gs.numPoints; i++){force[i]  = data[i]*(-1.0);}
     if(gs.ihave_kx0==1){
       double rad2i = 1.0/sqrt(2.0);
       for(int i=gs.kx0_strt; i<gs.kx0_end; i++){force[i] *= rad2i;}
     }//endif
   }//endif
  }else{
    for(int i=0; i<gs.numPoints; i++){
       force[i].re -= 0.5*data[i].re;
       force[i].im -= 0.5*data[i].im;
    }//endfor
  }//endif
  delete msg;  

//==============================================================================
// Compute the mag square of the forces

  double force_sq_sum=0.0;
  for(int i=0; i<gs.numPoints; i++){
    force_sq_sum+= force[i].getMagSqr();
  }//endfor
  contribute(sizeof(double), &force_sq_sum, CkReduction::sum_double, 
	     CkCallback(printMagForcePsi, NULL));
  
//==============================================================================
// Retrieve Non-orthog psi

  if(cp_min_opt==0){
     int ncoef          = gs.numPoints;
     complex *psi_g     = gs.packedPlaneData;
     complex *psi_g_scr = gs.packedPlaneDataScr;
     memcpy(psi_g,psi_g_scr,sizeof(complex)*ncoef);
  }//endif

//==============================================================================
// Resume 

  RTH_Runtime_resume(run_thread);

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

    fovlap_old = fovlap;  // every state gspaceplane[] has the same value

    double fovlap_loc = 0.0;
    double fovlap_now;
    int istate      = gs.istate_ind;
    int ncoef       = gs.numPoints;
    complex *psi_g  = gs.packedPlaneData;
    complex *forces = gs.packedForceData;
    CPINTEGRATE::CP_fovlap_calc(ncoef,istate,forces,&fovlap_now);
    fovlap_loc += fovlap_now;

    gs.fovlap_loc = fovlap_loc;  // save a local copy of your overlap
    contribute(sizeof(double),&fovlap_loc,CkReduction::sum_double,
	     CkCallback(CkIndex_CP_State_GSpacePlane::psiCgOvlap(NULL),thisProxy));

//----------------------------------------------------------------------------
}// end routine : computeCgOverlap
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::integrateModForce() {
//==============================================================================
    
  int cp_min_opt  = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  int cp_min_cg   = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_cg;

//==========================================================================
// (I) Set up the input to the integrator

  double gamma_conj_grad = 0.0;
  if( (cp_min_opt==1) && (cp_min_cg == 1) && (ireset==0)){
     gamma_conj_grad = fovlap/fovlap_old;
  }//endif
  ireset=0;

  int istate         = gs.istate_ind;
  int ncoef          = gs.numPoints;
  complex *psi_g     = gs.packedPlaneData; 
  complex *forces    = gs.packedForceData; 
  complex *vpsi_g    = gs.packedVelData; // for cp not minimization
  complex *forcesold = gs.packedVelData; // for miniziation not cp

//==========================================================================
// (II) Evolve the states using the forces/conjugate direction

  //---------------------------------------------------------------
  // (A) Debug output before integration
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
  gs.fictEke_ret = 0.0;
  CPINTEGRATE::CP_integrate(ncoef,istate,iteration,forces,forcesold,psi_g,
        k_x, k_y, k_z,coef_mass,gamma_conj_grad,&gs.fictEke_ret);
  //---------------------------------------------------------------
  // (C) Debug output after integration
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
#endif

//==========================================================================
// Contribute FictKe : screen output non-orthogonal psi

  double sendme[2];
  sendme[0] = gs.fictEke_ret;
  sendme[1] = (double)cp_min_opt;
  contribute(2*sizeof(double),sendme,CkReduction::sum_double, 
             CkCallback(printFictEke, NULL));

//==========================================================================
// Launch the atoms 

  if(thisIndex.x==0 && thisIndex.y==0){
    atomsGrpProxy.StartRealspaceForces();
  }//endif

//------------------------------------------------------------------------------
   } // end CP_State_GSpacePlane::integrateModForce
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendPsi() {
//==============================================================================
// Error checking

  acceptedPsi =false;
  if(config.gSpaceNumChunks!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, while gSpaceNumChunk!=1 is cool, I'm \n");
    CkPrintf("afraid you'll have to do the implementation yourself!\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//==============================================================================
// Prepare the data : If cp dynamics is going, save the non-orthogonal puppies.

  complex *data=gs.packedPlaneData;
  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
     int ncoef     = gs.numPoints;
     complex *scr  = gs.packedPlaneDataScr; //save non-orthog psi
     memcpy(scr,data,sizeof(complex)*ncoef);
  }//endif

  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){data[i] *= rad2i;}
  }//endif

//==============================================================================

  int s, c;
  int idx = (thisIndex.x/gs.S_grainSize) * gs.S_grainSize;
  int numPoints = gs.numPoints / config.gSpaceNumChunks;
  int dataCovered = gs.numPoints;

  c = 0;
  int toSend = (c == config.gSpaceNumChunks - 1) ? dataCovered : numPoints;

  startPairCalcLeft(&gpairCalcID1, toSend, data + c * numPoints, 
		    thisIndex.x, thisIndex.y, false);

//----------------------------------------------------------------------------
    }// end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendPsiV() {
//==============================================================================

  if(config.gSpaceNumChunks!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, while gSpaceNumChunk!=1 is cool, I'm \n");
    CkPrintf("afraid you'll have to do the implementation yourself!\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  complex *data=gs.packedVelData;
  /*  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
     int ncoef     = gs.numPoints;
     complex *scr  = gs.packedPlaneDataScr; //save non-orthog psi
     memcpy(scr,data,sizeof(complex)*ncoef);
  }//endif
  */

  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){data[i] *= rad2i;}
  }//endif

  int s, c;
  int idx = (thisIndex.x/gs.S_grainSize) * gs.S_grainSize;
  int numPoints = gs.numPoints / config.gSpaceNumChunks;
  int dataCovered = gs.numPoints;

  c = 0;
  int toSend = (c == config.gSpaceNumChunks - 1) ? dataCovered : numPoints;

  startPairCalcLeft(&gpairCalcID1, toSend, data + c * numPoints, 
		    thisIndex.x, thisIndex.y, true);

//----------------------------------------------------------------------------
}// end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNewPsi(CkReductionMsg *msg){
//=============================================================================
// (I) Unpack the contribution to newpsi (orthonormal psi)

  int N         = msg->getSize()/sizeof(complex);
  complex *data = (complex *)msg->getData();
  complex *psi  = gs.packedPlaneData;
  complex *vpsi = gs.packedVelData;

  if(partialCount<1){
    for(int i=0; i<N; i++){psi[i] = data[i];}
  }else{
    for(int i=0; i<N; i++){psi[i] += data[i];}
  }//endif
  delete msg;

//=============================================================================
// (II) If you have got it all : Rescale it, produce some output

  partialCount++;//psi arrives in as many as 2 reductions
  if(partialCount==AllExpected){ 
    //--------------------------------------------------------------------
    // (A) Reset counters and rescale the kx=0 stuff
    acceptedPsi=true;
    partialCount=0;
    if(gs.ihave_kx0==1){
      double rad2 = sqrt(2.0);
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i] *= rad2;}
    }//endif
    //--------------------------------------------------------------------
    // (B) Generate some screen output of orthogonal psi
    screenOutputPsi();
    //--------------------------------------------------------------------
    // (D) Go back to the top or exit
    if(iteration==config.maxIter){
      int i;
      contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
    }//endif
    RTH_Runtime_resume(run_thread);
  } // if partialCount

//----------------------------------------------------------------------------
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNewPsiV(CkReductionMsg *msg){
//=============================================================================
// (I) Unpack the contribution to newpsi (orthonormal psi)

  int N         = msg->getSize()/sizeof(complex);
  complex *data = (complex *)msg->getData();
  complex *vpsi = gs.packedVelData;

  if(partialCount<1){
    for(int i=0; i<N; i++){vpsi[i] = data[i];}
  }else{
    for(int i=0; i<N; i++){vpsi[i] += data[i];}
  }//endif
  
  delete msg;

//=============================================================================
// (II) If you have got it all : Rescale it, produce some output

  partialCount++;//psi arrives in as many as 2 reductions

  if(partialCount==AllExpected){ 
    //--------------------------------------------------------------------
    // (A) Reset counters and rescale the kx=0 stuff
    partialCount=0;
    if(gs.ihave_kx0==1){
      double rad2 = sqrt(2.0);
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){vpsi[i] *= rad2;}
    }//endif
    /* debugging barrier
    int wehaveours=1;
    contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
      CkCallback(CkIndex_CP_State_GSpacePlane::gdonePsiV(NULL),gSpacePlaneProxy));
    */
    needPsiV=false;
    RTH_Runtime_resume(run_thread);
  } // if partialCount

//----------------------------------------------------------------------------
}//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::screenOutputPsi(){
//==============================================================================

  int cp_min_opt  = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  int nstates     = scProxy.ckLocalBranch()->cpcharmParaInfo->nstates;
  complex *vpsi = gs.packedVelData;
  complex *fpsi = gs.packedForceData;
  complex *psi  = gs.packedPlaneData; //orthogonal psi
  if(cp_min_opt==0){psi=gs.packedPlaneDataScr;} //non-orthogonal psi

#ifdef _CP_DEBUG_COEF_SCREEN_
    if(gs.istate_ind==0 || gs.istate_ind==nstates-1){
      for(int i = 0; i < gs.numPoints; i++){
	if(k_x[i]==0 && k_y[i]==1 && k_z[i]==4 ){
	  CkPrintf("------------------------------------------------------\n");
	  CkPrintf(" Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],psi[i].re,psi[i].im);
          if(cp_min_opt==0){
            double vre=vpsi[i].re;
            double vim=vpsi[i].im;
            if(iteration>1){
	    }
 	    CkPrintf("VPsi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],vre,vim);
	  }//endif
	  CkPrintf("------------------------------------------------------\n");
	}//endif
	if(k_x[i]==2 && k_y[i]==1 && k_z[i]==3){
          double vre=vpsi[i].re;
          double vim=vpsi[i].im;
	  CkPrintf("------------------------------------------------------\n");
	  CkPrintf(" Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],psi[i].re,psi[i].im);
          if(cp_min_opt==0){
 	    CkPrintf("VPsi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],vre,vim);
	  }//endif
	  CkPrintf("------------------------------------------------------\n");
	}//endif
      }//endfor
    }//endif
    //--------------------------------------------------------------------
#endif

//----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::isAtSync(int numIter) {
//==============================================================================

#ifndef CMK_OPTIMIZE
    if(numIter == 2 * LOAD_BALANCE_STEP)
      traceBegin();
    if(numIter == 3 * LOAD_BALANCE_STEP)
      traceEnd();
#endif
    /* needs to be rewritten completely
    if(config.lbgspace)
      {

	StructFactCache *sfCache = sfCacheProxy.ckLocalBranch();
	CkAssert(sfCache != NULL);
	sfCache->removeAll();

      }
    */
    //    CkPrintf("G %d %d atsync\n",thisIndex.x, thisIndex.y);

    if(thisIndex.x==0 && thisIndex.y==0){
      isAtSyncPairCalc(&gpairCalcID1);
      isAtSyncPairCalc(&gpairCalcID2);
    }//endif
    AtSync();

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::ResumeFromSync() {
//==============================================================================
//    CmiPrintf("ResumeFromSync calls resume\n");

  // reset commlib proxies
  if(config.useCommlib)
      ComlibResetProxy(&real_proxy);

  // reset lambda PC proxies
  CkMulticastMgr *mcastGrp = 
        CProxy_CkMulticastMgr(gpairCalcID2.mCastGrpId).ckLocalBranch();         
  mcastGrp->resetSection(lambdaproxy);
  setResultProxy(&lambdaproxy, thisIndex.x, gpairCalcID2.GrainSize, 
          gpairCalcID2.mCastGrpId,true,CkCallback(CkIndex_Ortho::lbresume(NULL),
          orthoProxy));
  gpairCalcID2.resetProxy();

  LBTurnInstrumentOff();
//  CmiPrintf("G ResumeFromSync %d %d!\n",thisIndex.x, thisIndex.y);

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::syncpsi(){
//==============================================================================

  CkMulticastMgr *mcastGrp = 
     CProxy_CkMulticastMgr(gpairCalcID2.mCastGrpId).ckLocalBranch();         

  mcastGrp->resetSection(psiproxy);
  setResultProxy(&psiproxy,thisIndex.x,gpairCalcID1.GrainSize, gpairCalcID1.mCastGrpId, 
                  true, CkCallback(CkIndex_Ortho::lbresume(NULL),orthoProxy));

  if(AllExpected>1){
      mcastGrp->resetSection(psiproxyother);
      setResultProxy(&psiproxyother, thisIndex.x, gpairCalcID1.GrainSize, 
                   gpairCalcID1.mCastGrpId,true,CkCallback(CkIndex_Ortho::lbresume(NULL),
                   orthoProxy));
  }//endif

  //takes care of the paircalc result proxies which we own via the pairCalcID
  gpairCalcID1.resetProxy();

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::computeEnergies(int param, double d){
//==============================================================================

  switch(param){
      
    case ENERGY_EHART : 
      ehart_total = d;
      total_energy += d;
      ecount++;
      break;
      
    case ENERGY_ENL :
      enl_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EKE : 
      eke_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EGGA :
      egga_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EEXC :
      eexc_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EEXT :
      eext_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_EWD :
      ewd_total = d;
      total_energy += d;
      ecount++;
      break;
        
    case ENERGY_FICTEKE : 
      fictEke_total = d;
      total_energy += d;
      ecount++;
      break;
      
    case ENERGY_FMAG : 
      fmagPsi_total = d;
      ecount++;
      break;
      
    default :
      CkAbort("unknown energy");
      break;
  }//end switch

//==============================================================================

  if(ecount == NUM_ENERGIES){
    allEnergiesReceived = 1;        

    EnergyStruct estruct;
    estruct.enl          = enl_total;
    estruct.eke          = eke_total;
    estruct.eext         = eext_total;
    estruct.ehart        = ehart_total;
    estruct.eewald_recip = ewd_total;
    estruct.egga         = egga_total;
    estruct.eexc         = eexc_total;
    estruct.fictEke      = fictEke_total;
    estruct.totalEnergy  = total_energy;
    estruct.fmagPsi      = fmagPsi_total;

    egroupProxy.updateEnergies(estruct); // broadcast the complete set of energies
                                         // to the group so that all procs have them.

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

   for(int i = 0; i < ncoef; i++){
     
     gx = tpi*(k_x[i]*hmati[1] + k_y[i]*hmati[2] + k_z[i]*hmati[3]);
     gy = tpi*(k_x[i]*hmati[4] + k_y[i]*hmati[5] + k_z[i]*hmati[6]);
     gz = tpi*(k_x[i]*hmati[7] + k_y[i]*hmati[8] + k_z[i]*hmati[9]);
     g2 = gx*gx + gy*gy + gz*gz;

     if(g2<=ecut){
       double wght_now = 2.0;
       if(k_x[i]==0 && k_y[i]<0){wght_now=0.0;}
       if(k_x[i]==0 && k_y[i]==0 && k_z[i]<0){wght_now=0.0;}
       if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){wght_now=1.0;}
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
   CkPrintf("%.12g %.12g %.12g %.12g: %d : true eke\n",eke,eke2,norm,norm2,index);

//-----------------------------------------------------------------------------
   }// end routine : testeke
//==============================================================================

//==============================================================================
// Deprecated routine
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptAllLambda(CkReductionMsg *msg) {
    delete msg;
    CkAbort("GSP do not call acceptAllLambda\n");
}// end routine
//==============================================================================

void CP_State_GSpacePlane::requirePsiV() {
  needPsiV=true;
  int foo=1;
  // when everyone is ready, restart Ortho's backward path
  contribute(sizeof(int), &foo, CkReduction::min_int, CkCallback(CkIndex_Ortho::resumeV(NULL), orthoProxy));
}

bool CP_State_GSpacePlane::weneedPsiV() { return needPsiV;}


