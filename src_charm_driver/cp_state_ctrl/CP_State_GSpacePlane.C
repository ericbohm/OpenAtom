//======================================================
// Things to do : 
//    move resetiterstate
//======================================================

//#define _CP_DEBUG_WARN_SUSPEND_
//#define _CP_DEBUG_NONLOC_BARRIER_
//#define _CP_DEBUG_ORTHO_OFF_


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
#include "cpaimd.h"
#include "groups.h"
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
extern ComlibInstanceHandle mcastInstancePP;
extern CProxy_FFTcache fftCacheProxy;
extern CProxy_StructFactCache sfCacheProxy;

void cleanExit(void *, void *);
void allDoneCPForces(void *, void *);
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
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==1 || c->first_step!=1){
    //------------------------------------------------------------------------
    // (A) Start new iteration : Reset counters
#ifdef _CP_DEBUG_PSI_OFF_  // only move the atoms 
       c->myenergy_reduc_flag  =0; 
       c->myatom_integrate_flag=0; 
#else                     // move everything
       c->startNewIter();
    //------------------------------------------------------------------------
    // (B) Start SF/computeZ of NL forces. If you want, wait until they finish
#ifndef _CP_DEBUG_SFNL_OFF_ // non-local is allowed
       c->releaseSFComputeZ();
#ifdef _CP_DEBUG_NONLOC_BARRIER_
       RTH_Suspend();  // resume is called by acceptNLForces()
#endif
#endif  // non-local is allowed
    //------------------------------------------------------------------------
    // (C) FFT psi(gx,gy,gz)->psi(gx,gy,z), Send psi to real, 
    //      complete IFFT of F(gx,gy,z) to F(gx,gy,gz) when it comes back.
       c->doFFT(); 
       c->sendFFTData();
       RTH_Suspend();  // wait for (psi*vks)=F[gx,gy,z] to arive from RealSpace
       c->doIFFT();    // Message from realspace arrives : doifft(msg) resumes
    //------------------------------------------------------------------------
    // (D) If NL-pseudo forces are not done, wait for them.
#ifndef _CP_DEBUG_SFNL_OFF_ // non-local is allowed
       if(!c->doneNLForces()){
         RTH_Suspend(); // resume called in acceptNLforces when NL forces done
       }//endif
#endif // non-local is allowed
    //------------------------------------------------------------------------
    // (E) If NL-pseudo forces done, suspend.
       c->combineForcesGetEke();
#ifdef GIFFT_BARRIER  // pause for every single chare to finish
       if(!(c->allDoneIFFT())){
	  RTH_Suspend(); // wait for broadcast that all gspace is done  
       }//endif
#endif               //end pause
#endif // only move the atoms
    //------------------------------------------------------------------------
    // E) The atoms can't go until all cp forces are accounted for.
    //    However, the atoms can overlap with all this lambda, psi stuff.
       c->launchAtoms();
    //------------------------------------------------------------------------
    // (F) Add contraint forces (rotate forces to non-orthogonal frame)
#ifndef _CP_DEBUG_PSI_OFF_  // you are moving evertying
       c->sendLambda();        
#ifndef _CP_DEBUG_ORTHO_OFF_
       RTH_Suspend(); // wait for forces to be fixed up : acceptLambda resumes
#endif
    //------------------------------------------------------------------------
    // Get sum sq forces (even under dynamics its good to have) : also cg thingy
	c->computeCgOverlap();
	RTH_Suspend(); // wait for cg reduction : psiCgOvlap resumes
    //------------------------------------------------------------------------
    // (H) Output the states for cp dynamics
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0){
         c->writeStateDumpFile();// wait for output : psiwritecomplete resumes
  	 if(c->iwrite_now==1){RTH_Suspend();}
      }//endif
    //------------------------------------------------------------------------
    // (I) Evolve the electrons to the next step
      c->integrateModForce();
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0){
        c->sendRedPsi();  // Sync Redundant psi entries
        if(c->finishedRedPsi==0){
          RTH_Suspend();  // Resume is called in acceptRedPsi
	}//endif
        c->doneRedPsiIntegrate(); // after integrate AND acceptRedPsi
      }//endif
#endif  // you are moving everyting
    }// endif determine entry point
 //==========================================================================
 // (II) Orthogonalization and output code block
   //------------------------------------------------------------------------
   // (A) Orthogonalize
#ifndef _CP_DEBUG_PSI_OFF_   // move everything
    c->sendPsi();      // send to Pair Calculator
#ifndef _CP_DEBUG_ORTHO_OFF_
    RTH_Suspend();     // Wait for new Psi : resume is called in acceptNewPsi
#endif
   //------------------------------------------------------------------------
   // (B) Output the states for minimization
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 1){
       c->writeStateDumpFile();
       if(c->iwrite_now==1){RTH_Suspend();}// wait : psiwritecomplete resumes
    }//endif
   //------------------------------------------------------------------------
   // (C) Norb rotation
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0){
      if(c->weneedPsiV()){ //ortho will have told us this
         c->sendPsiV();
         RTH_Suspend();    // Wait for new PsiV : resume is called in acceptNewPsiV
      }//endif
    }//endif
   //------------------------------------------------------------------------
   // (d) Check for Energy reduction completion : should just be a safety
    if(c->myenergy_reduc_flag==0 && c->iteration>0){
      c->isuspend_energy=1;
#ifdef _CP_DEBUG_WARN_SUSPEND_
      CkPrintf("Suspend energy on proc %d : too many of these==bad scheme\n",CkMyPe());
#endif
      RTH_Suspend();     // resume called in acceptEnergy
    }//endif
#endif // you are moving everyting
   //------------------------------------------------------------------------
   // (E) Check for atom integration : should just be a safety
    if(c->myatom_integrate_flag==0 && c->iteration>0){
#ifdef _CP_DEBUG_WARN_SUSPEND_
      CkPrintf("Suspend atms on proc %d : too many of these==bad scheme\n",CkMyPe());
#endif
      c->isuspend_atms=1;
      RTH_Suspend();     // resume called in acceptAtoms
    }//endif
    c->first_step = 0;   // its not the first step anymore!
//--------------------------------------------------------------------------
   }//end while: Go back to top of the loop now (no suspending : no pausing)
//============================================================================

//--------------------------------------------------------------------------
  } RTH_Routine_end(CP_State_GSpacePlane,run)
//============================================================================


//============================================================================
//      Start the thread that controls execution of GSpacePlane object.
// It is invoked by main::doneInit which is a reduction client invoked by all
// (state,plane) objects after completion of CP_State_GspacePlane::initGSpace
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::run () {
  run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_State_GSpacePlane,run),this);
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
//     All GSpace objects finished back FFT : For debugging only
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
//     All GSpace objects have finished velocity rotation : Debugging only
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
//     All GSpace objects have finished writing coefs : NECESSARY
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::psiWriteComplete(CkReductionMsg *msg){
  delete msg;
  RTH_Runtime_resume(run_thread);
}
//============================================================================

//============================================================================
// Print out Quantum KE and put all the energies into the message group 
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void printEnergyEke(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d = ((double *)m->getData())[0];
  delete m;
#ifdef _CP_DEBUG_SFNL_OFF_
  CkPrintf("ENL         = OFF FOR DEBUGGING\n");
#endif
  CkPrintf("EKE         = %5.8lf\n", d);
  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_EKE, d);
}
//============================================================================


//============================================================================
// Print out Fict CP KE and send it to the energy group
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
  void printFictEke(void *param, void *msg){
  
  CkReductionMsg *m=(CkReductionMsg *)msg;
  double d0   = ((double *)m->getData())[0];
  double d1   = ((double *)m->getData())[1];
  double d2   = ((double *)m->getData())[2];
  delete m;

  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
    CkPrintf("Fict Temp   =  %.10g K\n", d0); // constant=kT
    CkPrintf("Fict NhcTemp=  %.10g K\n", d1);
    CkPrintf("Fict PotNhc =  %.10g K\n", d2);
    CkPrintf("Fict EConv  =  %.10g K\n", d0+d1+d2);
  }//endif
  gSpacePlaneProxy(0,0).computeEnergies(ENERGY_FICTEKE, d0);  

}
//============================================================================


//============================================================================
//  When ALL the cp forces are done, you can integrate the atoms.
//  Could be made slicker by using knowledge of which gpsace chares
//  are on each processor and coordinating that with atom launch.
//  cklocalbranch->count_psi_forc++; zeroed in atms grp constructor
//  and incremented could be incremented in launch atoms.
//  If(cklocalbranch->count_psi_forc==numgspacethisproc){
//    atomsgrpproxy(myid).startrealspaceforces();
//  }
//  That is, you only need to know that all the forces are complete
//  on each PROCECESSOR and you can go ahead to the atoms. The atoms
//  are a group which needs to contribute atom forces to a reduction.
//  If all the forces aren't in there, you can't start the reduction
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void allDoneCPForces(void *param, void *msg){
  CkReductionMsg *m=(CkReductionMsg *)msg;
  delete m;
  CkPrintf("All done CP forces\n");
  atomsGrpProxy.StartRealspaceForces();
}
//============================================================================

//============================================================================
// When the simulation is done, make a clean exit
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
// entry method to resume execution : debugging
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::resumePsiV (CkReductionMsg *msg) {
  delete msg;
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
// Entry method to resume execution after computing reduction over all planes
// and states to form psiCgOvlap (cg only) and magforPsi
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::psiCgOvlap(CkReductionMsg *msg){
//============================================================================
// Unpack

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt    = sim->cp_min_opt;
  double tol_cp_min = sim->tol_cp_min;
  double tol_cp_dyn = sim->tol_cp_dyn;
  double d0         = ((double *)msg->getData())[0];
  double d1         = ((double *)msg->getData())[1];
         d1         = sqrt(d1); // piny convention
  delete msg;  

  fovlap_old        = fovlap; // CG ovlap (all chares need it)
  fovlap            = d0;  
  fmagPsi_total_old = fmagPsi_total;  
  fmagPsi_total     = d1;    // mag of psi force (all chares need it)

//============================================================================
// Output the mag force, send to the energy group, set the exit flag 

  if(thisIndex.x==0 && thisIndex.y==0){
    CkPrintf("MagForPsi   =  %5.8lf\n", d1);
    CkPrintf("Memory      =  %d\n",CmiMemoryUsage());
    computeEnergies(ENERGY_FMAG, d1);
  }//endif

  if(cp_min_opt==0 && fmagPsi_total>tol_cp_dyn){
    exitFlag=1;
    if(thisIndex.x==0 && thisIndex.y==0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Mag psi force too large for CP dynamics. Caio! \n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    }//endif
  }//endif

#ifndef _CP_DEBUG_ORTHO_OFF_
  if(cp_min_opt==1 && fmagPsi_total<=tol_cp_min){
    exitFlag=1;
    if(thisIndex.x==0 && thisIndex.y==0){
      CkPrintf("----------------------------------------------\n");
      CkPrintf("   CP wavefunction force tolerence reached!   \n");
      CkPrintf("----------------------------------------------\n");
    }//endif
  }//endif
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
       CkPrintf("CpuTime(GSP)= %g\n",cpuTimeNow-cpuTimeOld);
       if(cp_min_opt==0){
         CkPrintf("Psi-Rotation= %d\n",nrotation);
       }//endif
     }//endif
  }//endif

  RTH_Runtime_resume(run_thread);

//============================================================================
  }// end routine
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt  = sim->cp_min_opt;
//============================================================================

  initialized     = false;
  first_step      = 1;
  iteration       = 0;
  nrotation       = 0;
  myatom_integrate_flag=0;
  myenergy_reduc_flag=0;
  isuspend_energy=0;
  isuspend_atms=0;

  total_energy   = 0.0;
  ehart_total    = 0.0;
  enl_total      = 0.0;
  eke_total      = 0.0;
  egga_total     = 0.0;
  eexc_total     = 0.0;
  eext_total     = 0.0;
  ewd_total      = 0.0;
  fovlap         = 0.0;
  fovlap_old     = 0.0;
  fmagPsi_total  = 0.0;
  fmagPsi_total0 = 0.0; // only chare(0,0) cares
  fmagPsi_total_old = 0.0;
  cpuTimeNow     = 0.0;
  fictEke_total  = 0.0;

  iwrite_now          = 0;
  ireset_cg           = 1;
  numReset_cg         = 0;
  exitFlag            = 0;
  finishedRedPsi      = 1;  
  finishedCpIntegrate = 0;
  if(cp_min_opt==0){finishedCpIntegrate = 1;}// alternate entry point
  doneDoingIFFT       = false;
  allgdoneifft        = false;
  acceptedPsi         = true;    // we start out with a psi
  acceptedVPsi        = true;    // we start out with a vpsi
  acceptedLambda      = false;   // no forces yet
  needPsiV            = false;   // don't need tolerance check in first step
  doneNewIter         = false;

  countFileOut    = 0;
  countRedPsi     = 0;
  countPsi        = 0;
  countVPsi       = 0;
  countIFFT       = 0;
  ecount          = 0; //No energies have been received.

  numRecvRedPsi   = sim->RCommPkg[thisIndex.x].num_recv_tot;
  int ourgrain    = thisIndex.x/config.sGrainSize*config.sGrainSize; 
  if(nstates == config.sGrainSize)
    AllPsiExpected=1;
  else if(ourgrain<(nstates-config.sGrainSize)) // corner has no extras
    AllPsiExpected=2;
  else
    AllPsiExpected=1;

//============================================================================

  gSpaceNumPoints = 0;
  tpsi           = NULL;
  tvpsi          = NULL;
  tk_x           = NULL;
  tk_y           = NULL;
  tk_z           = NULL;
  ffttempdata    = NULL;
  ffttempdataGrp = NULL;
  k_x            = NULL;
  k_y            = NULL;
  k_z            = NULL;
  coef_mass      = NULL;
  initGStateSlab(&gs,sizeX,size,gSpaceUnits,realSpaceUnits,s_grain,
                 thisIndex.y,thisIndex.x);

//============================================================================

  usesAtSync = CmiTrue;
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
void CP_State_GSpacePlane::pup(PUP::er &p) {
//============================================================================
  ArrayElement2D::pup(p);

  p|first_step; //control flags and functions reference by thread are public
  p|initialized;
  p|istart_typ_cp;
  p|iteration;
  p|nrotation;
  p|exitFlag;
  p|isuspend_energy;
  p|isuspend_atms;
  p|myatom_integrate_flag;
  p|myenergy_reduc_flag;
  p|finishedCpIntegrate;
  p|finishedRedPsi;
  p|iwrite_now;
  p|ireset_cg;
  p|numReset_cg;
  p|countPsi;
  p|countVPsi;
  p|countIFFT;
  p|countFileOut;
  p|ecount;
  p|countRedPsi;
  p|numRecvRedPsi;
  p|AllPsiExpected;
  p|needPsiV;
  p|doneDoingIFFT;
  p|doneNewIter;
  p|allgdoneifft;
  p|acceptedPsi;
  p|acceptedVPsi;
  p|acceptedLambda;
  p|itemp; // 2 temporary variables for debugging in scope
  p|jtemp;

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

  p|real_proxy;   // should I pup more proxies (lambda,psi,other?)
  p|sfCompSectionProxy;
  p|gpairCalcID1;
  p|gpairCalcID2;

  // temp variables like tk_x,ffttempdata probably should be pupped.
  // Need to be careful because they could have zero length
  // or be in use with a size. They are used to catch communication.
  // Might need to add a current size variable for them.
  gs.pup(p);
  p|gSpaceNumPoints;
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
    
//-------------------------------------------------------
   }// end routine : pup
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int numData       = config.numData;
  int ibinary_opt   = sim->ibinary_opt;
  int istart_typ_cp = sim->istart_typ_cp;
  CkVec <RunDescriptor> *sortedRunDescriptors = sim->sortedRunDescriptors;
  int *npts_lgrp    = sim->npts_per_chareG;
  int *nline_lgrp   = sim->nlines_per_chareG;
  int *istrt_lgrp   = NULL;
  int *iend_lgrp    = NULL;

//============================================================================
// Set the file name using the config path and state number

  char fname[1024];
  int ind_state=thisIndex.x;
  //------------------------------------------------------------------
  // Get the complex data, Psi(g) and the run descriptor (z-lines in g-space)

  complex *complexPoints  = new complex[numData];
  complex *vcomplexPoints = NULL;
  if(istart_typ_cp>=3){vcomplexPoints = new complex[numData];}

  int *kx=  new int[numData];
  int *ky=  new int[numData];
  int *kz=  new int[numData];
  int nlines_tot,nplane,nx,ny,nz;

  if(istart_typ_cp>=3){
    sprintf(fname, "%s/vState%d.out", config.dataPath, ind_state + 1);
    readState(numData,vcomplexPoints,fname,ibinary_opt,&nlines_tot,&nplane, 
            kx,ky,kz,&nx,&ny,&nz,istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,0,1);
  }//endif

  sprintf(fname, "%s/state%d.out", config.dataPath, ind_state + 1);
  readState(numData,complexPoints,fname,ibinary_opt,&nlines_tot,&nplane, 
            kx,ky,kz,&nx,&ny,&nz,istrt_lgrp,iend_lgrp,npts_lgrp,nline_lgrp,0,0);

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
      complex *temp          = complexPoints+ioff;
      CmiMemcpy(dataToBeSent,temp,(sizeof(complex) * numPoints));

      int numPointsV;
      complex *vdataToBeSent;
      if(istart_typ_cp>=3){
         numPointsV          = numPoints;
         vdataToBeSent       = new complex[numPoints];
         complex *vtemp      = vcomplexPoints+ioff;
         CmiMemcpy(vdataToBeSent,vtemp,(sizeof(complex) * numPoints));
      }else{
         numPointsV          = 1;
         vdataToBeSent       = new complex[numPointsV];
         vdataToBeSent[0].re = 0.0;
         vdataToBeSent[0].im = 0.0;
      }//endif

      RunDescriptor *runDesc = new RunDescriptor[runsToBeSent];
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
                                                numPoints,dataToBeSent,
                                                numPointsV,vdataToBeSent,
                                                nx,ny,nz,istart_typ_cp);
      delete [] dataToBeSent;
      delete [] vdataToBeSent;
      delete [] runDesc;

      ioff += numPoints;
  }//endfor : loop over all possible chares in g-space (pencils)

  CkAssert(numData==ioff);

//============================================================================
// Clean up

  delete [] complexPoints;
  if(istart_typ_cp>=3){delete [] vcomplexPoints;}
  delete [] kx;
  delete [] ky;
  delete [] kz;

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
void CP_State_GSpacePlane::initGSpace(int            runDescSize, 
                                      RunDescriptor* runs, 
                                      int            size, 
                                      complex*       points,
                                      int            vsize, 
                                      complex*       vpoints,
                                      int nx, int ny, int nz, int istart_cp) 
//============================================================================
   { //begin routine
//============================================================================
//#ifdef _CP_DEBUG_STATEG_VERBOSE_
//  CkPrintf("initGSpace %d.%d %d\n",thisIndex.x,thisIndex.y,size);
//#endif
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt  = sim->cp_min_opt;

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

  if(istart_typ_cp>=3 && cp_min_opt==0){
    CkAssert(vsize == size);
    CmiMemcpy(gs.packedVelData, vpoints, sizeof(complex)*gs.numPoints);
  }else{
    memset(gs.packedVelData, 0, sizeof(complex)*gs.numPoints);
  }//endif

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
//  gs.eke_ret = 0;

//============================================================================
// Init NHC, Sample velocities 

  int ncoef_use = gs.numPoints-gs.nkx0_red;
  int maxLenNHC = LEN_NHC_CP;
  int maxNumNHC = NUM_NHC_CP;
  CPINTEGRATE::initCPNHC(ncoef,maxLenNHC,maxNumNHC,&gs.len_nhc_cp,&gs.num_nhc_cp,
                         &gs.kTCP,&gs.tauNHCCP,&gs.mNHC,&gs.degfree,&gs.degfreeNHC,
                         &gs.gammaNHC,ncoef_use,gs.nkx0_zero);


  if(cp_min_opt==0){
   gs.xNHC=0.0;
   CPINTEGRATE::CPSmplVel(gs.numPoints,coef_mass,gs.packedVelData,
                         gs.len_nhc_cp,gs.num_nhc_cp,gs.mNHC,gs.vNHC,gs.kTCP,
                         istart_typ_cp,gs.nkx0_red,gs.nkx0_uni,gs.nkx0_zero,
                         gs.degfree,gs.degfreeNHC,gs.gammaNHC);
  }//endif

#ifdef _CP_DEBUG_PSI_OFF_
  memset(gs.packedPlaneData, 0, sizeof(complex)*gs.numPoints);
#endif
//#ifdef _CP_DEBUG_DYNAMICS_
  memset(gs.packedVelData, 0, sizeof(complex)*gs.numPoints);
//#endif

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
// This reduction is done to signal the end of initialization to main

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
void CP_State_GSpacePlane::makePCproxies(){
  lambdaproxy=makeOneResultSection_asym(&gpairCalcID2, thisIndex.x, thisIndex.y);
  psiproxy=makeOneResultSection_sym1(&gpairCalcID1, thisIndex.x, thisIndex.y);
  if(AllPsiExpected>1)
    psiproxyother=makeOneResultSection_sym2(&gpairCalcID1, thisIndex.x, thisIndex.y);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::startNewIter ()  {
//============================================================================
// Check for flow of control errors :

  if(iteration>0){
   if(egroupProxy.ckLocalBranch()->iteration_gsp != iteration || 
     atomsGrpProxy.ckLocalBranch()->iteration  != iteration){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Starting new iter before\n");
      CkPrintf("finishing atom integrate or iteration mismatch.\n");
      CkPrintf("iter_gsp %d iter_energy %d iter_atm %d\n",
 	        iteration,
	        egroupProxy.ckLocalBranch()->iteration_gsp,
                atomsGrpProxy.ckLocalBranch()->iteration);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
  } //endif

  doneNewIter = true;
  myenergy_reduc_flag  =0; // energies not reduced for new iter
  myatom_integrate_flag=0; // atoms not integrated on new iter

//============================================================================
// Reset all the counters that need to be reset (not more not less)
// otherwise race conditions can leak in.  Rely on the constructor
// for initialization.  Reset set your flags as soon as you are done
// with the tests that require them.

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int cp_min_opt = sim->cp_min_opt;

  // Finished integrate and red psi are safe.
  // You can't get to these until you get passed through this routine
  finishedCpIntegrate = 0;
  finishedRedPsi      = 1;   if(cp_min_opt==0){finishedRedPsi=0;}

  iteration++;   // my iteration # : not exactly in sync with other chares
                 //                  but should agree when chares meet.

//============================================================================
// Check Load Balancing, Increment counter, set done flags equal to false.

#ifndef CMK_OPTIMIZE
  if(iteration==TRACE_ON_STEP ){traceBegin();}
  if(iteration==TRACE_OFF_STEP){traceEnd();}
#endif

    if(config.lbgspace || config.lbpaircalc){
	if((iteration % (FIRST_BALANCE_STEP - PRE_BALANCE_STEP) == 0)  || 
           (iteration % (LOAD_BALANCE_STEP - PRE_BALANCE_STEP) == 0)){
	    LBTurnInstrumentOn();
	}//endif
    }//endif

//============================================================================
// Output psi at start of minimization for debugging

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
        msg->iteration_src = iteration;
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

#ifndef CMK_OPTIMIZE    
  double StartTime=CmiWallTimer();
#endif

// Do fft in forward direction, 1-D, in z direction
// A local function not a message : get pointer to memory for fft group

  ffttempdataGrp = fftCacheProxy.ckLocalBranch()->doGSRealFwFFT(gs.packedPlaneData, 
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
    for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = ffttempdataGrp[j];}
    real_proxy(thisIndex.x, z).doFFT(msg);  // same state, realspace char[z]

  }//endfor
    
  if (config.useCommlib){mssInstance.endIteration();}
    
  ffttempdataGrp = NULL; // its memory from the FFTgroup : don't touch it

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
  countIFFT++;
  if (countIFFT == gs.planeSize[1]) {
    countIFFT = 0;
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
  fftw_free(ffttempdata); // its mine, I can free it if I like
  ffttempdata = NULL;
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(GspaceBwFFT_, StartTime, CmiWallTimer());
#endif    

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

  doneDoingIFFT = true;
//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
bool CP_State_GSpacePlane::doneNLForces(){
//==============================================================================
return particlePlaneProxy(thisIndex.x,thisIndex.y).ckLocal()->doneGettingForces;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNLForces(){
//================================================================================
//  I. You have called resume until these guys report. You have not started FFT
#ifdef _CP_DEBUG_NONLOC_BARRIER_
  RTH_Runtime_resume(run_thread);
#endif
//================================================================================
//  I. If the fft forces finished first, you need to resume(they are waiting).
// II. If the fft forces are not done, hang loose till they finish.
#ifndef _CP_DEBUG_NONLOC_BARRIER_
  if(doneDoingIFFT){RTH_Runtime_resume(run_thread);}
#endif
//-------------------------------------------------------------------------------
  }//end routine   
//================================================================================

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
  myenergy_reduc_flag=0;

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
#ifdef _CP_DEBUG_PSI_OFF_
  iteration++;
  if(iteration==config.maxIter+1){
    int i=0;
    contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
  }else{
#endif
   int i=0;
   contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(allDoneCPForces,NULL));
#ifdef _CP_DEBUG_PSI_OFF_
  }//endif
#endif

}//end routine
//===============================================================================



//==============================================================================
// After MY Cp forces have arrived : sendLambda
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendLambda() {
//==============================================================================
// Reset set lambda (not done) and force counters (not done for NEXT step):

   CP_State_ParticlePlane *pp=particlePlaneProxy(thisIndex.x,thisIndex.y).ckLocal();
   acceptedLambda        = false;
   pp->doneGettingForces = false;
   doneDoingIFFT         = false;
   doneNewIter           = false;

//==============================================================================
// Scale the variables and launch lambda

  complex *psi   = gs.packedPlaneData;
  complex *force = gs.packedForceData;
  int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;

#ifndef _CP_DEBUG_ORTHO_OFF_
  if(gs.ihave_kx0==1 && cp_min_opt==0){
    double rad2i = 1.0/sqrt(2.0);
    double rad2  = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i]   *= rad2i;}
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){force[i] *= rad2;}
  }//endif
#endif

  if(config.gSpaceNumChunks!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, while gSpaceNumChunk!=1 is cool, I'm \n");
    CkPrintf("afraid you'll have to do the implementation yourself!\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  int numPoints   = gs.numPoints / config.gSpaceNumChunks;
  int dataCovered = gs.numPoints;
#ifndef _CP_DEBUG_ORTHO_OFF_
  int c = 0;
  int toSend = (c == config.gSpaceNumChunks - 1) ? dataCovered : numPoints;
  startPairCalcLeft(&gpairCalcID2, toSend, psi + c * numPoints, 
		    thisIndex.x, thisIndex.y, false);
  startPairCalcRight(&gpairCalcID2, toSend, force + c * numPoints, 
		     thisIndex.x, thisIndex.y);
#else
  acceptedLambda=true;
  memset(force, 0, sizeof(complex)*numPoints);
#endif

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
// Retrieve Non-orthog psi

  if(cp_min_opt==0){
     int ncoef          = gs.numPoints;
     complex *psi_g     = gs.packedPlaneData;
     complex *psi_g_scr = gs.packedPlaneDataScr;
     memcpy(psi_g,psi_g_scr,sizeof(complex)*ncoef);
  }//endif

//==============================================================================
// Resume 

  acceptedLambda=true;
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
// Flow of control check and local variables

  if(!acceptedLambda){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Flow of Control Error : Attempting to Cg ovlap\n");
     CkPrintf("without lambda correction\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif

   int istate      = gs.istate_ind;
   int ncoef       = gs.numPoints;
   complex *psi_g  = gs.packedPlaneData;
   complex *forces = gs.packedForceData;
   int cp_min_opt  = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
   int cp_min_cg   = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_cg;

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

//----------------------------------------------------------------------------
  }// end routine : computeCgOverlap
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// In this function data is written to files the simpliest way possible
//============================================================================
void CP_State_GSpacePlane::writeStateDumpFile() {
//============================================================================
// Local pointers, variables and error checking

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int cp_min_opt = (sim->cp_min_opt);
  int ndump_frq  = (sim->ndump_frq);

  int ind_state  = (thisIndex.x+1);
  int ind_chare  = (thisIndex.y+1);
  int ncoef      = gSpaceNumPoints;

  complex *psi      = gs.packedPlaneData;
  complex *vpsi     = gs.packedPlaneDataScr;  // switch definitions
  complex *vpsi_old = gs.packedVelData;       // to evolve psi to time, t.
  complex *forces   = gs.packedForceData;
  double **vNHC     = gs.vNHC_scr;            // switch to evolve to time,t.
  double **fNHC     = gs.fNHC;          
  double mNHC       = gs.mNHC;          
  int len_nhc_cp    = gs.len_nhc_cp;
  int num_nhc_cp    = gs.num_nhc_cp;
  int nkx0_red      = gs.nkx0_red;
  int nkx0_uni      = gs.nkx0_uni;
  int nkx0_zero     = gs.nkx0_zero;
  double kTCP       = gs.kTCP;
  double xNHC       = 0.0;

  if(!acceptedPsi || !acceptedLambda || !acceptedVPsi){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Flow of Control Error : Attempting to write states\n");
     CkPrintf("without completing psi, vpsi and Lambda\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// Set the file names and write the files

  iwrite_now = 0;
  if(config.stateOutputOn==1){

    if( ((iteration % ndump_frq)==0) || (iteration==config.maxIter) || 
          exitFlag==1){
    //------------------------------------------------------------------
    // Set the flag and tell the world you are writing
      iwrite_now = 1;
      if(ind_state==1 && ind_chare==1){
        CkPrintf("-----------------------------------\n");
        CkPrintf("Writing states to disk on step %d\n",iteration);
        CkPrintf("-----------------------------------\n");
      }//endif
    //------------------------------------------------------------------
    // Update the velocities into scratch as we are between steps
      if(cp_min_opt==0){
        memcpy(vpsi,vpsi_old,sizeof(complex)*ncoef);
        gs.copyVNHC();
        if(iteration>1){
          CPINTEGRATE::cp_evolve_vel(ncoef,forces,vpsi,coef_mass,
                       len_nhc_cp,num_nhc_cp,fNHC,vNHC,&xNHC,mNHC,kTCP,
                       nkx0_red,nkx0_uni,nkx0_zero,2,iteration,
                       gs.degfree,gs.degfreeNHC,gs.gammaNHC);
	}//endif
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
      for (int i=0;i<ncoef; i++){
        data[i]  = psi[i];  vdata[i] = vpsi[i];
        mk_x[i]  = k_x[i];  mk_y[i]  = k_y[i];  mk_z[i]  = k_z[i];
      }//endfor
      gSpacePlaneProxy(thisIndex.x,0).collectFileOutput(msg);
    //------------------------------------------------------------------
    // If you are not plane 0, you are done. Invoke the correct reduction.
      if(thisIndex.y!=0){
        int i = 0;
        if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==1){
          contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
	}else{
          contribute(sizeof(int),&i,CkReduction::sum_int,
                   CkCallback(CkIndex_CP_State_GSpacePlane::psiWriteComplete(NULL),
                              gSpacePlaneProxy));
	}//endif
      }//endif
    }//endif : its time to write

  }//endif : it is useful to write
  
//---------------------------------------------------------------------------
   }//write the file
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::collectFileOutput(GStateOutMsg *msg){
//============================================================================

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int sizeX        = (sim->sizeX);
  int sizeY        = (sim->sizeY);
  int sizeZ        = (sim->sizeZ);
  int npts_tot     = (sim->npts_tot);
  int *ipacked_off = (sim->index_output_off);
  int cp_min_opt   = (sim->cp_min_opt);
  int ibinary_write_opt = sim->ibinary_write_opt;

//============================================================================
// Receive the message

  if(countFileOut==0){
    tpsi  = new complex[npts_tot];
    tvpsi = new complex[npts_tot];
    tk_x  = new int[npts_tot];
    tk_y  = new int[npts_tot];
    tk_z  = new int[npts_tot];
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

//============================================================================
// If you've got the whole state, write it out and then invoke the reduction.

  if(countFileOut==nchareG){
     countFileOut = 0;
     int ind_state = thisIndex.x+1;
     char psiName[200]; char vpsiName[200];
       sprintf(psiName, "%s/newState%d.out", config.dataPath,ind_state);
       sprintf(vpsiName,"%s/newVstate%d.out",config.dataPath,ind_state);
       writeStateFile(npts_tot,tpsi,tvpsi,tk_x,tk_y,tk_z,cp_min_opt,
                      sizeX,sizeY,sizeZ,psiName,vpsiName,ibinary_write_opt);
     delete [] tpsi; tpsi  = NULL;
     delete [] tvpsi;tvpsi = NULL;
     delete [] tk_x; tk_x  = NULL;
     delete [] tk_y; tk_y  = NULL;
     delete [] tk_z; tk_z  = NULL;
     int i=0;
     if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==1){
       contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
     }else{
       contribute(sizeof(int),&i,CkReduction::sum_int,
                CkCallback(CkIndex_CP_State_GSpacePlane::psiWriteComplete(NULL),
                           gSpacePlaneProxy));
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

//==============================================================================
// II. Local pointers    

  int cp_min_opt     = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  int cp_min_cg      = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_cg;

  int istate         = gs.istate_ind;
  int ncoef          = gs.numPoints;
  int len_nhc        = gs.len_nhc_cp; 
  int num_nhc        = gs.num_nhc_cp;
  complex *psi_g     = gs.packedPlaneData; 
  complex *forces    = gs.packedForceData; 
  complex *vpsi_g    = gs.packedVelData; // for cp not minimization
  complex *forcesold = gs.packedVelData; // for miniziation not cp
  double **fNHC      = gs.fNHC;
  double **vNHC      = gs.vNHC;
  double *xNHC       = &gs.xNHC;
  double mNHC        = gs.mNHC;
  double kTCP        = gs.kTCP;
  int nkx0_red       = gs.nkx0_red;
  int nkx0_uni       = gs.nkx0_uni;
  int nkx0_zero      = gs.nkx0_zero;
  double degfree     = gs.degfree;
  double degfreeNHC  = gs.degfreeNHC;
  double gammaNHC    = gs.gammaNHC;
  double fictEke=0.0;
  double ekeNhc=0.0;

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

#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif

  fictEke = 0.0; ekeNhc=0.0;
  CPINTEGRATE::CP_integrate(ncoef,istate,iteration,forces,forcesold,psi_g,
               coef_mass,k_x,k_y,k_z,len_nhc,num_nhc,fNHC,vNHC,xNHC,mNHC,kTCP,
               gamma_conj_grad,&fictEke,nkx0_red,nkx0_uni,nkx0_zero,
               &ekeNhc,degfree,degfreeNHC,gammaNHC);

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(IntegrateModForces_, StartTime, CmiWallTimer());
#endif

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
// V. Contribute FictKe : Output and store in energy group

  double redSize = ((double) (nchareG*nstates));
  double sendme[3];
  sendme[0] = 2.0*BOLTZ*(fictEke/degfree)/redSize;
  sendme[1] = 2.0*BOLTZ*(ekeNhc/degfreeNHC)/redSize;
  sendme[2] = 2.0*BOLTZ*(xNHC[0]/degfreeNHC)/redSize;
//  sendme[0] = fictEke;
  contribute(3*sizeof(double),sendme,CkReduction::sum_double, 
             CkCallback(printFictEke, NULL));

//==========================================================================
// VI. Pack up and set the flag that indicating you've finished integrating.

  gs.fictEke_ret      = fictEke;
  gs.ekeNhc_ret       = ekeNhc;
  finishedCpIntegrate = 1;
  memset(forces, 0, sizeof(complex)*ncoef);

//------------------------------------------------------------------------------
   }// end CP_State_GSpacePlane::integrateModForce
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::sendRedPsi() {
//==============================================================================

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  complex *sendData = gs.packedPlaneData;
  int isend         = thisIndex.y;       // my g-space chare index
  int  *num_send    = RCommPkg[isend].num_send;
  int **lst_send    = RCommPkg[isend].lst_send;
  int num_recv_tot  = RCommPkg[isend].num_recv_tot;
  int num_send_tot  = RCommPkg[isend].num_send_tot;
  numRecvRedPsi     = num_recv_tot;

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
         gSpacePlaneProxy(thisIndex.x,irecv).acceptRedPsi(msg);
      }//endif

    }//endfor

  }//endif : no one to which I have to send

//==============================================================================
// Check for errors 

  if(iii!=num_send_tot){
    CkPrintf("Error in GSchare %d %d : %d %d\n",thisIndex.x,thisIndex.y,
 	                                           num_send_tot,iii);
    CkExit();
  }//endif
  if(numRecvRedPsi==0 && gs.nkx0_red>0){
    CkPrintf("Error in GSchare %d %d : %d %d\n",thisIndex.x,thisIndex.y,
	                                        numRecvRedPsi,gs.nkx0_red);
    CkExit();
  }//endif
  if(jjj != gs.nkx0_uni-gs.nkx0_zero){
    CkPrintf("Error in GSchare %d %d : %d %d\n",thisIndex.x,thisIndex.y,
	                                        jjj,gs.nkx0_uni);
    CkExit();
  }//endif

//==============================================================================
// If you don't have any receiving to do, you are done.

 if(num_recv_tot==0){finishedRedPsi=1;}

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptRedPsi(GSRedPsiMsg *msg) {
//==============================================================================

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  int ncoef         = msg->size;
  int isend         = msg->senderIndex;
  complex *msgData  = msg->data;
  complex *recvData = gs.packedRedPsi;
  int irecv         = thisIndex.y;       // my g-space chare index
  int  *num_recv    = RCommPkg[irecv].num_recv;
  int **lst_recv    = RCommPkg[irecv].lst_recv;
  int num_recv_tot  = RCommPkg[irecv].num_recv_tot;

//==============================================================================
// unpack

  if(num_recv[isend]!=ncoef){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Number sent not equal to number reciever expected \n");
    CkPrintf("Sender %d size %d Reciever %d size %d \n",isend,ncoef,irecv,num_recv[isend]);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(countRedPsi==0){itemp=0;jtemp=0;}
  jtemp+= ncoef;
  if(ncoef>0){itemp++;}

  for(int i=0;i<ncoef;i++){
    recvData[lst_recv[isend][i]].re= msgData[i].re;
    recvData[lst_recv[isend][i]].im=-msgData[i].im;
  }//endfor

  delete msg;

//==============================================================================
// Done

  countRedPsi++;
  if(countRedPsi==num_recv_tot){
    countRedPsi=0;
    finishedRedPsi++;
    if(itemp!=num_recv_tot){
      CkPrintf("Error in GSchare recv cnt %d %d : %d %d\n",thisIndex.x,thisIndex.y,
 	                                                   num_recv_tot,itemp);
      CkExit();
    }//endif
    if(jtemp!=gs.nkx0_red){
      CkPrintf("Error in GSchare recv cnt %d %d : %d %d\n",thisIndex.x,thisIndex.y,
 	                                                   gs.nkx0_red,jtemp);
      CkExit();
    }//endif
    if(finishedCpIntegrate==1){       // Only suspend if integrate done and Red not
      RTH_Runtime_resume(run_thread); // Only resume if Red just finished & integrate done
    }//endif
  }//endif

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doneRedPsiIntegrate() {

    finishedRedPsi++; //extra points for getting here
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt = sim->cp_min_opt;

  if(config.gSpaceNumChunks!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, while gSpaceNumChunk!=1 is cool, I'm \n");
    CkPrintf("afraid you'll have to do the implementation yourself!\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(finishedRedPsi==0 || finishedCpIntegrate==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't sendPsi without completing integrate\n");
    if(cp_min_opt==0){
      CkPrintf("and sending the Redundant psi values around\n");
    }//endif
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(cp_min_opt==0 && iteration>0){
    int ncoef  = gs.nkx0_red;
    if(ncoef>0 && finishedRedPsi!=2){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you can't sendPsi without receiving Redpsi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//==============================================================================
// Prepare the data : If cp dynamics is going, save the non-orthogonal puppies.

  acceptedPsi =false;

  complex *data=gs.packedPlaneData;

  /*  if(tpsi==NULL)
    tpsi  = new complex[gs.numPoints];
  memcpy(tpsi,data,sizeof(complex)*gs.numPoints);
  */
  if(cp_min_opt==0){
    //#ifdef  _CP_DEBUG_UPDATE_OFF_
    //    CmiMemcpy(gs.packedPlaneData,gs.packedPlaneDataTemp,
    //                sizeof(complex)*gs.numPoints);
    //    memset(gs.packedVelData,0,sizeof(complex)*gs.numPoints);
    //#endif      
     int ncoef     = gs.numPoints;
     complex *scr  = gs.packedPlaneDataScr; //save non-orthog psi
     memcpy(scr,data,sizeof(complex)*ncoef);
  }//endif

#ifndef _CP_DEBUG_ORTHO_OFF_
  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){data[i] *= rad2i;}
  }//endif
#endif

//==============================================================================

  int s, c;
  int idx = (thisIndex.x/gs.S_grainSize) * gs.S_grainSize;
  int numPoints = gs.numPoints / config.gSpaceNumChunks;
  int dataCovered = gs.numPoints;

#ifndef _CP_DEBUG_ORTHO_OFF_
  c = 0;
  int toSend = (c == config.gSpaceNumChunks - 1) ? dataCovered : numPoints;
  startPairCalcLeft(&gpairCalcID1, toSend, data + c * numPoints, 
		    thisIndex.x, thisIndex.y, false);
#else
  acceptedPsi=true;
  if((iteration==config.maxIter || exitFlag==1) && cp_min_opt==1 && 
      config.stateOutputOn==0){
      int i;
      contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
  }///endif
#endif

//----------------------------------------------------------------------------
    }// end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void  CP_State_GSpacePlane::sendPsiV() {
//==============================================================================

  nrotation++;
  acceptedVPsi=false;
  if(config.gSpaceNumChunks!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, while gSpaceNumChunk!=1 is cool, I'm \n");
    CkPrintf("afraid you'll have to do the implementation yourself!\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  complex *data=gs.packedVelData;
  //  complex *data=tpsi;

  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
     int ncoef     = gs.numPoints;
     complex *scr  = gs.packedPlaneDataScr;  // replace no-ortho psi 
     complex *psi  = gs.packedPlaneData;     // by orthonormal psi
     memcpy(scr,psi,sizeof(complex)*ncoef);
  }//endif

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

  if(countPsi<1){
    for(int i=0; i<N; i++){psi[i] = data[i];}
  }else{
    for(int i=0; i<N; i++){psi[i] += data[i];}
  }//endif
  delete msg;

//=============================================================================
// (II) If you have got it all : Rescale it, produce some output

  countPsi++;//psi arrives in as many as 2 reductions
  if(countPsi==AllPsiExpected){ 
    //--------------------------------------------------------------------
    // (A) Reset counters and rescale the kx=0 stuff
    acceptedPsi=true;
    countPsi=0;
    if(gs.ihave_kx0==1){
      double rad2 = sqrt(2.0);
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i] *= rad2;}
    }//endif
    //--------------------------------------------------------------------
    // (B) Generate some screen output of orthogonal psi
    screenOutputPsi();
    //--------------------------------------------------------------------
    // (D) Go back to the top or exit
    int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
    if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==0){
      int i;
      contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
    }//endif
    if((iteration==config.maxIter || exitFlag==1) && cp_min_opt==1 && 
        config.stateOutputOn==0){
      int i;
      contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
    }///endif
    //--------------------------------------------------------------------
    // (E) Reset psi if debugging
#ifdef  _CP_DEBUG_UPDATE_OFF_
     if(cp_min_opt==1){
       CmiMemcpy(gs.packedPlaneData,gs.packedPlaneDataTemp,
               sizeof(complex)*gs.numPoints);
       memset(gs.packedVelData,0,sizeof(complex)*gs.numPoints);
     }//endif
#endif      
    RTH_Runtime_resume(run_thread);
  } // if countPsi

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

  if(countPsi<1){
    for(int i=0; i<N; i++){vpsi[i] = data[i];}
  }else{
    for(int i=0; i<N; i++){vpsi[i] += data[i];}
  }//endif
  
  delete msg;

//=============================================================================
// (II) If you have got it all : Rescale it, produce some output

  countPsi++;//psi arrives in as many as 2 reductions

  if(countPsi==AllPsiExpected){ 
    //--------------------------------------------------------------------
    // (A) Reset counters and rescale the kx=0 stuff
    acceptedVPsi=true;
    countPsi=0;
    if(gs.ihave_kx0==1){
      double rad2 = sqrt(2.0);
      for(int i=gs.kx0_strt; i<gs.kx0_end; i++){vpsi[i] *= rad2;}
    }//endif
    /* debugging barrier
    int wehaveours=1;
    contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
      CkCallback(CkIndex_CP_State_GSpacePlane::gdonePsiV(NULL),gSpacePlaneProxy));
    */
    if(!acceptedPsi){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Attempting to v-psi\n");
      CkPrintf("resume without completing psi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    needPsiV=false;
    RTH_Runtime_resume(run_thread);
  } // if countPsi

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
// Probe for atom completion : Could have atom group invoke this function directly
//                             on all chares for which it is responsible.
//                             There is an atom group on each proc. Some
//                             number of gspaceplane objects reside on each proc.
//                             This information is known by main and could
//                             be put into atom constructor. The atom group
//                             would then know which chares upon which to 
//                             invoke acceptAtoms. Careful, careful with migration
//                             with this alternative scheme.
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptAtoms(GSAtmMsg *msg) {
   delete msg;
   myatom_integrate_flag=1;
   if(atomsGrpProxy.ckLocalBranch()->iteration != iteration){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Iteration\n");
      CkPrintf("mismatch between atoms and g-space planes\n");
      CkPrintf("suspend %d iteration_atm %d iteration_gsp %d\n",
                isuspend_atms,iteration,atomsGrpProxy.ckLocalBranch()->iteration);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
   if(isuspend_atms==1){ // I suspended to wait for atoms, resume me baby.
     isuspend_atms=0;
     RTH_Runtime_resume(run_thread);
   }//endif
}//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptEnergy(GSAtmMsg *msg) {
   delete msg;
   myenergy_reduc_flag=1;
   if(egroupProxy.ckLocalBranch()->iteration_gsp != iteration){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Iteration\n");
      CkPrintf("mismatch between energy and g-space planes\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
   if(isuspend_energy==1){ // I suspended to wait for energy, resume me baby.
     isuspend_energy=0;
     RTH_Runtime_resume(run_thread);
   }//endif
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

  if(AllPsiExpected>1){
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
      ecount++;
      break;
      
    case ENERGY_FMAG : 
      fmagPsi_total0 = d;
      ecount++;
      break;
      
    default :
      CkAbort("unknown energy");
      break;
  }//end switch

//==============================================================================
// if you debuggin you get fewer energies

  int isub =0;
#ifdef _CP_DEBUG_SFNL_OFF_
  isub++;
#endif
#ifdef  _CP_DEBUG_RHO_OFF_
  isub+=5;
#endif
#ifdef _CP_DEBUG_HARTEEXT_OFF_
#ifndef  _CP_DEBUG_RHO_OFF_
  isub+=3;
#endif
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
    egroupProxy.updateEnergiesFromGS(estruct); // broadcast the electronic energies
                                               //  so that all procs have them
    total_energy        = 0.0;
    ecount              = 0;

  }// got all the energies

//-----------------------------------------------------------------------------
   }// end routine : computenergies
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::requirePsiV() {
  needPsiV=true;
  int foo=1;
  // when everyone is ready, restart Ortho's backward path
  contribute(sizeof(int), &foo, CkReduction::min_int, 
             CkCallback(CkIndex_Ortho::resumeV(NULL), orthoProxy));
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
bool CP_State_GSpacePlane::weneedPsiV() { return needPsiV;}
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
