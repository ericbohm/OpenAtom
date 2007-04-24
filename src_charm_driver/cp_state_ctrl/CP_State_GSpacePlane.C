//======================================================
// Things to do : 
//    move resetiterstate
//======================================================

//#define _CP_DEBUG_WARN_SUSPEND_
//#define _CP_DEBUG_NONLOC_BARRIER_
//#define _CP_DEBUG_ORTHO_OFF_
//#define _CP_DEBUG_PSI_OFF_
//#define GPSI_BARRIER
//#define GIFFT_BARRIER
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
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "eesCache.h"
#include "fftCacheSlab.h"
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

extern CProxy_main                    mainProxy;
extern CProxy_TimeKeeper              TimeKeeperProxy;
extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CP_State_GSpacePlane    gSpacePlaneProxy;
extern CProxy_Ortho                   orthoProxy;
extern CProxy_CP_State_ParticlePlane  particlePlaneProxy;
extern CProxy_CPcharmParaInfoGrp      scProxy;
extern CProxy_AtomsGrp                atomsGrpProxy;
extern CProxy_StructureFactor         sfCompProxy;
extern CProxy_EnergyGroup             egroupProxy; //energy group proxy
extern CProxy_FFTcache                fftCacheProxy;
extern CProxy_StructFactCache         sfCacheProxy;
extern CProxy_eesCache                eesCacheProxy;

extern CProxy_ComlibManager mgrProxy;
extern ComlibInstanceHandle gssInstance;
extern ComlibInstanceHandle mcastInstancePP;

extern int nstates;
extern int sizeX;
extern int nchareG;              // number of g-space chares <= sizeX and >=nplane_x

void cleanExit(void *, void *);
void allDoneCPForces(void *, void *);
void printFictEke(void *, void *);
void printEnergyEke(void *, void *);
void testeke(int ,complex *,int *,int *,int *, int ,int);

//#define _CP_DEBUG_STATEG_VERBOSE_
//#define _CP_DEBUG_WARN_SUSPEND_

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
       if(c->ees_nonlocal==0){
         c->releaseSFComputeZ();
#ifdef _CP_DEBUG_NONLOC_BARRIER_
         RTH_Suspend();  // resume is called by acceptNLForces()
#endif  // Barrier for non-local
       }//endif
#endif // non-local is allowed
    //------------------------------------------------------------------------
    // (C) Before starting any state related comps wait for everyone to arrive
#ifdef GPSI_BARRIER  // pause for every single chare to finish
       if(!(c->allAcceptedPsiDone())){
	  RTH_Suspend(); // wait for broadcast that all psi is done  
       }//endif
#endif               //end pause
    //------------------------------------------------------------------------
    // (D) FFT psi(gx,gy,gz)->psi(gx,gy,z), Send psi to real, 
#ifndef _CP_DEBUG_VKS_OFF_ // if vks forces are allowed do the following and
                           // EES-NL launch is in realstate
       c->doFFT(); 
       c->sendFFTData();
#else
#ifndef _CP_DEBUG_SFNL_OFF_
       if(c->ees_nonlocal==1){c->startNLEes(true);}  // EES-NL launch must be here or later
#endif
#endif
    //------------------------------------------------------------------------
    // (F) When Psi forces come back to us, do back FFT
#ifndef _CP_DEBUG_VKS_OFF_ // if vks forces are allowed
       RTH_Suspend();  // wait for (psi*vks)=F[gx,gy,z] to arive from RealSpace
       c->doIFFT();    // Message from realspace arrives : doifft(msg) resumes
#else
       c->doneDoingIFFT = true;
#endif
    //------------------------------------------------------------------------
    // (G) If NL-pseudo forces are not done, wait for them.
#ifndef _CP_DEBUG_SFNL_OFF_ // non-local is allowed
       c->isuspendNLForces = 0;
       if(!c->doneNLForces()){
         c->isuspendNLForces = 1;
         RTH_Suspend(); // resume called in acceptNLforces or acceptNLForcesEes
                        // when NL forces are ALL done (all channels/atm types).
                        // The latter invoked by CP_State_ParticlePlane.
       }//endif
#endif // non-local is allowed
    //------------------------------------------------------------------------
    // (H) Added up all force contributions
       c->combineForcesGetEke();
#ifdef GIFFT_BARRIER  // pause for every single chare to finish
       if(!(c->allDoneIFFT())){
	  RTH_Suspend(); // wait for broadcast that all gspace is done  
       }//endif
#endif //end pause
#endif // only move the atoms
    //------------------------------------------------------------------------
    // I) The atoms can't go until cp forces are completely finished (all chares)
    //    However, the atoms can overlap with all this lambda, psi stuff.
       c->launchAtoms();
    //------------------------------------------------------------------------
    // (G) Add contraint forces (rotate forces to non-orthogonal frame)
#ifndef _CP_DEBUG_PSI_OFF_  // you are moving everything
       c->sendLambda();
#ifndef _CP_DEBUG_ORTHO_OFF_
       RTH_Suspend(); // wait for forces to be fixed up 
                      // acceptLambda resumes
#endif
    //------------------------------------------------------------------------
    // (H) Get sum sq forces (even under dynamics its good to have) : also cg thingy
       c->computeCgOverlap();
       RTH_Suspend(); // wait for cg reduction : psiCgOvlap resumes
    //------------------------------------------------------------------------
    // (I) Output the states for cp dynamics  : iter==1 nothing has moved
       if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0 &&
          c->iteration<= config.maxIter && c->iteration>1){
          c->writeStateDumpFile();   // wait for output : psiwritecomplete resumes
          if(c->iwrite_now==1){RTH_Suspend();}
       }//endif
    //------------------------------------------------------------------------
    // (J) Evolve the electrons to the next step    
       c->integrateModForce();
       if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0){
          c->sendRedPsi();  // Sync Redundant psi entries
          if(c->iRecvRedPsi==0){
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
    c->sendPsi();   // send to Pair Calculator
#ifndef _CP_DEBUG_ORTHO_OFF_
    RTH_Suspend();  // Wait for new Psi : 
                    // resume is called in acceptNewPsi
#endif
   //------------------------------------------------------------------------
   // (B) Output the states for minimization
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 1 &&
       c->iteration<= config.maxIter){
       c->writeStateDumpFile();
       if(c->iwrite_now==1){RTH_Suspend();}// wait : psiwritecomplete resumes
    }//endif
    //------------------------------------------------------------------------
    // (C) Velocity Norb rotation
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt == 0 && c->iteration>1 && 
       c->iteration<= config.maxIter){
       if(c->weneedPsiV()){//ortho will have told us this
          c->sendRedPsiV();      // Sync psiV at time t : resume called in acceptRedPsiV
          if(c->iRecvRedPsiV==0){
             RTH_Suspend();       // Resume called in acceptRedPsiV
	  }//endif
          c->doneRedPsiVIntegrate(); // tuck away your g=0 plane velocities
          c->sendPsiV();         // Rotate yourself 
          RTH_Suspend();         // Wait for new PsiV : resume is called in acceptNewPsiV
       }//endif
    }//endif
   //------------------------------------------------------------------------
   // (D) Check for Energy reduction completion : should just be a safety
    if(c->myenergy_reduc_flag==0 && c->iteration>0 && c->isuspend_energy==0){
       c->isuspend_energy=1;
    }//endif
#endif // you are moving everything
   //------------------------------------------------------------------------
   // (E) Check for atom integration : should just be a safety
    if(c->myatom_integrate_flag==0 && c->iteration>0 && c->isuspend_atms==0){
      c->isuspend_atms=1;
    }//endif
   //------------------------------------------------------------------------
   // (F) If the atom or energy stuff is slow, relax for a bit
    if(c->isuspend_atms==1 || c->isuspend_energy==1){
#ifdef _CP_DEBUG_WARN_SUSPEND_
      CkPrintf("Suspend atm/energy on proc %d : chare %d %d : %d %d\n",
             CkMyPe(),c->istate_ind,c->iplane_ind,c->isuspend_atms,c->isuspend_energy);
#endif
      RTH_Suspend();     // resume called in acceptEnergy or in acceptAtoms
    }//endif
   //------------------------------------------------------------------------
   // (G) If you have triggered an exit condition just chill until ckexit
    if(c->cleanExitCalled==1){RTH_Suspend();} 
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
//     All GSpace objects have finished psi : Debugging only
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::gdonePsi(CkReductionMsg *msg){
      delete msg;
      //let my ffts go!
      allAcceptedPsi=true;
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
  double d3   = ((double *)m->getData())[3];
  double d4   = ((double *)m->getData())[4];
  delete m;

  if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0){
    CkPrintf("Fict Temp   =  %.10g K\n", d0); // per g-chare temp
    CkPrintf("Fict Eke    =  %.10g K\n", d2); // total kinetic energy
    CkPrintf("Fict TempNHC=  %.10g K\n", d1); // per g-chare tempNHC
    CkPrintf("Fict EkeNHC =  %.10g K\n", d3); // total NHC kinetic energy
    CkPrintf("Fict PotNHC =  %.10g K\n", d4); // total potNHC
    CkPrintf("Fict EConv  =  %.10g K\n", d2+d3+d4);
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
  PRINT_LINE_STAR; CkPrintf("\n"); CkPrintf("\n");

  PRINT_LINE_STAR;
  CkPrintf("         Open Atom Simulation Complete                \n");
  PRINT_LINE_STAR;
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
  AtomsGrp *ag         = atomsGrpProxy.ckLocalBranch(); // find me the local copy

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
    CkPrintf("MagForPsi   =  %5.8lf | %5.8lf per entity\n", d1,d1/rnatm);
    CkPrintf("Memory      =  %ld\n",CmiMemoryUsage());
    computeEnergies(ENERGY_FMAG, d1);
  }//endif

  if(cp_min_opt==0 && fmagPsi_total>rnatm*tol_cp_dyn){
    exitFlag=1;
    if(thisIndex.x==0 && thisIndex.y==0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Mag psi force %.10g > %.10g too large for CP dynamics. Caio! \n",
	       fmagPsi_total/rnatm,tol_cp_dyn);
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
         int heavyside = 1-(iteration-iterRotation >= 1 ? 1 : 0);
         CkPrintf("Step = %d : Step Last Rot = %d : Interval Rot = %d : Num Rot = %d : %d\n",
                   iteration,iterRotation,iteration-iterRotation,nrotation,heavyside);
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
                                           int    s_grain,
					   int    _numChunks,
					   int   _gforward,
					   int   _gbackward
					   ) 
//============================================================================
   {//begin routine
//============================================================================
//  ckout << "State G Space Constructor : "
//	<< thisIndex.x << " " << thisIndex.y << " " <<CkMyPe() << endl;
//============================================================================

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt  = sim->cp_min_opt;

//============================================================================


  istate_ind           = thisIndex.x;
  iplane_ind           = thisIndex.y;
  ees_nonlocal         = sim->ees_nloc_on;  
  numChunks            = _numChunks;
  initialized          = false;
  first_step           = 1;
  iteration            = 0;
  nrotation            = 0;
  iterRotation         = 0;
  myatom_integrate_flag= 0;
  myenergy_reduc_flag  = 0;
  isuspend_energy      = 0;
  isuspend_atms        = 0;
  forwardTimeKeep = _gforward;
  backwardTimeKeep = _gbackward;
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
  iwrite_now          = 0;
  ireset_cg           = 1;
  numReset_cg         = 0;
  exitFlag            = 0;
  iRecvRedPsi         = 1;  
  iSentRedPsi         = 1;
  iRecvRedPsiV        = 0;
  iSentRedPsiV        = 0;

  finishedCpIntegrate = 0;
  if(cp_min_opt==0){finishedCpIntegrate = 1;}// alternate entry point
  doneDoingIFFT       = false;
  allgdoneifft        = false;
  triggerNL           = false;
  NLready             = false;
  acceptedPsi         = true;    // we start out with a psi
  allAcceptedPsi      = true;    // we start out with a psi
  acceptedVPsi        = true;    // we start out with a vpsi
  acceptedLambda      = false;   // no forces yet
  needPsiV            = false;   // don't need tolerance check in first step
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
 // Symm PC accounting 
  int ourgrain    = thisIndex.x/s_grain*s_grain; 
  if(nstates == s_grain){
     AllPsiExpected=1;
  }else{ 
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
	AllLambdaExpected=(2*numGrains-1)*config.numChunksAsym;
      else
	AllLambdaExpected=numGrains*config.numChunksAsym*AllLambdaExpected;

  }//endif
  else
    {
      AllLambdaExpected*=config.numChunksAsym;

    }
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
  initGStateSlab(&gs,sizeX,size,gSpaceUnits,realSpaceUnits,s_grain,
                 thisIndex.y,thisIndex.x,len_nhc_cp,num_nhc_cp,nck_nhc_cp);

//============================================================================
// Load Balancing etc

  usesAtSync = CmiTrue;
  if(config.lbgspace){
    setMigratable(true);
  }else{
    setMigratable(false);
  }//endif

//============================================================================
// Section Reductions and SF proxy creation

  real_proxy = realSpacePlaneProxy;
  if (config.useGssInsRealP){
     ComlibAssociateProxy(&gssInstance,real_proxy);
  }//endif

  // create structure factor proxy
  if(thisIndex.x==0 && ees_nonlocal==0){
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

//============================================================================
// Register with the cache : Eric's multiple reduction schemes ensure its done
//                           before we need it.

   registrationFlag  = 1;
   eesCache *eesData = eesCacheProxy.ckLocalBranch ();
   eesData->registerCacheGSP(thisIndex.x,thisIndex.y);
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

  redPlane = 1;
  if(nchareG>1){
#ifdef _FANCY_PLANES_
    if(CkNumPes()>2*nstates){
      CkVec <int> usedVec;
      CkVec <int> peUsedByNLZ;
      CkVec <int> planeUsedByNLZ;
      FILE *fp;
      if(thisIndex.x+1==nstates){fp=fopen();}
      for(int state=0; state<thisIndex.x;state++){
        redPlane=nchareG-1;
        while(redPlane>=0){
          bool used=false;
          int thisstateplaneproc=GSImaptable.get(state,redPlane);
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

//============================================================================
// Contribute to the reduction telling main we are done

  int constructed=1;
  contribute(sizeof(int), &constructed, CkReduction::sum_int, 
	     CkCallback(CkIndex_main::doneInit(NULL),mainProxy));

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

  p|first_step; //control flags and functions reference by thread are public
  p|istate_ind;
  p|iplane_ind;
  p|registrationFlag;
  p|isuspendNLForces;
  p|initialized;
  p|istart_typ_cp;
  p|iteration;
  p|nrotation;
  p|exitFlag;
  p|isuspend_energy;
  p|isuspend_atms;
  p|cleanExitCalled;
  p|myatom_integrate_flag;
  p|myenergy_reduc_flag;
  p|finishedCpIntegrate;
  p|iRecvRedPsi;
  p|iSentRedPsi;
  p|iRecvRedPsiV;
  p|iSentRedPsiV;
  p|iwrite_now;
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
  p|needPsiV;
  p|triggerNL;
  p|NLready;
  p|doneDoingIFFT;
  p|doneNewIter;
  p|allgdoneifft;
  p|acceptedPsi;
  p|allAcceptedPsi;
  p|acceptedVPsi;
  p|acceptedLambda;
  p|itemp; // 2 temporary variables for debugging in scope
  p|jtemp;

  p|ees_nonlocal;

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
  p|sfCompSectionProxy;
  p|gpairCalcID1;
  p|gpairCalcID2;
  p| numChunks;
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
  if(p.isUnpacking())
    {
      run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_State_GSpacePlane,run),this);
    }
  RTH_Runtime_pup(run_thread,p,this);
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

  int ngridaNL      =  sim->ngrid_nloc_a;
  int ngridbNL      =  sim->ngrid_nloc_b;
  int ngridcNL      =  sim->ngrid_nloc_c;

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

      int numPoints    = 0;
      for (int j = 0; j < sortedRunDescriptors[x].size(); j++){
	numPoints += sortedRunDescriptors[x][j].length;
      }//endfor

      complex *dataToBeSent  = (complex *)fftw_malloc(numPoints*sizeof(complex));
      complex *temp          = complexPoints+ioff;
      memcpy(dataToBeSent,temp,(sizeof(complex) * numPoints));

      int numPointsV;
      complex *vdataToBeSent;
      if(istart_typ_cp>=3){
         numPointsV          = numPoints;
         vdataToBeSent       = (complex *)fftw_malloc(numPointsV*sizeof(complex));
         complex *vtemp      = vcomplexPoints+ioff;
         memcpy(vdataToBeSent,vtemp,(sizeof(complex) * numPoints));
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

      gSpacePlaneProxy(ind_state, x).initGSpace(
                                numPoints,dataToBeSent,numPointsV,vdataToBeSent,
				nx,ny,nz,ngridaNL,ngridbNL,ngridcNL,istart_typ_cp);
      fftw_free(dataToBeSent);
      fftw_free(vdataToBeSent);

      ioff += numPoints;
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
    CkPrintf("initGSpace %d.%d %d\n",thisIndex.x,thisIndex.y,size);
#endif

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  eesCache *eesData     = eesCacheProxy.ckLocalBranch ();

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

  gs.numRuns     = eesData->GspData[iplane_ind].numRuns;
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

  gs.numPoints   = eesData->GspData[iplane_ind].ncoef;;
  CkAssert(gs.numPoints == size);

  gs.ees_nonlocal        = ees_nonlocal;

  gs.packedPlaneData     = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));
  gs.packedForceData     = (complex *)fftw_malloc(gs.numFull*sizeof(complex));
  gs.packedVelData       = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));
  memcpy(gs.packedPlaneData, points, sizeof(complex)*gs.numPoints);
  bzero(gs.packedForceData,sizeof(complex)*gs.numFull);

  if(cp_min_opt==0){
    gs.packedPlaneDataScr   = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));
    gs.packedPlaneDataTemp  = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));  
    memset(gs.packedPlaneDataScr, 0, sizeof(complex)*gs.numPoints);
    memcpy(gs.packedPlaneDataTemp, points, sizeof(complex)*gs.numPoints);
  }//endif

  if(cp_min_opt==1 && cp_min_update==0){
    gs.packedPlaneDataTemp = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));  
    memcpy(gs.packedPlaneDataTemp, points, sizeof(complex)*gs.numPoints);
  }//endif

#ifdef _CP_DEBUG_SCALC_ONLY_
  if(cp_min_opt==0){
    gs.packedPlaneDataTemp2 = (complex *)fftw_malloc(gs.numPoints*sizeof(complex));  
    bzero(gs.packedPlaneDataTemp2,sizeof(complex)*gs.numPoints);
  }//endif
#endif

  // Under cp_min veldata is the conjugate gradient : always need it.
  if(istart_typ_cp>=3 && cp_min_opt==0){
    CkAssert(vsize == size);
    memcpy(gs.packedVelData, vpoints, sizeof(complex)*gs.numPoints);
  }else{
    memset(gs.packedVelData, 0, sizeof(complex)*gs.numPoints);
  }//endif

//============================================================================
// Setup gpspaceplane and particle plane

  CmiAssert(gSpacePlaneProxy(thisIndex.x, thisIndex.y).ckLocal());
  CmiAssert(particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal());
  particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal()->initKVectors(&gs);

//============================================================================
// Setup k-vector ranges, masses and zero the force overlap

  int *k_x          = eesData->GspData[iplane_ind].ka;
  int *k_y          = eesData->GspData[iplane_ind].kb;
  int *k_z          = eesData->GspData[iplane_ind].kc;
  double *coef_mass = eesData->GspData[iplane_ind].coef_mass;
  int mycoef        = eesData->GspData[iplane_ind].ncoef;
  gSpaceNumPoints   = gs.numPoints;

  if(eesData->allowedGspChares[iplane_ind]==0 || mycoef != gSpaceNumPoints){
    CkPrintf("Plane %d of state %d toasy %d %d\n",iplane_ind,thisIndex.x,
              mycoef,gSpaceNumPoints);
    CkExit();
  }//endif

  gs.setKRange(gSpaceNumPoints,k_x,k_y,k_z);

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

  if(ees_nonlocal==0){
   for(int atm=0;atm<config.numSfGrps; atm++){ //each atm
    for(int dup=0;dup<config.numSfDups;dup++){ //each dup
      if(dup==thisIndex.x){
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("GSP [%d,%d] on PE %d sending KVectors to SF[%d,%d,%d]\n",
                    thisIndex.x, thisIndex.y, CkMyPe(), atm, thisIndex.y, dup);
#endif
        sfCompProxy(atm,thisIndex.y,dup).acceptKVectors(gSpaceNumPoints, k_x, k_y, k_z);
      }//endif
    }//endfor
   }//endfor
  }//endif

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
  lambdaproxy=new CProxySection_PairCalculator[config.numChunksAsym];
  lambdaproxyother=new CProxySection_PairCalculator[config.numChunksAsym];
  psiproxy=new CProxySection_PairCalculator[config.numChunksSym];
  psiproxyother=new CProxySection_PairCalculator[config.numChunksSym];
  //need one proxy per chunk
  if(!config.gSpaceSum){
      for(int chunk=0;chunk<config.numChunksAsym;chunk++){
	  lambdaproxy[chunk]=makeOneResultSection_asym(&gpairCalcID2, 
                                                       thisIndex.x, thisIndex.y,chunk);
	  if(AllLambdaExpected/config.numChunksAsym == 2)//additional col. red. in dynamics
	    lambdaproxyother[chunk]=makeOneResultSection_asym_column(&gpairCalcID2, 
                                                        thisIndex.x, thisIndex.y,chunk);
      }//endfor chunk
      for(int chunk=0; chunk < config.numChunksSym ;chunk++){
	  psiproxy[chunk]=makeOneResultSection_sym1(&gpairCalcID1, 
                                                     thisIndex.x, thisIndex.y,chunk);
	  if(AllPsiExpected / config.numChunksSym > 1)
	    psiproxyother[chunk]=makeOneResultSection_sym2(&gpairCalcID1, 
                                                           thisIndex.x, thisIndex.y,chunk);
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

  if(iteration>0){
   if(egroupProxy.ckLocalBranch()->iteration_gsp != iteration || 
     atomsGrpProxy.ckLocalBranch()->iteration  != iteration){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Starting new iter before\n");
      CkPrintf("finishing atom integrate or iteration mismatch.\n");
      CkPrintf("iter_gsp %d iter_energy %d iter_atm %d and %d and %d\n",
 	        iteration,egroupProxy.ckLocalBranch()->iteration_gsp,
                atomsGrpProxy.ckLocalBranch()->iteration,
	       config.maxIter,cleanExitCalled);
      CkPrintf("chare %d %d\n",thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
  } //endif

  if(!acceptedVPsi){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Starting new iter before\n");
      CkPrintf("finishing Vpsi on chare %d %d\n",thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
  }//endif

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
  iRecvRedPsi      = 1;   if(cp_min_opt==0 && numRecvRedPsi>0){iRecvRedPsi  = 0;}
  iRecvRedPsiV     = 1;   if(cp_min_opt==0 && numRecvRedPsi>0){iRecvRedPsiV = 0;}
  iSentRedPsi      = 1;   if(cp_min_opt==0){iSentRedPsi  = 0;}
  iSentRedPsiV     = 1;   if(cp_min_opt==0){iSentRedPsiV = 0;}

  iteration++;   // my iteration # : not exactly in sync with other chares
                 //                  but should agree when chares meet.
  //  if(!config.launchNLeesFromRho)
  //    triggerNL=true;
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("GSP [%d,%d] StartNewIter : %d\n",thisIndex.x, thisIndex.y,iteration);
#endif

//============================================================================
// Check Load Balancing, Increment counter, set done flags equal to false.

#ifndef CMK_OPTIMIZE
  if(iteration==TRACE_ON_STEP ){traceBegin();}
  if(iteration==TRACE_OFF_STEP){traceEnd();}
#endif

    if(config.lbgspace || config.lbpaircalc ||config.lbdensity){
	if((iteration % (FIRST_BALANCE_STEP - PRE_BALANCE_STEP) == 0)  || 
           (iteration % (LOAD_BALANCE_STEP - PRE_BALANCE_STEP) == 0)){
	    LBTurnInstrumentOn();
	}//endif
    }//endif

//============================================================================
// Output psi at start of minimization for debugging

  if(iteration==1 && cp_min_opt==1){screenOutputPsi();}
#ifdef _CP_SUBSTEP_TIMING_
  if(forwardTimeKeep>0){
      double gstart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),TimeKeeperProxy);
      contribute(sizeof(double),&gstart,CkReduction::min_double, cb , forwardTimeKeep);
  }//endif
#endif


//---------------------------------------------------------------------------
}//end routine
//============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::releaseSFComputeZ() {
//==============================================================================
#ifdef _CP_DEBUG_SF_CACHE_
    CkPrintf("GSP [%d,%d] releases SFComp\n",thisIndex.x, thisIndex.y);
#endif
#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("GSP [%d,%d] releases SFComp\n",thisIndex.x, thisIndex.y);
#endif

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

#ifndef CMK_OPTIMIZE    
  double StartTime=CmiWallTimer();
#endif

// Do fft in forward direction, 1-D, in z direction
// A local function not a message : get pointer to memory for fft group

  eesCache *eesData   = eesCacheProxy.ckLocalBranch ();
  RunDescriptor *runs = eesData->GspData[iplane_ind].runs;

  fftCacheProxy.ckLocalBranch()->doStpFFTGtoR_Gchare(gs.packedPlaneData,gs.packedForceData, 
	           gs.numFull,gs.numPoints,gs.numLines,gs.numRuns,runs,gs.zdim,iplane_ind);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(GspaceFwFFT_, StartTime, CmiWallTimer());
#endif   
  NLready=true;  
 
}
//============================================================================


//============================================================================
// Send result to realSpacePlane : perform the transpose
// Force data cannot be overwritten due to all to all nature of comm.
// Until everyone sends, no one gets anything back
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CP_State_GSpacePlane::sendFFTData () {
#ifdef _CP_DEBUG_STATEG_VERBOSE_
    CkPrintf("sendfft %d.%d \n",thisIndex.x,thisIndex.y);
#endif

  complex *data_out = gs.packedForceData;
  int numLines = gs.numLines; // same amount of data to each realspace chare puppy
  int sizeZ    = gs.planeSize[1];

//============================================================================
// Do a Comlib Dance

  if (config.useGssInsRealP){gssInstance.beginIteration();}

//============================================================================
// Send your (x,y,z) to processors z.

  for(int z=0; z < sizeZ; z++) {

   // Malloc and prio the message
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
    complex *data    = msg->data;
    for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = data_out[j];}
    real_proxy(thisIndex.x, z).doFFT(msg);  // same state, realspace char[z]

   // progress engine baby
    CmiNetworkProgress();

  }//endfor

//============================================================================    
// Finish up 

  if (config.useGssInsRealP){gssInstance.endIteration();}
#ifdef _CP_SUBSTEP_TIMING_
  if(forwardTimeKeep>0)
    {
      double gend=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),TimeKeeperProxy);
      contribute(sizeof(double),&gend,CkReduction::max_double, cb , forwardTimeKeep);
    }
#endif

//----------------------------------------------------------------------
  }//end routine 
//============================================================================


//============================================================================
/**
 *
 * This is used to send data to the CP_State_GSpacePlanes, which  do the
 * inverse ffts (upon receiving data from all the corresponding
 * realSpacePlanes)
 *
 * Force cannot be overwitten because we can't receive until all chares send.
 * The beauty of all to all comm
 *
 * Forces are initialized HERE. No need to zero them etc. elsewhere.
 *
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_GSpacePlane::doIFFT(GSIFFTMsg *msg) {
//============================================================================
#ifdef _CP_SUBSTEP_TIMING_
  if(backwardTimeKeep>0)
    {
      double gstart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),TimeKeeperProxy);
      contribute(sizeof(double),&gstart,CkReduction::min_double, cb , backwardTimeKeep);
    }
#endif

#ifdef _CP_DEBUG_STATEG_VERBOSE_
    CkPrintf("doIfft %d.%d \n",thisIndex.x,thisIndex.y);
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
// Now do the IFFT in place 

#ifndef CMK_OPTIMIZE    
  double StartTime=CmiWallTimer();
#endif

  eesCache *eesData   = eesCacheProxy.ckLocalBranch ();
  RunDescriptor *runs = eesData->GspData[iplane_ind].runs;
  FFTcache *fftcache        = fftCacheProxy.ckLocalBranch();

  fftcache->getCacheMem("CP_State_GSpacePlane::doIFFT");
  complex *forcTmp = fftcache->tmpData;
  fftcache->doStpFFTRtoG_Gchare(gs.packedForceData,forcTmp,
            gs.numFull,gs.numPoints,gs.numLines,gs.numRuns,runs,gs.zdim,iplane_ind);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(GspaceBwFFT_, StartTime, CmiWallTimer());
#endif    

//============================================================================
// Finish up by multiplying the the FFT scaling factor

  double scaleFactor = -2.0/double(scProxy.ckLocalBranch()->cpcharmParaInfo->sizeX * 
                                   scProxy.ckLocalBranch()->cpcharmParaInfo->sizeY * 
                                   scProxy.ckLocalBranch()->cpcharmParaInfo->sizeZ);

  complex *forces = gs.packedForceData;
  for(int index=0; index<gs.numPoints; index++){
    forces[index] = forcTmp[index]*scaleFactor;
  }/*endfor*/
  fftcache->freeCacheMem("CP_State_GSpacePlane::doIFFT");

  CmiNetworkProgress();

//============================================================================
// Report that you are done

#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
     CkPrintf("Done doing Ifft gsp : %d %d\n",thisIndex.x,thisIndex.y);
#endif

  doneDoingIFFT = true;

//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
bool CP_State_GSpacePlane::doneNLForces(){
//==============================================================================

#ifdef _CP_DEBUG_STATEG_VERBOSE_
  if(!particlePlaneProxy(thisIndex.x,thisIndex.y).ckLocal()->doneGettingForces){
    if(thisIndex.x==0)
       CkPrintf("suspend me in gsp \n");
  }//endif
#endif
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
  if(doneDoingIFFT){
    if(isuspendNLForces == 0){
      CkPrintf("I never suspended NL forces %d %d\n",istate_ind,iplane_ind);     
      CkExit();
    }//endif
    isuspendNLForces = 0;
    RTH_Runtime_resume(run_thread);
  }//endif
#endif
//-------------------------------------------------------------------------------
  }//end routine   
//================================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptNLForcesEes(){
//================================================================================
//  I. You have called resume until these guys report. You have not started FFT
#ifdef _CP_DEBUG_NONLOC_BARRIER_
  RTH_Runtime_resume(run_thread);
#endif
//================================================================================
//  I. If the fft forces finished first, you need to resume(they are waiting).
// II. If the fft forces are not done, hang loose till they finish.
#ifndef _CP_DEBUG_NONLOC_BARRIER_
  if(doneDoingIFFT){
    RTH_Runtime_resume(run_thread);
  }//endif
#endif
//-------------------------------------------------------------------------------
  }//end routine   
//================================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::combineForcesGetEke(){
//================================================================================

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

  CP_State_ParticlePlane *pp = particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
  eesCache *eesData          = eesCacheProxy.ckLocalBranch ();
  complex *ppForces = pp->myForces;
  int ncoef       = gs.numPoints;
  int *k_x          = eesData->GspData[iplane_ind].ka;
  int *k_y          = eesData->GspData[iplane_ind].kb;
  int *k_z          = eesData->GspData[iplane_ind].kc;
  double *g2        = eesData->GspData[iplane_ind].g2;

//================================================================================
// Add forces from particle plane to forces from IFFT then zero them

#ifdef _CP_DEBUG_VKS_FORCES_
  if(thisIndex.x==0 && thisIndex.y==0){
    FILE *fp = fopen("vks_forces_s0_p0.out","w");
    int ncoef       = gs.numPoints;
    complex *forces = gs.packedForceData;
    for(int i=0;i<ncoef;i++){
      fprintf(fp,"%d %d %d : %g %g\n",k_x[i],k_y[i],k_z[i],
	      forces[i].re,forces[i].im);
    }//endfor
    fclose(fp);
  }//endif
#endif

  complex *forces = gs.packedForceData;
#ifdef _CP_DEBUG_SFNL_OFF_
  bzero(ppForces,ncoef*sizeof(complex));
#endif
#ifdef _CP_DEBUG_VKS_OFF_ 
  bzero(forces,ncoef*sizeof(complex));
  bzero(ppForces,ncoef*sizeof(complex));
#endif


#ifdef _CP_GS_DUMP_VKS_
    dumpMatrixDouble("vksBf",(double *)ppForces, 1, 
                gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    dumpMatrixDouble("forceBf",(double *)forces, 1, 
                gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_VKS_
  if(savedvksBf==NULL){ // load it
      savedvksBf= new complex[gs.numPoints];
      loadMatrixDouble("vksBf",(double *)savedvksBf, 1, 
                 gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  if(savedforceBf==NULL){ // load it
      savedforceBf= new complex[gs.numPoints];
      loadMatrixDouble("forceBf",(double *)savedforceBf, 1, 
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

  CPNONLOCAL::CP_eke_calc(ncoef,istate,forces,psi_g,k_x,k_y,k_z,g2,eke_ret,
			  config.doublePack,nkx0);
  contribute(sizeof(double), &gs.eke_ret, CkReduction::sum_double, 
	     CkCallback(printEnergyEke, NULL));
  myenergy_reduc_flag=0;

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
#ifdef _CP_DEBUG_PSI_OFF_
  iteration++;
  if(iteration==config.maxIter+1){
    int i=0;
    contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
    cleanExitCalled = 1;
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

  if(cp_min_opt==0){
    int ncoef           = gs.numPoints;
    complex *psi_g      = gs.packedPlaneData;
    complex *psi_g_tmp  = gs.packedPlaneDataTemp;
    memcpy(psi_g_tmp,psi_g,sizeof(complex)*ncoef);
  }//endif

#ifndef _CP_DEBUG_ORTHO_OFF_
  if(gs.ihave_kx0==1 && cp_min_opt==0){
    double rad2i = 1.0/sqrt(2.0);
    double rad2  = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i]   *= rad2i;}
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){force[i] *= rad2;}
  }//endif
#endif
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  for(int i=0;i<config.numChunksAsym;i++)
    CkAssert(countLambdaO[i]==0);
#endif
#ifdef _CP_GS_DUMP_LAMBDA_
    dumpMatrixDouble("lambdaBf",(double *)force, 1, 
                     gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
    dumpMatrixDouble("psiBf",(double *)psi, 1, 
                     gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  if(savedlambdaBf==NULL){ // load it
      savedlambdaBf= new complex[gs.numPoints];
      loadMatrixDouble("lambdaBf",(double *)savedlambdaBf, 1, 
                       gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  if(savedpsiBf==NULL){ // load it
      savedpsiBf= new complex[gs.numPoints];
      loadMatrixDouble("psiBf",(double *)savedpsiBf, 1, 
                        gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif

  for(int i=0;i<gs.numPoints;i++){
      if(fabs(force[i].re-savedlambdaBf[i].re)>0.0001){
	  fprintf(stderr, "GSP [%d,%d] %d element lambda  %.10g not %.10g\n",
                  thisIndex.x, thisIndex.y,i, force[i].re, savedlambdaBf[i].re);
      }//endif
      CkAssert(fabs(force[i].re-savedlambdaBf[i].re)<0.0001);
      CkAssert(fabs(force[i].im-savedlambdaBf[i].im)<0.0001);
  }//endfor
  for(int i=0;i<gs.numPoints;i++){
      CkAssert(fabs(psi[i].re-savedpsiBf[i].re)<0.0001);
      CkAssert(fabs(psi[i].im-savedpsiBf[i].im)<0.0001);
  }//endfor
#endif

  int numPoints   = gs.numPoints;
#ifndef _CP_DEBUG_ORTHO_OFF_
  int toSend = numPoints;
  startPairCalcLeft(&gpairCalcID2,toSend,psi,thisIndex.x,thisIndex.y,false);
  CmiNetworkProgress();
  startPairCalcRight(&gpairCalcID2,toSend,force,thisIndex.x, thisIndex.y);
#else
  acceptedLambda=true;
  bzero(force,sizeof(complex)*numPoints);
#endif

#ifdef _CP_DEBUG_STATEG_VERBOSE_ 
  if(thisIndex.x==0){CkPrintf("Sent Lambda %d %d\n",thisIndex.y,cleanExitCalled);}
#endif
#ifdef _CP_SUBSTEP_TIMING_
  if(backwardTimeKeep>0){
      double gend=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),TimeKeeperProxy);
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

  int cp_min_opt    = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind].ka;

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
    dumpMatrixDouble("lambdab4",(double *)force, 1, gs.numPoints*2,
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
#ifdef CMK_VERSION_BLUEGENE
#pragma disjoint(*force, *data)
      __alignx(16,force);
      __alignx(16,data);
#endif

  //---------------------------------------------------
  // B) Double Pack
  if(config.doublePack==1){
   if(cp_min_opt==1){
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
     for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       double wght  = (k_x[idest]==0 ? 0.5 : 1);
       force[idest].re -= wght*data[i].re;
       force[idest].im -= wght*data[i].im;
     }//endfor
   }else{
     if(countLambdaO[offset]<1){
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
     }else{
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
       for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
     }// coutlambda
   }//endif : cpmin
 
  }//endif : double pack

  //---------------------------------------------------
  // C) Double Pack
  if(config.doublePack==0){
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
    for(int i=0,idest=chunkoffset; i<N; i++,idest){
       force[idest].re -= 0.5*data[i].re;
       force[idest].im -= 0.5*data[i].im;
    }//endfor
  }//endif

  //---------------------------------------------------
  // D) Cleanup
  delete msg;  

//==============================================================================
// Do we have everything?

  countLambdaO[offset]++;
  if(countLambda==AllLambdaExpected){ 
    doLambda();
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

  int cp_min_opt    = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind].ka;

  complex *force    = gs.packedForceData;
  int chunksize     = gs.numPoints/config.numChunksAsym;
  int chunkoffset   = offset*chunksize; // how far into the points this contribution lies

//==============================================================================
// I) Increment the counter and do some checking

  countLambda++;  //lambda arrives in as many as 2 * numChunks reductions
  //  CkPrintf("[%d %d] accepts lambda off %d %d of %d\n", thisIndex.x, thisIndex.y,offset, countLambda, AllLambdaExpected);
/*
  if(thisIndex.y==0){
    dumpMatrixDouble("lambdab4",(double *)force, 1, gs.numPoints*2,
                     thisIndex.y,thisIndex.x,thisIndex.x,0,false);
  }//endif

  if(thisIndex.y==0){
      CkPrintf("LAMBDA [%d %d], offset %d chunkoffset %d N %d countLambdao %d\n", 
                thisIndex.x, thisIndex.y, offset, chunkoffset, N, countLambdaO[offset]);
  }//endif
*/

//=============================================================================
// (II) Add it in to our forces : Careful about offsets, doublepack and cpmin/cp

 //----------------------------------------------------------
 //A) BlueGene nonsense

#ifdef CMK_VERSION_BLUEGENE
#pragma disjoint(*force, *data)
      __alignx(16,force);
      __alignx(16,data);
#endif
 //----------------------------------------------------------
 //B) Double Pack

  if(config.doublePack==1){
   if(cp_min_opt==1){
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
     for(int i=0,idest=chunkoffset; i<N; i++,idest++){
       double wght  = (k_x[idest]==0 ? 0.5 : 1);
       force[idest].re -= wght*data[i].re;
       force[idest].im -= wght*data[i].im;
     }//endfor
   }else{
     if(countLambdaO[offset]<1){
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
        for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  = data[i]*(-1.0);}
     }else{
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
        for(int i=0,idest=chunkoffset; i<N; i++,idest++){force[idest]  += data[i]*(-1.0);}
     }//endif : off set thingy
   }//endif : cp_min_on
  }//endif : double pack

 //----------------------------------------------------------
 //C) Single pack

  if(config.doublePack==0){
#ifdef CMK_VERSION_BLUEGENE
#pragma unroll(10)
#endif
    for(int i=0,idest=chunkoffset; i<N; i++,idest){
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
    doLambda();
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

  CkAssert(countLambda==AllLambdaExpected);
  int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  complex *force = gs.packedForceData;
#ifdef _NAN_CHECK_
  for(int i=0;i<gs.numPoints ;i++)
    {
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
    memcpy(psi_g,psi_g_scr,sizeof(complex)*ncoef); // overwrite ortho with non-ortho
  }//endif
    
//==============================================================================
// Resume 

  acceptedLambda=true;
  countLambda=0;
  bzero(countLambdaO,config.numChunksAsym*sizeof(int));

//==============================================================================
// Debug

#ifdef _CP_GS_DUMP_LAMBDA_
    dumpMatrixDouble("lambdaAf",(double *)force, 1, gs.numPoints*2,
                      thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  if(savedlambdaAf==NULL){ // load it
      savedlambdaAf= new complex[gs.numPoints];
      loadMatrixDouble("lambdaAf",(double *)savedlambdaAf, 1, 
                        gs.numPoints*2,thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  for(int i=0;i<gs.numPoints;i++){
      CkAssert(fabs(force[i].re-savedlambdaAf[i].re)<0.0001);
      CkAssert(fabs(force[i].im-savedlambdaAf[i].im)<0.0001);
  }//endfor
#endif

#ifdef _CP_DEBUG_STATEG_VERBOSE_
  if(thisIndex.x==0)
   CkPrintf("doLambda %d %d\n",thisIndex.y,cleanExitCalled);
#endif

//==============================================================================
// Back to the threaded loop

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
void CP_State_GSpacePlane::writeStateDumpFile() {
//============================================================================
// Local pointers, variables and error checking

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int cp_min_opt = (sim->cp_min_opt);
  int ndump_frq  = (sim->ndump_frq);

  int ind_state  = (thisIndex.x+1);
  int ind_chare  = (thisIndex.y+1);
  int ncoef      = gSpaceNumPoints;

  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind].ka;
  int *k_y          = eesData->GspData[iplane_ind].kb;
  int *k_z          = eesData->GspData[iplane_ind].kc;
  double *coef_mass = eesData->GspData[iplane_ind].coef_mass;

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

  if(!acceptedPsi || !acceptedLambda || !acceptedVPsi){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Flow of Control Error : Attempting to write states\n");
     CkPrintf("without completing psi, vpsi and Lambda\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  int myiteration = iteration;
  if(cp_min_opt==0){myiteration=iteration-1;}

//============================================================================
// Set the file names and write the files

  iwrite_now     = 0;
  if(config.stateOutputOn==1){

    if( ((myiteration % ndump_frq)==0) || (iteration==config.maxIter) || 
          exitFlag==1){
    //------------------------------------------------------------------
    // Set the flag and tell the world you are writing
      iwrite_now = 1;
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
      gSpacePlaneProxy(thisIndex.x,redPlane).collectFileOutput(msg);
    //------------------------------------------------------------------
    // If you are not plane redPlane, you are done. Invoke the correct reduction.
      if(thisIndex.y!=redPlane){
        int i = 0;
        if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==1){
          if(myatom_integrate_flag==1 && myenergy_reduc_flag==1){
            contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
	  }//endif
          cleanExitCalled = 1;
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
     char psiName[200]; char vpsiName[200];
     sprintf(psiName, "%s/state%d.out", config.dataPathOut,ind_state);
     sprintf(vpsiName,"%s/vState%d.out",config.dataPathOut,ind_state);
     writeStateFile(npts_tot,tpsi,tvpsi,tk_x,tk_y,tk_z,cp_min_opt,
                    sizeX,sizeY,sizeZ,psiName,vpsiName,ibinary_write_opt,
                    myiteration,ind_state);
     fftw_free(tpsi); tpsi  = NULL;
     fftw_free(tvpsi);tvpsi = NULL;
     fftw_free(tk_x); tk_x  = NULL;
     fftw_free(tk_y); tk_y  = NULL;
     fftw_free(tk_z); tk_z  = NULL;
     int i=0;
     if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==1){
       if(myatom_integrate_flag==1 && myenergy_reduc_flag==1){
         contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
         cleanExitCalled = 1;
       }//endif
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

#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("integrate %d %d\n",thisIndex.y,cleanExitCalled);
#endif

//==============================================================================
// II. Local pointers    

  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind].ka;
  int *k_y          = eesData->GspData[iplane_ind].kb;
  int *k_z          = eesData->GspData[iplane_ind].kc;
  double *coef_mass = eesData->GspData[iplane_ind].coef_mass;

  int cp_min_opt     = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  int cp_min_cg      = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_cg;

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

#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif

//#define _GLENN_CHECK_DYNAMICS_
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
  psi_g = gs.packedPlaneDataTemp2;
#endif

  fictEke = 0.0; ekeNhc=0.0; potNHC=0.0;
  CPINTEGRATE::CP_integrate(ncoef,istate,iteration,forces,forcesold,psi_g,
               coef_mass,k_x,k_y,k_z,len_nhc,num_nhc,nck_nhc,fNHC,vNHC,xNHC,xNHCP,
   	       mNHC,v0NHC,a2NHC,a4NHC,kTCP,gamma_conj_grad,&fictEke,
               nkx0_red,nkx0_uni,nkx0_zero,&ekeNhc,&potNHC,degfree,degfreeNHC,
	       degFreeSplt,istrNHC,iendNHC,halfStepEvolve);
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
             CkCallback(printFictEke, NULL));

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


  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
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
      RTH_Runtime_resume(run_thread); 
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt       = sim->cp_min_opt;

  if(finishedCpIntegrate==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't sendPsi without completing integrate\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if(cp_min_opt==0 && iteration>0){
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
  allAcceptedPsi = false;

  complex *psi   = gs.packedPlaneData;
  int numPoints  = gs.numPoints;

  if(cp_min_opt==0){
     int ncoef     = gs.numPoints;
     complex *scr  = gs.packedPlaneDataScr; //save non-orthog psi
     memcpy(scr,psi,sizeof(complex)*ncoef);
  }//endif

#ifndef _CP_DEBUG_ORTHO_OFF_
  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i] *= rad2i;}
  }//endif
#endif

//==============================================================================
// Debugging 

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  for(int i=0;i<config.numChunksSym;i++){CkAssert(countPsiO[i]==0);}
#endif

#ifdef _CP_GS_DUMP_PSI_
    dumpMatrixDouble("psiBfp",(double *)psi, 1, gs.numPoints*2,
                     thisIndex.y,thisIndex.x,thisIndex.x,0,false);     
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  if(savedpsiBfp==NULL){ // load it
      savedpsiBfp= new complex[gs.numPoints];
      loadMatrixDouble("psiBfp",(double *)savedpsiBfp, 1, gs.numPoints*2,
                        thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  for(int i=0;i<gs.numPoints;i++){
      if(fabs(psi[i].re-savedpsiBfp[i].re)>0.0001){
	  fprintf(stderr, "GSP [%d,%d] %d element psi  %.10g not %.10g\n",
                  thisIndex.x, thisIndex.y,i, psi[i].re, savedpsiBfp[i].re);
      }//endif
      CkAssert(fabs(psi[i].re-savedpsiBfp[i].re)<0.0001);
      CkAssert(fabs(psi[i].im-savedpsiBfp[i].im)<0.0001);
  }//endfor
#endif

//==============================================================================
// Start the calculator

#ifndef _CP_DEBUG_ORTHO_OFF_
  startPairCalcLeft(&gpairCalcID1, numPoints, psi, thisIndex.x, thisIndex.y, false);
#else
  acceptedPsi=true;
  if((iteration==config.maxIter || exitFlag==1) && cp_min_opt==1 && 
      config.stateOutputOn==0){
      if(myatom_integrate_flag==1 && myenergy_reduc_flag==1){
        int i;
        contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
      }//endif
      cleanExitCalled = 1;
  }///endif
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt       = sim->cp_min_opt;
  if(cp_min_opt==0 && iteration>0){
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
    //memcpy(&(psi[idest]), &(data[0]), N*sizeof(complex)); //slower?
    for(int i=0; i<N; i++,idest++){psi[idest] = data[i];}
  }else{
    for(int i=0; i<N; i++,idest++){psi[idest] += data[i];}
    //    fastAdd((double *) psi,(double *)data,N*2);
  }//endif

  delete msg;

//=============================================================================
// (II) If you have got it all : Rescale it, produce some output

  countPsi++;//psi arrives in as many as 2 *numblock reductions
  countPsiO[offset]++;//psi arrives in as many as 2 
  if(countPsi==AllPsiExpected){ 
    doNewPsi();
#ifdef _CP_DEBUG_STATEG_VERBOSE_
   if(thisIndex.x==0)
    CkPrintf("aceeptpsi %d %d\n",thisIndex.y,cleanExitCalled);
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt       = sim->cp_min_opt;
  if(cp_min_opt==0 && iteration>0){
    if(iRecvRedPsi!=1 || iSentRedPsi!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, you can't acceptPsi without receiving Redpsi\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//=============================================================================
// (0) Check for Nans

//  CkPrintf("GSP [%d,%d] acceptNewPsi\n",thisIndex.x,thisIndex.y);

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
    //memcpy(&(psi[idest]), &(data[0]), N*sizeof(complex)); //slower?
    for(int i=0; i<N; i++,idest++){psi[idest] = data[i];}
  }else{
    for(int i=0; i<N; i++,idest++){psi[idest] += data[i];}
    //    fastAdd((double *) psi,(double *)data,N*2);
  }//endif

  delete msg;

//==============================================================================
// When you are done, continue

  countPsi++;         //psi arrives in as many as 2 *numblock * numgrain reductions
  countPsiO[offset]++;//psi arrives in as many as 2 * numgrain
  if(countPsi==AllPsiExpected){ 
    doNewPsi();
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

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int cp_min_opt       = sim->cp_min_opt;
  int cp_min_update    = sim->cp_min_update;

  if(cp_min_opt==0 && iteration>0){
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

#ifdef GPSI_BARRIER
  int wehaveours=1;
  contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
	     CkCallback(CkIndex_CP_State_GSpacePlane::gdonePsi(NULL),gSpacePlaneProxy));
#endif

  acceptedPsi = true;
  countPsi    = 0;
  bzero(countPsiO,config.numChunksSym*sizeof(int));
  if(gs.ihave_kx0==1){
    double rad2 = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){psi[i] *= rad2;}
  }//endif

//=============================================================================
// (B) Generate some screen output of orthogonal psi

  screenOutputPsi();

//=============================================================================
// (D) Go back to the top or exit

  if((iteration==config.maxIter || exitFlag==1)&& cp_min_opt==0){
     if(myatom_integrate_flag==1 && myenergy_reduc_flag==1){
      int i;
      contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
     }//endif
     cleanExitCalled = 1;
  }//endifw

  if((iteration==config.maxIter || exitFlag==1) && cp_min_opt==1 && 
     config.stateOutputOn==0){
     if(myatom_integrate_flag==1 && myenergy_reduc_flag==1){
       int i;
       contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
     }//endif
     cleanExitCalled = 1;
  }///endif

//=============================================================================
// (E) Debug psi

#ifdef _CP_GS_DUMP_PSI_
  dumpMatrixDouble("psiAf",(double *)psi, 1, gs.numPoints*2,thisIndex.y,thisIndex.x,
                    thisIndex.x,0,false);    
#endif

#ifdef _CP_GS_DEBUG_COMPARE_PSI_
  if(savedpsiAf==NULL){
      savedpsiAf= new complex[gs.numPoints];
      loadMatrixDouble("psiAf",(double *)savedpsiAf, 1, gs.numPoints*2,
                        thisIndex.y,thisIndex.x,thisIndex.x,0,false);    
  }//endif
  for(int i=0;i<gs.numPoints;i++){
      if(fabs(psi[i].re-savedpsiAf[i].re)>0.0001){
        fprintf(stderr, "GSP [%d,%d] %d element psi  %.10g not %.10g\n",
        thisIndex.x, thisIndex.y,i, psi[i].re, savedpsiAf[i].re);
      }//endif
      CkAssert(fabs(psi[i].re-savedpsiAf[i].re)<0.0001);
      CkAssert(fabs(psi[i].im-savedpsiAf[i].im)<0.0001);
  }//endfor
#endif

//=============================================================================
// (E) Reset psi 

  if(cp_min_opt==1 && cp_min_update==0){
    memcpy(gs.packedPlaneData,gs.packedPlaneDataTemp,
	      sizeof(complex)*gs.numPoints);
    memset(gs.packedVelData,0,sizeof(complex)*gs.numPoints);
  }//endif

//==============================================================================
// Back to the threaded loop.

  RTH_Runtime_resume(run_thread);

//----------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::sendRedPsiV(){
//==============================================================================
// I) Local Pointers

  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  double *coef_mass = eesData->GspData[iplane_ind].coef_mass;

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

  CPcharmParaInfo *sim       = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
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
         gSpacePlaneProxy(thisIndex.x,irecv).acceptRedPsiV(msg);
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
// I send the stuff and I need a new velocity

 iSentRedPsiV  = 1;
 acceptedVPsi = false;

//==============================================================================
  }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptRedPsiV(GSRedPsiMsg *msg) {
//==============================================================================

  CPcharmParaInfo *sim       = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  int ncoef         = msg->size;
  int isend         = msg->senderIndex;
  complex *msgData  = msg->data;
  complex *recvData = gs.packedRedPsiV;

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
    countRedPsiV = 0;
    iRecvRedPsiV  = 1;
    if(jtemp!=gs.nkx0_red){
      CkPrintf("Error in GSchare recv cnt %d %d : %d %d\n",thisIndex.x,thisIndex.y,
 	                                                   gs.nkx0_red,jtemp);
      CkExit();
    }//endif
    // If sent before I received then I resume
    if(iSentRedPsiV == 1){
        RTH_Runtime_resume(run_thread);
    }//endif
  }//endif

//-----------------------------------------------------------------------------
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::doneRedPsiVIntegrate() {

  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  double *coef_mass = eesData->GspData[iplane_ind].coef_mass;
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

  CPcharmParaInfo *sim       = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  RedundantCommPkg *RCommPkg = sim->RCommPkg;

  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't sendPsiV without\n");
    CkPrintf("sending the Redundant psi values around\n");
    CkPrintf("chare %d %d : finished %d %d : %d %d\n",
	     thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  acceptedVPsi = false;

//==============================================================================
//
  nrotation++;
  iterRotation = iteration+1;

  int ncoef     = gs.numPoints;
  complex *data = gs.packedVelData;
  complex *scr  = gs.packedPlaneDataScr;  // replace no-ortho psi 
  complex *psi  = gs.packedPlaneData;     // by orthonormal psi when norb rotating
  memcpy(scr,psi,sizeof(complex)*ncoef);

  if(gs.ihave_kx0==1){
    double rad2i = 1.0/sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){data[i] *= rad2i;}
  }//endif

  int numPoints = gs.numPoints;
  startPairCalcLeft(&gpairCalcID1,numPoints,data,thisIndex.x,thisIndex.y,true);

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

  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't accetPsiV without\n");
    CkPrintf("sending the Redundant psi values around\n");
    CkPrintf("chare %d %d : finished %d %d : %d %d\n",
	     thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//=============================================================================
// (II) Unpack the contribution to newpsi (orthonormal psi)

  int idest=chunkoffset;
  if(countVPsiO[offset]<1){
    // memcpy(&(vpsi[idest]), &(data[0]), N*sizeof(complex));//slower?
    for(int i=0; i<N; i++,idest++){vpsi[idest] = data[i];}
  }else{
    for(int i=0; i<N; i++,idest++){vpsi[idest] += data[i];}
  }//endif  

  delete msg;

//=============================================================================
// (III) When all has arrive, onward to victory

  countVPsi++;         //psi arrives in as many as 2 reductions
  countVPsiO[offset]++;//psi arrives in as many as 2 reductions

  if(countVPsi==AllPsiExpected){ 
    doNewPsiV();
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

  complex *vpsi   = gs.packedVelData;
  int chunksize   = gs.numPoints/config.numChunksSym;
  int chunkoffset = offset*chunksize;; // how far into the points this contribution lies

  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't accetPsiV without\n");
    CkPrintf("sending the Redundant psi values around\n");
    CkPrintf("chare %d %d : finished %d %d : %d %d\n",
	     thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//=============================================================================
// (I) Unpack the contribution to newpsi (orthonormal psi)

  int idest=chunkoffset;
  if(countVPsiO[offset]<1){
    // memcpy(&(vpsi[idest]), &(data[0]), N*sizeof(complex));//slower?
    for(int i=0; i<N; i++,idest++){vpsi[idest] = data[i];}
  }else{
    for(int i=0; i<N; i++,idest++){vpsi[idest] += data[i];}
  }//endif
  
  delete msg;

//=============================================================================
// (II) Continue

  countVPsi++;//psi arrives in as many as 2 reductions
  countVPsiO[offset]++;//psi arrives in as many as 2 reductions

  if(countVPsi==AllPsiExpected){ 
    doNewPsiV();
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

  CkAssert(countVPsi==AllPsiExpected); 

  if(iRecvRedPsiV!=1 || iSentRedPsiV!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude, you can't doNewPsiV without\n");
    CkPrintf("sending the Redundant psi values around\n");
    CkPrintf("chare %d %d : finished %d %d : %d %d\n",
	     thisIndex.x,thisIndex.y,iRecvRedPsiV,iSentRedPsiV,numRecvRedPsi,gs.nkx0_red);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//=============================================================================
// (I) Reset counters and rescale the kx=0 stuff

  acceptedVPsi    = true;
  needPsiV        = false;
  countVPsi       = 0;

  complex *vpsi = gs.packedVelData;

  for(int i=0;i<config.numChunksSym;i++){countVPsiO[i]=0;}
  if(gs.ihave_kx0==1){
    double rad2 = sqrt(2.0);
    for(int i=gs.kx0_strt; i<gs.kx0_end; i++){vpsi[i] *= rad2;}
  }//endif

//=============================================================================
// II) A Barrier for debugging

  /* debugging barrier
     int wehaveours=1;
     contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
     CkCallback(CkIndex_CP_State_GSpacePlane::gdonePsiV(NULL),gSpacePlaneProxy));
  */

//=============================================================================
// III) Replace by finite difference until update is better

 CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
 eesCache *eesData    = eesCacheProxy.ckLocalBranch ();

 double *coef_mass    = eesData->GspData[iplane_ind].coef_mass;
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
// III) Back to the threaded loop

  RTH_Runtime_resume(run_thread);

//----------------------------------------------------------------------------
   }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::screenOutputPsi(){
//==============================================================================
#ifdef _CP_DEBUG_STATEG_VERBOSE_
  if(thisIndex.x==0){CkPrintf("output %d %d\n",thisIndex.y,cleanExitCalled);}
#endif

  eesCache *eesData = eesCacheProxy.ckLocalBranch ();
  int *k_x          = eesData->GspData[iplane_ind].ka;
  int *k_y          = eesData->GspData[iplane_ind].kb;
  int *k_z          = eesData->GspData[iplane_ind].kc;

  int cp_min_opt    = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  int nstates       = scProxy.ckLocalBranch()->cpcharmParaInfo->nstates;

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
	  PRINT_LINE_DASH;
	  CkPrintf(" Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],psi[i].re,psi[i].im);
          if(cp_min_opt==0){
            double vre=vpsi[i].re;
            double vim=vpsi[i].im;
 	    CkPrintf("VPsi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],vre,vim);
	  }//endif
	  PRINT_LINE_DASH;
	}//endif
	if(k_x[i]==2 && k_y[i]==1 && k_z[i]==3){
          double vre=vpsi[i].re;
          double vim=vpsi[i].im;
	  PRINT_LINE_DASH;
	  CkPrintf(" Psi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],psi[i].re,psi[i].im);
          if(cp_min_opt==0){
 	    CkPrintf("VPsi[is=%d ka=%d kb=%d kc=%d] : %.15g %.15g\n",
		   gs.istate_ind+1,k_x[i],k_y[i],k_z[i],vre,vim);
	  }//endif
	  PRINT_LINE_DASH;
	}//endif
      }//endfor
    }//endif
  }//endif
#endif

//==============================================================================
// II) Tell the world you are done with the output

#ifdef _CP_DEBUG_STATE_GPP_VERBOSE_
  CkPrintf("GSP [%d,%d] screenwrite: %d\n",thisIndex.x, thisIndex.y,iteration);
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
//
//  enable config.localAtomBarrier to trigger that behavior
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptAtoms(GSAtmMsg *msg) {
//==============================================================================
// Do not delete msg. Its a no keep.
//==============================================================================
// Flip my flag, check my iteration

   myatom_integrate_flag=1;
   if(atomsGrpProxy.ckLocalBranch()->iteration != iteration){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error GSP::acceptatoms : \n");
      CkPrintf("Iteration mismatch between atoms and g-space planes\n");
      CkPrintf("suspend_atms %d suspend_energy %d,iteration_gsp %d iteration_atm %d\n",
                isuspend_atms,isuspend_energy,
                iteration,atomsGrpProxy.ckLocalBranch()->iteration);
      CkPrintf("state %d plane %d proc %d\n",thisIndex.x,thisIndex.y,CkMyPe());
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
#ifdef _CP_DEBUG_WARN_SUSPEND_
   CkPrintf("Atoms on proc %d GSP chare %d %d : %d %d\n",
             CkMyPe(),thisIndex.x,thisIndex.y,myatom_integrate_flag,iteration);
#endif

//==============================================================================
// Exit if the time is right

   if(myenergy_reduc_flag==1 && cleanExitCalled==1){
     int i;
     contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
   }//endif

//==============================================================================
// Lift the suspension if my energy is ready to rock and roll

   if(isuspend_atms==1){ // I suspended to wait for atoms, resume me baby.
     isuspend_atms=0;
     if(isuspend_energy==0){RTH_Runtime_resume(run_thread);}
   }//endif

//-----------------------------------------------------------------------------
  }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::acceptEnergy(GSAtmMsg *msg) {
//==============================================================================
//   do not delete message : Its a nokeep 
//==============================================================================
// Flip my flag, check my iteration

   myenergy_reduc_flag=1;
   if(egroupProxy.ckLocalBranch()->iteration_gsp != iteration){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error : Iteration\n");
      CkPrintf("mismatch between energy and g-space planes\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
   }//endif
#ifdef _CP_DEBUG_WARN_SUSPEND_
   CkPrintf("Energy on proc %d GSP chare %d %d : %d %d\n",
             CkMyPe(),thisIndex.x,thisIndex.y,myenergy_reduc_flag,iteration);
#endif

//==============================================================================
// Exit if the time is right

   if(myatom_integrate_flag==1 && cleanExitCalled==1){
     int i;
     contribute(sizeof(int),&i,CkReduction::sum_int,CkCallback(cleanExit,NULL));
   }//endif

//==============================================================================
// Lift the suspension if my atoms are ready to rock and roll

   if(isuspend_energy==1){ // I suspended to wait for energy, resume me baby.
     isuspend_energy=0;
     if(isuspend_atms==0){RTH_Runtime_resume(run_thread);}
   }//endif

//------------------------------------------------------------------------------
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

//==============================================================================
// reset commlib proxies

  if(config.useGssInsRealP){ComlibResetProxy(&real_proxy);}
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int cp_min_opt = sim->cp_min_opt;

//==============================================================================
// reset lambda PC proxies

  CkMulticastMgr *mcastGrp = 
        CProxy_CkMulticastMgr(gpairCalcID2.mCastGrpId[thisIndex.y]).ckLocalBranch();         
  for(int chunk=0;chunk<config.numChunksAsym;chunk++){
      mcastGrp->resetSection(lambdaproxy[chunk]);
      setResultProxy(&lambdaproxy[chunk], thisIndex.x, gpairCalcID2.GrainSize, 
		     gpairCalcID2.mCastGrpId[thisIndex.y],true,
                     CkCallback(CkIndex_Ortho::lbresume(NULL),orthoProxy));
      if(cp_min_opt==0){
	  mcastGrp->resetSection(lambdaproxyother[chunk]);
	  setResultProxy(&lambdaproxyother[chunk], thisIndex.x, gpairCalcID2.GrainSize, 
			 gpairCalcID2.mCastGrpId[thisIndex.y], true,
  		         CkCallback(CkIndex_Ortho::lbresume(NULL),orthoProxy));
      }//endif
  }//endfo
  gpairCalcID2.resetProxy();

//==============================================================================
// turn off : output

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
     CProxy_CkMulticastMgr(gpairCalcID2.mCastGrpId[thisIndex.y]).ckLocalBranch();         

  for(int chunk=0;chunk<config.numChunksSym;chunk++){
      mcastGrp->resetSection(psiproxy[chunk]);
      setResultProxy(&psiproxy[chunk],thisIndex.x,gpairCalcID1.GrainSize, 
                     gpairCalcID1.mCastGrpId[thisIndex.y],true, 
                     CkCallback(CkIndex_Ortho::lbresume(NULL),orthoProxy));
      if(AllPsiExpected>1){
	mcastGrp->resetSection(psiproxyother[chunk]);
	setResultProxy(&psiproxyother[chunk], thisIndex.x, gpairCalcID1.GrainSize, 
		       gpairCalcID1.mCastGrpId[thisIndex.y],true,
                       CkCallback(CkIndex_Ortho::lbresume(NULL),orthoProxy));
      }//endif
  }//endfor

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
//==============================================================================

  needPsiV     = true;
  // when everyone is ready, restart Ortho's backward path
  int foo=1;
  contribute(sizeof(int), &foo, CkReduction::min_int, 
             CkCallback(CkIndex_Ortho::resumeV(NULL), orthoProxy));

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
bool CP_State_GSpacePlane::weneedPsiV() {
   if(needPsiV){acceptedVPsi=false;}
   return needPsiV;
}
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


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void CP_State_GSpacePlane::startNLEes(bool local,int iteration_loc){
//==============================================================================
// check and set constants

  if(!local){triggerNL=true;}
  if(iteration_loc!=iteration){CkPrintf("startNLees broken\n");CkExit();}

//==============================================================================
// I) Make sure we don't start this before we're ready

  if(NLready && triggerNL){
      //CkPrintf("GS[%d,%d] triggering NL\n");
#define _NLEES_PRIO_START_
#ifdef _NLEES_PRIO_START_OFF_
      CP_State_ParticlePlane *pp = particlePlaneProxy(thisIndex.x, thisIndex.y).ckLocal();
      pp->startNLEes(iteration);
#else
      NLDummyMsg *msg = new(8*sizeof(int)) NLDummyMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = config.sfpriority;
      msg->iteration=iteration;
      particlePlaneProxy(thisIndex.x, thisIndex.y).lPrioStartNLEes(msg);
#endif
      triggerNL=false;
      NLready=false;
  }//endif

//-----------------------------------------------------------------------------
  }//end routine
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
