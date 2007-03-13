//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file groups.C
 * 
 *           Processor group class Functions : Atoms and parainfo
 *
 */
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

//==============================================================================
#include "charm++.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include <math.h>
#include "fftCacheSlab.h"
#include "eesCache.h"
#include "CP_State_Plane.h"

//----------------------------------------------------------------------------
#define CHARM_ON
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "../../src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"

//----------------------------------------------------------------------------

extern CProxy_EnergyGroup          egroupProxy;
extern Config                      config;
extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern CProxy_AtomsGrp             atomsGrpProxy;
extern CProxy_EnergyGroup          egroupProxy;
extern CProxy_StructFactCache      sfCacheProxy;
extern CProxy_eesCache             eesCacheProxy;

//----------------------------------------------------------------------------

void IntegrationComplete(void *, void *);

//#define _CP_DEBUG_PSI_OFF_
//#define _CP_ENERGY_GRP_VERBOSE_
//#define _CP_DEBUG_ATMS_
//#define _CP_DEBUG_ATMS_EXIT_

//==============================================================================




//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** Constructor
 *
 *
 */
//==============================================================================
AtomsGrp::AtomsGrp(int n, int n_nl, int len_nhc_, int iextended_on_,int cp_min_opt_,
                   int cp_wave_opt_, int isokin_opt_,double kT_, Atom* a, AtomNHC *aNHC){
//==============================================================================
// parameters, options and energies

    natm            = n;
    natm_nl         = n_nl;
    len_nhc         = len_nhc_;
    iextended_on    = iextended_on_;
    cp_min_opt      = cp_min_opt_;
    cp_wave_opt     = cp_wave_opt_;
    isokin_opt      = isokin_opt_;
    iteration       = 0;
    kT              = kT_;
    pot_ewd_rs      = 0.0;
    eKinetic        = 0.0;
    eKineticNhc     = 0.0;
    potNhc          = 0.0;    
    countAtm        = 0;

//==============================================================================
// Initial positions, forces, velocities 

    atoms           = new Atom[natm];
    atomsNHC        = new AtomNHC[natm];
    CmiMemcpy(atoms, a, natm * sizeof(Atom));           // atoms has no vectors
    for(int i=0;i<natm;i++){atomsNHC[i].Init(&aNHC[i]);}

    zeroforces();
    if(iextended_on==1 && cp_min_opt==0){
       zeronhc();
    }//endif

//==============================================================================
// A copy of the atoms for fast routines

    fastAtoms.natm = natm;
    fastAtoms.q    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.x    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.y    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.z    = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.fx   = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.fy   = (double *)fftw_malloc(natm*sizeof(double));
    fastAtoms.fz   = (double *)fftw_malloc(natm*sizeof(double));

    copySlowToFast();
    double *qq = fastAtoms.q;
    for(int i=0;i<natm;i++){qq[i]=atoms[i].q;}

//-----------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* Destructor
 *
 *
 */
//==============================================================================
AtomsGrp::~AtomsGrp(){
     delete [] atoms;
     delete [] atomsNHC;
     fftw_free(fastAtoms.x);
     fftw_free(fastAtoms.y);
     fftw_free(fastAtoms.z);
     fftw_free(fastAtoms.fx);
     fftw_free(fastAtoms.fy);
     fftw_free(fastAtoms.fz);
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** Constructor
 *
 *
 */
//==============================================================================
CPcharmParaInfoGrp::CPcharmParaInfoGrp(CPcharmParaInfo &sim){
	cpcharmParaInfo = new CPcharmParaInfo(sim);
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* Destructor
 *
 *
 */
//==============================================================================
CPcharmParaInfoGrp::~CPcharmParaInfoGrp(){
	delete cpcharmParaInfo;
}
//==============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::StartRealspaceForces(){
//==========================================================================
// Get the real space atom forces


   int myid   = CkMyPe();
   int nproc  = CkNumPes();
   pot_ewd_rs = 0.0;

#ifndef  _CP_DEBUG_PSI_OFF_
   if(myid<natm-1){
     CPRSPACEION::CP_getionforce(natm,&fastAtoms,myid,nproc,&pot_ewd_rs);
   }//endif
#endif

#ifdef _CP_DEBUG_ATMS_
   CkPrintf("GJM_DBG: calling contribute atm forces %d\n",myid);
#endif
   contributeforces();
   
//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::contributeforces(){
//==========================================================================

  int i,j;
  int myid     = CkMyPe();
  double *ftot = (double *)malloc((3*natm+1)*sizeof(double));

  copyFastToSlow();
  for(i=0,j=0; i<natm; i++,j+=3){
    ftot[j]   = atoms[i].fx;
    ftot[j+1] = atoms[i].fy;
    ftot[j+2] = atoms[i].fz;
  }//endfor
  ftot[3*natm]=pot_ewd_rs;

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: inside contribute forces %d : %d\n",myid,natm);
#endif
  CkCallback cb(CkIndex_AtomsGrp::recvContribute(NULL), atomsGrpProxy);
  contribute((3*natm+1)*sizeof(double),ftot,CkReduction::sum_double,cb);

//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::recvContribute(CkReductionMsg *msg) {
//==========================================================================
// Local pointers

  int i,j;
  double *ftot      = (double *) msg->getData();

  EnergyGroup *eg   = egroupProxy.ckLocalBranch();

  int nproc         = CkNumPes();
  int myid          = CkMyPe();
  int output_on     = config.atmOutputOn;

//============================================================
// Copy out the reduction of energy and forces

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: inside recv forces %d : %d\n",myid,natm);
#endif
  double pot_ewd_rs = ftot[3*natm];
  double fmag = 0.0;
  for(i=0,j=0;i<natm;i++,j+=3){
    atoms[i].fx = ftot[j];
    atoms[i].fy = ftot[j+1];
    atoms[i].fz = ftot[j+2];
    fmag += (ftot[j]*ftot[j]+ftot[j+1]*ftot[j+1]+ftot[j+2]*ftot[j+2]);
#ifdef _CP_DEBUG_ATMS_
    if(myid==0){
      CkPrintf("%d : %g %g %g\n",i,atoms[i].fx,atoms[i].fy,atoms[i].fz);
    }//endif
#endif
  }//endfor
  fmag /= (double)(3*natm);
  fmag  = sqrt(fmag);
  delete msg;

#ifdef _CP_DEBUG_ATMS_EXIT_
  if(myid==0){CkExit();}
#endif

//=================================================================
// Model forces 

#ifdef  _CP_DEBUG_PSI_OFF_
  double omega  = (0.0241888/15.0); // 15 fs^{-1}
  double omega2 = omega*omega;
  double pot_harm = 0.0;
  if(iteration==0){
    for(i=0;i<natm;i++){
      atoms[i].x = 0.0;
      atoms[i].y = 0.0;
      atoms[i].z = 0.0;
    }//endfor
  }//endif
  for(i=0;i<natm;i++){
    atoms[i].fx = -omega2*atoms[i].x*atoms[i].m;
    atoms[i].fy = -omega2*atoms[i].y*atoms[i].m;;
    atoms[i].fz = -omega2*atoms[i].z*atoms[i].m;;
    pot_harm   += (atoms[i].m*omega2*(atoms[i].x*atoms[i].x+
                                      atoms[i].y*atoms[i].y+
                                      atoms[i].z*atoms[i].z));
  }//endfor
  pot_harm *= 0.5;
#endif

//=================================================================
// Check the forces

#ifdef _GLENN_CHECK_FORCES
   FILE *fp = fopen("forces.save","w");
   for(i=0;i<natm;i++){
     fprintf(fp,"%.12g %.12g %.12g\n",atoms[i].fx,atoms[i].fy,atoms[i].fz);
   }
   fclose(fp);
   CkExit();
#endif

//============================================================
// Integrate the atoms

#ifdef _CP_DEBUG_ATMS_
  CkPrintf("GJM_DBG: Before atom integrate %d : %d\n",myid,natm);
#endif

   double eKinetic_loc   = 0.0;
   double eKineticNhc_loc= 0.0;
   double potNhc_loc     = 0.0;
   int iwrite_atm        = 0;
   int myoutput_on       = output_on;

   //   CkPrintf("entering atom integrate!\n");
   if(iteration+1>config.maxIter){myoutput_on = 0;}

   int div     = (natm / nproc);
   int rem     = (natm % nproc);
   int natmStr = div*myid;
   int natmNow = div;
   if(myid>=rem){natmStr += rem;}
   if(myid< rem){natmStr += myid;}
   if(myid< rem){natmNow++;}
   int natmEnd = natmNow+natmStr;

   ATOMINTEGRATE::ctrl_atom_integrate(iteration,natm,len_nhc,cp_min_opt,
                    cp_wave_opt,iextended_on,atoms,atomsNHC,myid,
                    &eKinetic_loc,&eKineticNhc_loc,&potNhc_loc,&iwrite_atm,
                    myoutput_on,natmNow,natmStr,natmEnd);

//============================================================
  // Debug output

#ifdef  _CP_DEBUG_PSI_OFF_
   double etot_atm;
   if(isokin_opt==0){
     etot_atm = eKinetic_loc+eKineticNhc_loc+potNhc_loc+pot_harm;
     CkPrintf("iteration %d : tot energy %.12g on %d\n",iteration,etot_atm,myid);
   }else{
     etot_atm = eKineticNhc_loc+potNhc_loc;
     CkPrintf("iteration %d : tot energy %.12g %.12g on %d\n",iteration,etot_atm,
                                                    (eKinetic_loc+pot_harm),myid);
   }//endif
#endif

//============================================================
// Tuck things that can be tucked.


  eg->estruct.eewald_real     = pot_ewd_rs;  
  eg->estruct.fmag_atm        = fmag;

//============================================================
// Get ready for the next iteration : 
//      zero forces, outputAtmEnergy, atomsDone

  zeroforces();

  if(cp_wave_opt==0 && cp_min_opt==0){
    sendAtoms(eKinetic_loc,eKineticNhc_loc,potNhc_loc,natmNow,natmStr,natmEnd);
  }else{
    eKinetic                    = 0.0;
    eKineticNhc                 = 0.0;
    potNhc                      = 0.0;
    eg->estruct.eKinetic_atm    = 0.0;
    eg->estruct.eKineticNhc_atm = 0.0;
    eg->estruct.potNhc_atm      = 0.0;

    copySlowToFast();
    outputAtmEnergy();

    i=0;
    CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),atomsGrpProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

//-------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
// Atom energy output
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::outputAtmEnergy() {
//==========================================================================

  EnergyGroup *eg       = egroupProxy.ckLocalBranch();
  int myid              = CkMyPe();  
  double eKinetic       = eg->estruct.eKinetic_atm;
  double eKineticNhc    = eg->estruct.eKineticNhc_atm;
  double fmag           = eg->estruct.fmag_atm;
  double pot_ewd_rs_now = eg->estruct.eewald_real;
  double potNhc         = eg->estruct.potNhc_atm;
  double free_atm       = 3*((double)natm);

  if(myid==0){
     CkPrintf("EWALD_REAL  = %5.8lf\n",pot_ewd_rs_now);
     if(cp_min_opt==0){
        CkPrintf("atm eKin    = %5.8lf\n",eKinetic);
        CkPrintf("atm Temp    = %5.8lf\n",(2.0*eKinetic*BOLTZ/free_atm));
        if(iextended_on==1){
          double free_Nhc;
          if(isokin_opt==0){
            free_Nhc    = free_atm*((double)len_nhc);
          }else{
            free_Nhc    = free_atm*((double)(len_nhc-1));
	  }//endif
          CkPrintf("atm eKinNhc = %5.8lf\n",eKineticNhc);
          CkPrintf("atm TempNHC = %5.8lf\n",(2.0*eKineticNhc*BOLTZ/free_Nhc));
          CkPrintf("atm potNHC  = %5.8lf\n",potNhc);
	}//endif
     }else{
        CkPrintf("atm fmag    = %5.8lf\n",fmag);
     }//endif
  }//endif

//-------------------------------------------------------------------------
   }//end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void AtomsGrp::atomsDone(CkReductionMsg *msg) {
//==========================================================================
  delete msg;
  atomsDone();
}
//==========================================================================


//==========================================================================
// Needs to have each proc invoke directly acceptatoms method of the
// gspaceplanes which are mapped to it. Without migration, we have that map
// at startup. With migration, one must write an enroll/dropout routine.
// All 
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::atomsDone() {
//==========================================================================
// Increment iteration counters

 int myid = CkMyPe();

 EnergyGroup *eg             = egroupProxy.ckLocalBranch();
 iteration++;
 eg->estruct.iteration_atm   = iteration;
 eg->iteration_atm           = iteration;

#ifdef _DEBUG_CHECK_ATOM_COMM_
 char fname[100];
 sprintf(fname,"atoms.out.%d.%d",iteration,myid);
 FILE *fp = fopen(fname,"w");
 for(int i=0;i<natm;i++){
   fprintf(fp,"%.10g %.10g %.10g\n",atoms[i].x,atoms[i].y,atoms[i].z);
 }//endfor
 fclose(fp);
#endif

//==========================================================================
// Use the cool new data caching system to say we're done.

 if(config.localAtomBarrier){

   eesCache *eesData = eesCacheProxy.ckLocalBranch ();
   int *indState     = eesData->gspStateInd;
   int *indPlane     = eesData->gspPlaneInd;
   int ngo           = eesData->nchareGSPProcT;

   GSAtmMsg *msg = new  GSAtmMsg;
   for(int i=0; i<ngo; i++){
     int iadd = gSpacePlaneProxy(indState[i],indPlane[i]).ckLocal()->registrationFlag;
     if(iadd!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("atom : Bad registration cache flag on proc %d %d %d %d\n",
                myid,iadd,indState[i],indPlane[i]);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
     }//endif
     gSpacePlaneProxy(indState[i],indPlane[i]).ckLocal()->acceptAtoms(msg); 
   }//endfor
   

 }else{

   if(myid==0){
      GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = config.sfpriority-10;
      gSpacePlaneProxy.acceptAtoms(msg);
   }//endif
 
 }//endif

//--------------------------------------------------------------------------
   }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void AtomsGrp::sendAtoms(double eKinetic_loc,double eKineticNhc_loc,double potNhc_loc,
                         int natmNow,int natmStr,int natmEnd){
//==========================================================================
// Malloc the message

  int nsize    = 9*natmNow+3;
  AtomMsg *msg = new (nsize,8*sizeof(int)) AtomMsg;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  *(int*)CkPriorityPtr(msg) = config.sfpriority-10;

  double *atmData = msg->data;
  msg->nsize      = nsize;
  msg->natmStr    = natmStr;
  msg->natmEnd    = natmEnd;

//==========================================================================
// pack atom positions : new for use : old for output

  for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
    atmData[(j)  ]=atoms[i].x;
    atmData[(j+1)]=atoms[i].y;
    atmData[(j+2)]=atoms[i].z;
    atmData[(j+3)]=atoms[i].xold;
    atmData[(j+4)]=atoms[i].yold;
    atmData[(j+5)]=atoms[i].zold;
    atmData[(j+6)]=atoms[i].vxold;
    atmData[(j+7)]=atoms[i].vyold;
    atmData[(j+8)]=atoms[i].vzold;
  }//endfor

//==========================================================================
// pack the 3 energies

   atmData[(nsize-3)] = eKinetic_loc;
   atmData[(nsize-2)] = eKineticNhc_loc;
   atmData[(nsize-1)] = potNhc_loc;

//==========================================================================
// Send the message to everyone in the group

  atomsGrpProxy.acceptAtoms(msg);

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void AtomsGrp::acceptAtoms(AtomMsg *msg) {
//==========================================================================

  AtomsGrp *ag      = atomsGrpProxy.ckLocalBranch();
  double *atmData   = msg->data;
  int    nsize      = msg->nsize;
  int    natmStr    = msg->natmStr;
  int    natmEnd    = msg->natmEnd;

//==========================================================================
// unpack atom position and velocity

  for(int i=natmStr,j=0;i<natmEnd;i++,j+=9){
    atoms[i].x     = atmData[(j)  ];
    atoms[i].y     = atmData[(j+1)];
    atoms[i].z     = atmData[(j+2)];
    atoms[i].xold  = atmData[(j+3)];
    atoms[i].yold  = atmData[(j+4)];
    atoms[i].zold  = atmData[(j+5)];
    atoms[i].vxold = atmData[(j+6)];
    atoms[i].vyold = atmData[(j+7)];
    atoms[i].vzold = atmData[(j+8)];
  }//endfor

//==========================================================================
// unpack energy

  countAtm++;

  if(countAtm==1){
    eKinetic    = 0;
    eKineticNhc = 0;
    potNhc      = 0;
  }//endif

  eKinetic    += atmData[(nsize-3)];
  eKineticNhc += atmData[(nsize-2)];  
  potNhc      += atmData[(nsize-1)];

//==========================================================================
// Delete the message

  delete msg;

//==========================================================================
// Copy to the fast vectors and phone home

  if(countAtm==CkNumPes()){
     countAtm = 0;

     EnergyGroup *eg             = egroupProxy.ckLocalBranch();
     eg->estruct.eKinetic_atm    = eKinetic;
     eg->estruct.eKineticNhc_atm = eKineticNhc;
     eg->estruct.potNhc_atm      = potNhc;

     copySlowToFast();
     outputAtmEnergy();

     //everybody has to have received all the atoms before continuing : not just me
     int i=0;
     CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),atomsGrpProxy);
     contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
//Energy group that can retrieve the energies from
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyGroup::EnergyGroup () {
    iteration_gsp = 0;
    iteration_atm = 0;

    //non local
    estruct.enl = 0;

    //local external energy
    estruct.eext = 0;
    estruct.eke = 0;

    //hartree energy
    estruct.ehart = 0;
    
    //ion-ion
    estruct.eewald_recip = 0;
    estruct.eewald_real = 0;
    
    //exchange correlation
    estruct.egga = 0;
    estruct.eexc = 0;
    estruct.fmagPsi = 0;

    //CP Fict KE
    estruct.fictEke = 0;

    // total electronic part
    estruct.totalElecEnergy = 0; // needs ewald_real to be physical
    estruct.iteration_gsp = 0;

 
    // atm stuff
    estruct.eKinetic_atm    = 0;    // classical kinetic energy
    estruct.eKineticNhc_atm = 0; // NHC kinetic energy
    estruct.potNhc_atm      = 0;      // NHC pot energy
    estruct.fmag_atm        = 0;        // magnitude of atm forces
    estruct.iteration_atm   = 0;

//-------------------------------------------------------------------------
  } //end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void EnergyGroup::updateEnergiesFromGS(EnergyStruct &es) {
//==========================================================================

      estruct.enl          = es.enl;
      estruct.eke          = es.eke;
      estruct.eext         = es.eext;
      estruct.ehart        = es.ehart;
      estruct.eewald_recip = es.eewald_recip;
      estruct.egga         = es.egga;
      estruct.eexc         = es.eexc;
      estruct.fictEke      = es.fictEke;
      estruct.totalElecEnergy  = es.totalElecEnergy;
      estruct.fmagPsi      = es.fmagPsi;
      estruct.iteration_gsp= es.iteration_gsp;
      iteration_gsp        = es.iteration_gsp;
#ifdef _DEBUG_ESTRUCT_
       CkPrintf("Energies received %lf, %lf, %lf, %lf, %lf\n", 
                 estruct.enl,estruct.eke,estruct.eext,estruct.ehart, 
                 estruct.eewald_recipt,estruct.egga,estruct.eexc,
                 estruct.fictEke,estruct.totalEnergy);
#endif
    int i=0;
    CkCallback cb(CkIndex_EnergyGroup::energyDone(NULL),egroupProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void EnergyGroup::energyDone(CkReductionMsg *msg) {
//==========================================================================
   delete msg;
   energyDone();
}
//==========================================================================


//==========================================================================
// Needs to have each proc invoke directly acceptenergy method of the
// gspaceplanes which are mapped to it. Without migration, we have that map
// at startup. With migration, one must write an enroll/dropout routine.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void EnergyGroup::energyDone(){
//==========================================================================
// Use the cool new data caching system

 int myid          = CkMyPe();
 if(config.localAtomBarrier){

   eesCache *eesData = eesCacheProxy.ckLocalBranch ();
   int *indState     = eesData->gspStateInd;
   int *indPlane     = eesData->gspPlaneInd;
   int ngo           = eesData->nchareGSPProcT;

   GSAtmMsg *msg = new  GSAtmMsg;
   for(int i=0; i<ngo; i++){
     int iadd = gSpacePlaneProxy(indState[i],indPlane[i]).ckLocal()->registrationFlag;
     if(iadd!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Energy : Bad registration cache flag on proc %d %d %d %d\n",
                myid,iadd,indState[i],indPlane[i]);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
     }//endif
     gSpacePlaneProxy(indState[i],indPlane[i]).ckLocal()->acceptEnergy(msg); 
   }//endfor

 }else{

   if(myid==0){
      GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = config.sfpriority-10;
      gSpacePlaneProxy.acceptEnergy(msg);
   }//endif
 
 }//endif

//-------------------------------------------------------------------------
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyStruct GetEnergyStruct() {
    return egroupProxy.ckLocalBranch()->getEnergyStruct();
}
//==========================================================================
