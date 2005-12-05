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

//#define DEBUG_GLENN

//==============================================================================
#include "charm++.h"
#include "groups.h"
#include "util.h"
#include "cpaimd.h"
#include <math.h>
#include "sim_subroutines.h"
#include "CP_State_Plane.h"

//----------------------------------------------------------------------------
#define CHARM_ON
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#include "../../src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"

//----------------------------------------------------------------------------
extern int atom_integrate_done;  // not a readonly global : a group of one element
extern CProxy_EnergyGroup egroupProxy;
extern Config config;
extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_EnergyGroup egroupProxy;	
void IntegrationComplete(void *, void *);

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
                   int cp_wave_opt_, double kT_, Atom* a, AtomNHC *aNHC){
	natm                = n;
        natm_nl             = n_nl;
        len_nhc             = len_nhc_;
        iextended_on        = iextended_on_;
        cp_min_opt          = cp_min_opt_;
        cp_wave_opt         = cp_wave_opt_;
        iteration           = 0;
        kT                  = kT_;
        pot_ewd_rs_loc      = 0.0;
        pot_ewd_rs          = 0.0;
        eKinetic            = 0.0;
        eKineticNhc         = 0.0;
        potNhc              = 0.0;    
	atoms               = new Atom[natm];
	atomsNHC            = new AtomNHC[natm];
	CmiMemcpy(atoms, a, natm * sizeof(Atom)); // atoms has no vectors
        for(int i=0;i<natm;i++){atomsNHC[i].Init(&aNHC[i]);}
	zeroforces();
        if(iextended_on==1 && cp_min_opt==0){
           zeronhc();
           computeFNhc();
	}//endif
        atom_integrate_done=1;
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
AtomsGrp::~AtomsGrp(){
	delete [] atoms;
	delete [] atomsNHC;
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

   AtomsGrp *ag      = atomsGrpProxy.ckLocalBranch();
   Atom *atoms       = ag->atoms;
   int natm          = ag->natm;
   int myid          = CkMyPe();
   double pot_ewd_rs = 0.0;

   atom_integrate_done = 0;
#ifndef DEBUG_GLENN
   if(myid==0){CPRSPACEION::CP_getionforce(natm,atoms,myid,&pot_ewd_rs);}
#endif

#ifdef GJM_DBG_ATMS
   CkPrintf("GJM_DBG: calling contribute atm forces %d\n",myid);
#endif
   atomsGrpProxy[myid].contributeforces(pot_ewd_rs);
   
//==========================================================================
  }//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================

void AtomsGrp::contributeforces(double pot_ewd_rs){

//==========================================================================

  int i,j;
  AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
  int myid     = CkMyPe();
  int natm     = ag->natm;
  double *ftot = (double *)malloc((3*natm+1)*sizeof(double));

  for(i=0,j=0; i<natm; i++,j+=3){
    ftot[j]   = ag->atoms[i].fx;
    ftot[j+1] = ag->atoms[i].fy;
    ftot[j+2] = ag->atoms[i].fz;
  }//endfor
  ftot[3*natm]=pot_ewd_rs;
  ag->pot_ewd_rs = pot_ewd_rs_loc;
#ifdef GJM_DBG_ATMS
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
  AtomsGrp *ag      = atomsGrpProxy.ckLocalBranch();
  EnergyGroup *eg   = egroupProxy.ckLocalBranch();
  int myid          = CkMyPe();
  int natm          = ag->natm;
  int len_nhc       = ag->len_nhc;
  int cp_min_opt    = ag->cp_min_opt;
  int cp_wave_opt   = ag->cp_wave_opt;
  int iextended_on  = ag->iextended_on;
  Atom *atoms       = ag->atoms;
  AtomNHC *atomsNHC = ag->atomsNHC;
  int output_on     = config.stateOutputOn;

//============================================================
// Copy out the reduction of energy and forces

#ifdef GJM_DBG_ATMS
  CkPrintf("GJM_DBG: inside recv forces %d : %d\n",myid,natm);
#endif
  double pot_ewd_rs = ftot[3*natm];
  double fmag = 0.0;
  for(i=0,j=0;i<natm;i++,j+=3){
    atoms[i].fx = ftot[j];
    atoms[i].fy = ftot[j+1];
    atoms[i].fz = ftot[j+2];
    fmag += (ftot[j]*ftot[j]+ftot[j+1]*ftot[j+1]+ftot[j+2]*ftot[j+2]);
#ifdef GJM_DEBUG_ATMS
    if(myid==0){
      CkPrintf("%d : %g %g %g\n",i,atoms[i].fx,atoms[i].fy,atoms[i].fz);
    }//endif
#endif
  }//endfor
  fmag /= (double)(3*natm);
  fmag  = sqrt(fmag);
  delete msg;

#ifdef GJM_DEBUG_ATMS_EXIT
  if(myid==0){CkExit();}
#endif

//============================================================
// Integrate the atoms

#ifdef GJM_DBG_ATMS
  CkPrintf("GJM_DBG: Before atom integrate %d : %d\n",myid,natm);
#endif
  double eKinetic_loc   =0.0;
  double eKineticNhc_loc=0.0;
  double potNhc_loc     =0.0;
  int iwrite_atm        =0;

#ifdef DEBUG_GLENN
  double dt     = (0.25/0.0241888);
  double omega  = 0.01/dt;
  double omega2 = omega*omega;
  double pot_harm = 0.0;
  if(iteration==0){
    for(i=0,j=0;i<natm;i++){
      atoms[i].x = 0.0;
      atoms[i].y = 0.0;
      atoms[i].z = 0.0;
    }
  }
  for(i=0,j=0;i<natm;i++){
    atoms[i].fx = -omega2*atoms[i].x*atoms[i].m;
    atoms[i].fy = -omega2*atoms[i].y*atoms[i].m;;
    atoms[i].fz = -omega2*atoms[i].z*atoms[i].m;;
    pot_harm   += (atoms[i].m*omega2*(atoms[i].x*atoms[i].x+
                                      atoms[i].y*atoms[i].y+
                                      atoms[i].z*atoms[i].z));
  }//endfor
#endif
  ATOMINTEGRATE::ctrl_atom_integrate(iteration,natm,len_nhc,cp_min_opt,
                    cp_wave_opt,iextended_on,atoms,atomsNHC,myid,
                    &eKinetic_loc,&eKineticNhc_loc,&potNhc_loc,&iwrite_atm,
                    output_on);
#ifdef DEBUG_GLENN
  double etot_atm = eKinetic_loc+eKineticNhc_loc+potNhc_loc+pot_harm;
  CkPrintf("iteration %d : tot energy %g on %d\n",iteration,etot_atm,myid);
#endif
  iteration++;

//============================================================
// Tuck things : At present no reduction required for kin,nhc stuff

  ag->pot_ewd_rs  = pot_ewd_rs;
  ag->eKinetic    = eKinetic_loc;
  ag->eKineticNhc = eKineticNhc_loc;
  ag->potNhc      = potNhc_loc;     
  ag->iteration   = iteration;

  eg->estruct.eewald_real     = pot_ewd_rs;  
  eg->estruct.fmag_atm        = fmag;
  eg->estruct.eKinetic_atm    = eKinetic_loc;
  eg->estruct.eKineticNhc_atm = eKineticNhc_loc;  
  eg->estruct.potNhc_atm      = potNhc_loc;  
  eg->estruct.iteration_atm   = iteration;
  eg->iteration_atm           = iteration;

//============================================================
// Get ready for the next iteration : 
//      zero forces, outputAtmEnergy, atomsDone

  zeroforces();
  outputAtmEnergy();

  if(iwrite_atm!=0){
     int i=0;
     CkCallback cb(CkIndex_AtomsGrp::atomsDone(NULL),atomsGrpProxy);
     contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }else{
     atom_integrate_done = 1;
     // since atoms is a group, should use point to point messages
     // sent to gsp's assigned to my proc rather than a single bcast from proc=0
     if(myid==0){ 
        GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
        CkSetQueueing(msg,CK_QUEUEING_IFIFO);
        *(int*)CkPriorityPtr(msg) = config.sfpriority+config.numSfGrps; 
        gSpacePlaneProxy.acceptAtoms(msg);
     }//endif
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
  AtomsGrp *ag       = atomsGrpProxy.ckLocalBranch();
  EnergyGroup *eg    = egroupProxy.ckLocalBranch();
  int cp_min_opt     = ag->cp_min_opt;
  int myid           = CkMyPe();  
  double eKinetic    = eg->estruct.eKinetic_atm;
  double eKineticNhc = eg->estruct.eKineticNhc_atm;
  double fmag        = eg->estruct.fmag_atm;
  double pot_ewd_rs  = eg->estruct.eewald_real;
  double potNhc      = eg->estruct.potNhc_atm;
  int natm           = ag->natm;
  int len_nhc        = ag->len_nhc;
  double free_atm    = 3*((double)natm);
  double free_Nhc    = free_atm*((double)len_nhc);

  if(myid==0){
     CkPrintf("EWALD_REAL  = %5.8lf\n",pot_ewd_rs);
     if(cp_min_opt==0){
        CkPrintf("atm eKin    = %5.8lf\n",eKinetic);
        CkPrintf("atm eKinNhc = %5.8lf\n",eKineticNhc);
        CkPrintf("atm Temp    = %5.8lf\n",(2.0*eKinetic*BOLTZ/free_atm));
        CkPrintf("atm TempNHC = %5.8lf\n",(2.0*eKineticNhc*BOLTZ/free_Nhc));
        CkPrintf("atm potNHC = %5.8lf\n",potNhc);
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
  int myid = CkMyPe();
  atom_integrate_done = 1;
  if(myid==0){
     GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
     CkSetQueueing(msg, CK_QUEUEING_IFIFO);
     *(int*)CkPriorityPtr(msg) = config.sfpriority+config.numSfGrps; 
     gSpacePlaneProxy.acceptAtoms(msg);
  }//endif
}
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
     // since energy is a group, should use point to point messages
     // sent to gsp's assigned to my proc rather than a single bcast from proc=0
      int myid = CkMyPe();
      if(myid==0){
         GSAtmMsg *msg = new (8*sizeof(int)) GSAtmMsg;
         CkSetQueueing(msg, CK_QUEUEING_IFIFO);
         *(int*)CkPriorityPtr(msg) = config.sfpriority+config.numSfGrps;
         gSpacePlaneProxy.acceptEnergy(msg);
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
