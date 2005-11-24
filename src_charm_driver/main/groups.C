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

#include "charm++.h"
#include "groups.h"
#include "util.h"
#include "cpaimd.h"
#include <math.h>

//----------------------------------------------------------------------------

#include "../../src_piny_physics_v1.0/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cprspaceion.h"

//----------------------------------------------------------------------------

extern int atom_integrate_done;
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
        atom_integrate_done = 1;  // initial conditions are `done'
	atoms               = new Atom[natm];
	atomsNHC            = new AtomNHC[natm];
	CmiMemcpy(atoms, a, natm * sizeof(Atom)); // atoms has no vectors
        for(int i=0;i<natm;i++){atomsNHC[i].Init(&aNHC[i]);}
	zeroforces();
        if(iextended_on==1 && cp_min_opt==0){
           zeronhc();
           computeFNhc();
	}//endif
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

   atom_integrate_done = 0; // flip the global flag
   if(myid==0){CPRSPACEION::CP_getionforce(natm,atoms,myid,&pot_ewd_rs);}

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
  int myid          = CkMyPe();
  int natm          = ag->natm;
  int len_nhc       = ag->len_nhc;
  int cp_min_opt    = ag->cp_min_opt;
  int cp_wave_opt   = ag->cp_wave_opt;
  int iextended_on  = ag->iextended_on;
  Atom *atoms       = ag->atoms;
  AtomNHC *atomsNHC = ag->atomsNHC;

//============================================================
// Copy out the reduction of energy and forces

#ifdef GJM_DBG_ATMS
  CkPrintf("GJM_DBG: inside recv forces %d : %d\n",myid,natm);
#endif
  for(i=0,j=0;i<natm;i++,j+=3){
    atoms[i].fx = ftot[j];
    atoms[i].fy = ftot[j+1];
    atoms[i].fz = ftot[j+2];
#ifdef GJM_DEBUG_ATMS
    if(myid==0){
      CkPrintf("%d : %g %g %g\n",i,atoms[i].fx,atoms[i].fy,atoms[i].fz);
    }//endif
#endif
  }//endfor
  EnergyGroup *eg = egroupProxy.ckLocalBranch();
  eg->estruct.eewald_real=ftot[3*natm];
  delete msg;

#ifdef GJM_DEBUG_ATMS_EXIT
  if(myid==0){CkExit();}
#endif

//============================================================
// Integrate the atoms

#ifdef GJM_DBG_ATMS
  CkPrintf("GJM_DBG: Before atom integrate %d : %d\n",myid,natm);
#endif
  double eKinetic_loc;
  double eKineticNhc_loc;
  double potNhc_loc;
  ATOMINTEGRATE::ctrl_atom_integrate(iteration,natm,len_nhc,cp_min_opt,
                    cp_wave_opt,iextended_on,atoms,atomsNHC,myid,
                    &eKinetic_loc,&eKineticNhc_loc,&potNhc_loc);
  ag->eKinetic    = eKinetic_loc;
  ag->eKineticNhc = eKineticNhc_loc;
  ag->potNhc      = potNhc_loc; // before evolution

//============================================================
// Get ready for the next iteration : 
//    increment counter, flip integrate flag, zero forces.

  ag->iteration++;
  atom_integrate_done = 1;
  zeroforces();

//-------------------------------------------------------------------------
   }//end routine
//==========================================================================


//==========================================================================
//Energy group that can retrieve the energies from
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyGroup::EnergyGroup () {
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

    estruct.fictEke = 0;
    estruct.totalEnergy = 0;
} //end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyStruct GetEnergyStruct() {
    return egroupProxy.ckLocalBranch()->getEnergyStruct();
}
//==========================================================================
