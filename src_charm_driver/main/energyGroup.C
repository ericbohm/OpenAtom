#include "energyGroup.h"
#include "eesCache.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "utility/util.h"
#include "CPcharmParaInfoGrp.h"
#include "load_balance/IntMap.h"
#include "charm++.h"
extern CkVec <CProxy_EnergyGroup>          UegroupProxy;
extern CProxy_TemperController temperControllerProxy;
extern CkVec <CProxy_eesCache>             UeesCacheProxy;
extern CkVec <CProxy_CP_State_GSpacePlane> UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>         UgSpaceDriverProxy;
extern CkVec <CProxy_AtomsCompute>         UatomsComputeProxy;
//==========================================================================
//Energy group for each Uber to hold the energies
// the energies are explicitely reduced over Ubers for spin and k-points
// but not for beads or tempers
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyGroup::EnergyGroup (UberCollection _thisInstance) : thisInstance(_thisInstance) {
  iteration_gsp = 0;
  iteration_atm = 0;
  kpointEnergyDoneCount=0;
  //non local
  estruct.enl = 0;
  //local external energy
  estruct.eext = 0;
  estruct.eke = 0;

  //hartree energy
  estruct.ehart = 0;

  //ion-ion
  estruct.eewald_recip = 0;
  estruct.eewald_real  = 0;
  estruct.grimmeVdw    = 0;

  //exchange correlation
  estruct.egga = 0;
  estruct.eexc = 0;
  estruct.fmagPsi = 0;

  //CP Fict KE
  estruct.fictEke = 0;

  // total electronic part
  estruct.totalElecEnergy = 0; // needs ewald_real to be physical
  estruct.iteration_gsp = 0;

  estruct.totalpotPIMDChain = 0.0;

  // atm stuff
  estruct.eKinetic_atm    = 0;    // classical kinetic energy
  estruct.eKineticNhc_atm = 0;    // NHC kinetic energy
  estruct.potNhc_atm      = 0;    // NHC pot energy
  estruct.fmag_atm        = 0;    // magnitude of atm forces
  estruct.iteration_atm   = 0;
  estruct.potPIMDChain    = 0;
  countOfEnergies         = 0;
  //-------------------------------------------------------------------------
} //end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * CP_StateGspacePlane(0,0) from this Uber Instance calls this to
 * replicate the energies everywhere for consistency and tolerance
 * checking.
 */
void EnergyGroup::updateEnergiesFromGS(EnergyStruct &es, UberCollection sender) {
  //==========================================================================
  countOfEnergies++;
  if(config.UberJmax>1 || config.UberMmax >1)
  {// need to sum enl and eke across kpoints and spin
    estruct.enl          += es.enl;
    estruct.eke          += es.eke;
    estruct.fictEke      += es.fictEke;
    estruct.fmagPsi      += es.fmagPsi;
  }
  else
  {
    estruct.enl          = es.enl;
    estruct.eke          = es.eke;
    estruct.fictEke      = es.fictEke;
    estruct.fmagPsi      = es.fmagPsi;
  }
  // these other stuff comes from rho and we only have one of them for
  // all kpoints
  estruct.eext         = es.eext;
  estruct.ehart        = es.ehart;
  estruct.eewald_recip = es.eewald_recip;
  estruct.egga         = es.egga;
  estruct.eexc         = es.eexc;
  // these are gspace things that I don't know what to do with in
  // the kpoint>1 case

  // we can construct this after all k-points have reported
  estruct.totalElecEnergy  = es.totalElecEnergy;
  estruct.iteration_gsp= es.iteration_gsp;
  iteration_gsp        = es.iteration_gsp;
#ifdef _DEBUG_ESTRUCT_
  CkPrintf("Energies received %lf, %lf, %lf, %lf, %lf\n", 
      estruct.enl,estruct.eke,estruct.eext,estruct.ehart, 
      estruct.eewald_recipt,estruct.egga,estruct.eexc,
      estruct.fictEke,estruct.totalEnergy);
#endif
  if(countOfEnergies == config.UberJmax * config.UberMmax)
    {
      int i=0;
      int myBeadIndex    = thisInstance.idxU.x;
      CkCallback cb(CkIndex_EnergyGroup::energyDone(NULL),UegroupProxy[thisInstance.proxyOffset]);
      contribute(sizeof(int),&i,CkReduction::sum_int,cb);
      countOfEnergies=0;
      estruct.totalElecEnergy  = estruct.enl + estruct.eke + estruct.eext + estruct.ehart
	+  estruct.egga+    estruct.eexc;
      CPcharmParaInfo *sim  = CPcharmParaInfo::get(); 
      if(CkMyPe()==0)
	{
	  if(iteration_gsp % sim->nscreen_frq==0 )
	    {
	      if(config.UberJmax>1 || config.UberMmax>1)
		{
		  CkPrintf("[b=%d] Iter [%d] ENL_TOT              = %5.8lf\n", myBeadIndex, iteration_gsp, estruct.enl);
		  CkPrintf("[b=%d] Iter [%d] EKE_TOT              = %5.8lf\n", myBeadIndex, iteration_gsp, estruct.eke);
		}
	      CkPrintf("[b=%d] Iter [%d] ELEC_ENERGY_TOT    = %5.8lf\n", myBeadIndex, iteration_gsp, estruct.totalElecEnergy);
	    }
	  UatomsComputeProxy[thisInstance.proxyOffset](0).energyReady();
	}
      estruct.enl          = 0.0;
      estruct.eke          = 0.0;
      estruct.fictEke      = 0.0;
      estruct.fmagPsi      = 0.0;

    }

  //-------------------------------------------------------------------------
}//end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * call GSpaceDriver->doneComputingEnergy() on all co-located gspace chares
 * which allows the new step to advance with new psi
 */
void EnergyGroup::energyDone(CkReductionMsg *msg) {
  //==========================================================================
  delete msg;
  // we will receive one of these per bead and temper
  if(config.UberKmax>1 && config.temperCycle >0 ) // its temper time, 
    { 
        // resumeFromTemper will reactivate us later
      // when all the energies are done, the 0th element calls sendToTemper
      int i=1;
      CkCallback cb(CkIndex_EnergyGroup::sendToTemper(NULL),0,  UegroupProxy[thisInstance.proxyOffset]);
      contribute(sizeof(int),&i,CkReduction::sum_int,cb);
    }
  else
    {
      energyDone();
    }
}
//==========================================================================


// called on 0th element, send our energies to tempercontroller
// when every instance has done so, it will switch them around
void EnergyGroup::sendToTemper(CkReductionMsg *m)
{
  delete m;
  //  CkPrintf("{%d,%d,%d} energyGroup %d sendToTemper iter %d at freq %d\n",thisInstance.idxU.x, thisInstance.idxU.y, thisInstance.idxU.z, CkMyPe(), iteration_atm, config.temperCycle);
  temperControllerProxy[0].acceptData(thisInstance.idxU.z, iteration_atm, estruct);
}


// triggered by instancecontroller
void EnergyGroup::resumeFromTemper()
{
  // you will receive 1, only release when all are ready
  ktemps++;

  if(ktemps==config.UberJmax*config.UberMmax)
  {
    //    CkPrintf("{%d} energyGroup [%d] resumeFromTemper, my iter %d\n",thisInstance.idxU.z, CkMyPe(), iteration_atm);
    energyDone();
    ktemps=0;
  }
}

//==========================================================================
// Needs to have each proc invoke directly the doneComputingEnergy method of the
// GSpaceDrivers which are mapped to it. Without migration, we have that map
// at startup. With migration, one must write an enroll/dropout routine.
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
void EnergyGroup::energyDone(){
  //==========================================================================
  // Use the cool new data caching system

  int myid          = CkMyPe();
  eesCache *eesData = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  int limit=eesData->gspKptSpinStatePlaneVec.length();

  for(int i=0; i<limit; i++)
    {
      UberCollection thisPoint=thisInstance;
      thisPoint.idxU.y=eesData->gspKptSpinStatePlaneVec[i].getKpoint();
      thisPoint.idxU.s=eesData->gspKptSpinStatePlaneVec[i].getSpin();
      thisPoint.setPO();
      int state=eesData->gspKptSpinStatePlaneVec[i].getState();
      int plane=eesData->gspKptSpinStatePlaneVec[i].getPlane();
      int iadd = UgSpacePlaneProxy[thisPoint.proxyOffset](state,plane).ckLocal()->registrationFlag;
      if(iadd!=1){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("atom : Bad registration cache flag on proc %d %d %d %d\n",
		 myid,iadd,state,plane);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
      }//endif
      //      if(state==0 && plane==0)      CkPrintf("{%d,%d,%d}[%d] energyGroup::energyDone calling gsp(%d,%d).doneComputingEnergy iteration %d thisPoint.idU.y %d thisPoint.idxU.s %d offset %d \n ", thisInstance.idxU.x, thisInstance.idxU.y, thisInstance.idxU.z, CkMyPe(), state,plane, iteration_atm, thisPoint.idxU.y, thisPoint.idxU.s, thisPoint.proxyOffset);     
      UgSpaceDriverProxy[thisPoint.proxyOffset](state,plane).doneComputingEnergy(iteration_atm); 
    }//endfor
}//end routine



/*//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
EnergyStruct GetEnergyStruct() {
return UegroupProxy[thisInstance.proxyOffset].ckLocalBranch()->getEnergyStruct();
}
//==========================================================================
 */
#include "EnergyGroup.def.h"
