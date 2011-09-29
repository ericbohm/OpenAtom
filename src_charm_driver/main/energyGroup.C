#include "energyGroup.h"
#include "eesCache.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "utility/util.h"
#include "CPcharmParaInfoGrp.h"
#include "load_balance/IntMap.h"
#include "charm++.h"


//==========================================================================
//Energy group that can retrieve the energies from
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
    estruct.eKineticNhc_atm = 0;    // NHC kinetic energy
    estruct.potNhc_atm      = 0;    // NHC pot energy
    estruct.fmag_atm        = 0;    // magnitude of atm forces
    estruct.iteration_atm   = 0;
    estruct.potPIMDChain    = 0;

//-------------------------------------------------------------------------
  } //end routine
//==========================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * CP_StateGspacePlane(0,0) calls this to replicate the energies everywhere for
 * consistency and tolerance checking.
 */
void EnergyGroup::updateEnergiesFromGS(EnergyStruct &es, UberCollection sender) {
//==========================================================================

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
      if(config.UberJmax>1)
	CkPrintf("GLENN!! you need to fix estruct.totalElecEnergy for k-points\n");
      estruct.totalElecEnergy  = es.totalElecEnergy;
      estruct.iteration_gsp= es.iteration_gsp;
      iteration_gsp        = es.iteration_gsp;
#ifdef _DEBUG_ESTRUCT_
       CkPrintf("Energies received %lf, %lf, %lf, %lf, %lf\n", 
                 estruct.enl,estruct.eke,estruct.eext,estruct.ehart, 
                 estruct.eewald_recipt,estruct.egga,estruct.eexc,
                 estruct.fictEke,estruct.totalEnergy);
#endif
    int i=0;
    CkCallback cb(CkIndex_EnergyGroup::energyDone(NULL),UegroupProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);

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
   // need to receive one of these for each k-point
   if(++kpointEnergyDoneCount==config.UberJmax)
     {
       kpointEnergyDoneCount=0;
       if(config.UberKmax>1 && config.temperCycle >0 && iteration_atm % config.temperCycle == 0) // its temper time, 
	 { 
	   ktemps=0;
	   // resumeFromTemper will reactivate us later
	   int i=1;
	   CkCallback cb(CkIndex_EnergyGroup::sendToTemper(NULL),0,  UegroupProxy[thisInstance.proxyOffset]);
	   contribute(sizeof(int),&i,CkReduction::sum_int,cb);
	 }
       else
	 {
	   energyDone();
	 }
     }
}
//==========================================================================


void EnergyGroup::sendToTemper(CkReductionMsg *m)
{
  delete m;
  temperControllerProxy[0].acceptData(thisInstance.idxU.x, estruct);
}

void EnergyGroup::resumeFromTemper()
{
 // you will receive 1 per kpoint, only release when all are ready
 ktemps++;
 if(ktemps==config.UberJmax)
   {
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
   for(int kpoint=0; kpoint< config.UberJmax; kpoint++){ //each
							 //k-point
							 //needs to be
							 //handled
     //     CkPrintf("{%d}[%d] EnergyGroup::energyDone.\n ", thisInstance.proxyOffset, CkMyPe());     
     UberCollection thisPoint=thisInstance;
     thisPoint.idxU.y=kpoint; // not at the gamma point
     thisPoint.setPO();
     
     eesCache *eesData = UeesCacheProxy[thisPoint.proxyOffset].ckLocalBranch ();
     int *indState     = eesData->gspStateInd;
     int *indPlane     = eesData->gspPlaneInd;
     int ngo           = eesData->nchareGSPProcT;
     for(int i=0; i<ngo; i++){
       int iadd = UgSpacePlaneProxy[thisPoint.proxyOffset](indState[i],indPlane[i]).ckLocal()->registrationFlag;
       if(iadd!=1){
	 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	 CkPrintf("Energy : Bad registration cache flag on proc %d %d %d %d\n",
		  myid,iadd,indState[i],indPlane[i]);
	 CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
	 CkExit();
       }//endif
       UgSpaceDriverProxy[thisPoint.proxyOffset](indState[i],indPlane[i]).doneComputingEnergy(iteration_atm); 
     }//endfor
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
