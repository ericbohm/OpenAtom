#include "energyGroup.h"
#include "eesCache.h"
#include "cp_state_ctrl/CP_State_GSpacePlane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "utility/util.h"
#include "CPcharmParaInfoGrp.h"
#include "load_balance/IntMap.h"
#include "charm++.h"
#include "EnergyCommMgr.h"
#include "ckmulticast.h"
extern CkVec < CkVec <int> > UberPes;

extern CkVec <CProxy_EnergyGroup>          UegroupProxy;
extern CProxy_TemperController temperControllerProxy;
extern CkVec <CProxy_eesCache>             UeesCacheProxy;
extern CkVec <CProxy_CP_State_GSpacePlane> UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>         UgSpaceDriverProxy;
extern CkVec <CProxy_AtomsCompute>         UatomsComputeProxy;
extern CkVec <MapType1> EnergyCommMgrImaptable;
extern CkGroupID mCastGrpId;
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
  ktemps                  = 0;
  //-------------------------------------------------------------------------
} //end routine
//==========================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Creates a spanning tree for all elements of the group in this
 * UberInstance that are within the PEs allocated for this UberInstance.

 * For use in localizing communication by instance.

 */

void EnergyGroup::createSpanningSection() {
  int numpes = UberPes[thisInstance.proxyOffset].length();

  if(EnergyCommMgrImaptable[thisInstance.proxyOffset].get(CkMyPe())>=0){
    int me = CkMyPe();
    for(int i=0;i<numpes;i++)
      {
	if(UberPes[thisInstance.proxyOffset][i]==CkMyPe()){
	  me=i;
	  break;
	}
      }
    //  int BRANCHING_FACTOR= (numPes>3) ? 3: 1;
    //   int BRANCHING_FACTOR= (numPes>1024) ? 5: 3;
    int BRANCHING_FACTOR= 3;

    int numchild = std::min(std::max(0, numpes - BRANCHING_FACTOR * me - 1), BRANCHING_FACTOR);
    int *children=NULL;
    if(numchild)
      children  = new int[numchild];
    for (int i=0; i<numchild; i++)
      {
	children[i] = UberPes[thisInstance.proxyOffset][BRANCHING_FACTOR*me+i+1];
	//	CkPrintf("{%d} %d th child is index %d pe %d \n",thisInstance.proxyOffset, i, BRANCHING_FACTOR*me+i+1, UberPes[thisInstance.proxyOffset][BRANCHING_FACTOR*me+i+1]);
      }
    /* parent of root should be root itself */
    int parent = std::max(0, UberPes[thisInstance.proxyOffset][(me-1)/BRANCHING_FACTOR]);	
    //    CkPrintf("{%d}[%d] me: %d, parent: %d, numchild: %d, numPes in section %d \n", thisInstance.proxyOffset, CkMyPe(), me, parent, numchild, numpes);
    CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    temperCookie = mCastGrp->addToGrpSection(thisgroup, 1, parent, children, numchild); 
    mCastGrp->setReductionClient(temperCookie, new CkCallback(CkIndex_EnergyGroup::sendToTemper(NULL),UberPes[thisInstance.proxyOffset][0],  UegroupProxy[thisInstance.proxyOffset]));
    /* Should be done at the root only, broadcast only supported from root */
    if(CkMyPe() == UberPes[thisInstance.proxyOffset][0]){   //root
      //      CkPrintf("{%d}[%d] section proxy init\n",thisInstance.proxyOffset,CkMyPe());
      secProxy = CProxySection_EnergyGroup(temperCookie, children, numchild); 
      secProxy.ckDelegate(mCastGrp);
      ECookieMsg *msg = new (8*sizeof(int)) ECookieMsg();
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = -1;
      secProxy.initTemperCookie(msg);

    }
  }

}

void EnergyGroup::initTemperCookie(ECookieMsg *msg){
  //  CkPrintf("{%d}[%d] initTemperCookie\n",thisInstance.proxyOffset, CkMyPe());
  CkGetSectionInfo(temperCookie, msg);
  CkCallback cb(CkIndex_EnergyGroup::sectionDone(NULL),UberPes[thisInstance.proxyOffset][0], UegroupProxy[thisInstance.proxyOffset]);
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  int i=1;
  mCastGrp->contribute(sizeof(int),&i,CkReduction::sum_int,temperCookie,cb);
  delete msg;
}

//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * CP_StateGspacePlane(0,0) from this Uber Instance calls this to
 * replicate the energies everywhere for consistency and tolerance
 * checking.
 */
void EnergyGroup::updateEnergiesFromGSSectBcast(EnergyStructMsg *msg)
{

  int numpes = UberPes[thisInstance.proxyOffset].length();
  EnergyStructMsg *emsg= new (8*sizeof(int)) EnergyStructMsg;
  emsg->es=msg->es;
  emsg->sender=msg->sender;
  CkSetQueueing(emsg, CK_QUEUEING_IFIFO);
  *(int*)CkPriorityPtr(emsg) = -1;
  secProxy.updateEnergiesFromGS(emsg);
  delete msg;
}
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Updateenergiesfromgssectbcast forwards this from
 * CP_StateGspacePlane(0,0) for this Uber Instance to replicate the
 * energies everywhere for consistency and tolerance checking.
 */
void EnergyGroup::updateEnergiesFromGS(EnergyStructMsg *emsg)
{
  //==========================================================================
  //  CkPrintf("{%d}[%d] update energies from GS\n",thisInstance.proxyOffset, CkMyPe());
  CkGetSectionInfo(temperCookie, emsg);
  EnergyStruct es=emsg->es;
  UberCollection sender=emsg->sender;
  delete emsg;
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
      // if(EnergyCommMgrImaptable[thisInstance.proxyOffset].get(CkMyPe())>=0){
      //	CkPrintf("{%d}[%d] contributing for energyDone\n",thisInstance.proxyOffset, CkMyPe());
	CkCallback cb(CkIndex_EnergyGroup::energyDone(NULL),UberPes[thisInstance.proxyOffset][0], UegroupProxy[thisInstance.proxyOffset]);
	CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
	mCastGrp->contribute(sizeof(int),&i,CkReduction::sum_int,temperCookie,cb);
	countOfEnergies=0;
	estruct.totalElecEnergy  = estruct.enl + estruct.eke + estruct.eext + estruct.ehart
	  +  estruct.egga+    estruct.eexc;
	CPcharmParaInfo *sim  = CPcharmParaInfo::get(); 
	if(UberPes[thisInstance.proxyOffset][0]==CkMyPe())
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
	//      }
	/*      else
	{
	  CkPrintf("{%d}[%d] make me stop arriving at this point by using section\n",thisInstance.proxyOffset, CkMyPe());
	}
	*/
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
  // CkPrintf("{%d}[%d] arrived in energyDone\n",thisInstance.proxyOffset, CkMyPe());
  if(config.UberKmax>1 && config.temperCycle >0 ) // its temper time, 
    { 
      // CkPrintf("{%d}[%d] energyDone(CkReduction *m) sending to temperControllerProxy[0].acceptData\n",thisInstance.proxyOffset, CkMyPe());
      temperControllerProxy[0].acceptData(thisInstance.idxU.z, iteration_atm, estruct);
        // resumeFromTemper will reactivate us later

    }
  else
    {
      //CkPrintf("{%d}[%d] energyDone(CkReduction *m) called energyDone()\n",thisInstance.proxyOffset, CkMyPe());
      resumeFromTemperSectBcast();
    }
}
//==========================================================================

// called on 0th element, send our energies to tempercontroller
// when every instance has done so, it will switch them around
void EnergyGroup::sendToTemper(CkReductionMsg *m)
{
  delete m;
  temperControllerProxy[0].acceptData(thisInstance.idxU.z, iteration_atm, estruct);
  //CkPrintf("{%d,%d,%d} energyGroup %d sendToTemper iter %d at freq %d\n",thisInstance.idxU.x, thisInstance.idxU.y, thisInstance.idxU.z, CkMyPe(), iteration_atm, config.temperCycle);
 
}

void EnergyGroup::sectionDone(CkReductionMsg *m)
{
  //  CkPrintf("{%d,%d,%d} energyGroup %d sectionDone count %d\n",thisInstance.idxU.x, thisInstance.idxU.y, thisInstance.idxU.z, CkMyPe(), ((int *) m->getData())[0]);
  delete m;
}

void EnergyGroup::resumeFromTemperSectBcast()
{
  //  CkPrintf("{%d}[%d] calling resumeFromTemperSectBcast starting mcast resumeFromTemper UberPes[thisInstance.proxyOffset][0] =%d\n",thisInstance.proxyOffset,CkMyPe(),UberPes[thisInstance.proxyOffset][0]);
  secProxy.resumeFromTemper(new ECookieMsg());
}

// triggered by instancecontroller
void EnergyGroup::resumeFromTemper(ECookieMsg *msg)
{
  CkGetSectionInfo(temperCookie, msg);
  delete msg;
  // you will receive 1, only release when all are ready
  ktemps++;
  //  CkPrintf("[%d] EnergyGroup::resumeFromTemper ktemps %d test %d\n",CkMyPe(), ktemps, config.UberJmax*config.UberMmax);
  if(ktemps==config.UberJmax*config.UberMmax)
  {
    //    CkPrintf("{%d} energyGroup [%d] resumeFromTemper calling energyDone(), my iter %d\n",thisInstance.idxU.z, CkMyPe(), iteration_atm);
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
  //  CkPrintf("{%d}[%d] EnergyGroup::energyDone() \n",thisInstance.proxyOffset, CkMyPe());
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
      //      CkPrintf("{%d,%d,%d}[%d] energyGroup::energyDone calling gsp(%d,%d).doneComputingEnergy iteration %d thisPoint.idU.y %d thisPoint.idxU.s %d offset %d \n ", thisInstance.idxU.x, thisInstance.idxU.y, thisInstance.idxU.z, CkMyPe(), state,plane, iteration_atm, thisPoint.idxU.y, thisPoint.idxU.s, thisPoint.proxyOffset);     
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

#include "energyMessages.def.h"
#include "EnergyGroup.def.h"

