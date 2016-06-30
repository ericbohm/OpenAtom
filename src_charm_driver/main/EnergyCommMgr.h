#ifndef _ENERGYCOMMMGR_H_
#define _ENERGYCOMMMGR_H_

extern CkVec <CProxy_EnergyGroup>                UegroupProxy;
/*****************************************************************************
 * Handles the instance local collective operations on the EnergyGroup

 * Necessary because groups exist on all PEs, though each instance only
 * uses the elements that overlap its UberPE footprint.
 * 
 * Which means that collectives on the group spam all PEs, resulting
 * in unacceptable overheads.
 * 
 * In a better future Charm++ that supports section reductions on
 * group sections, this thing should be removed
 */



class EnergyCommMgr : public CBase_EnergyCommMgr {
  public:
 EnergyCommMgr(UberCollection inst):thisInstance(inst){CkPrintf("{%d}[%d] EnergyCommMgr::EnergyCommMgr on pe %d\n",thisInstance.proxyOffset, thisIndex, CkMyPe());}
  EnergyCommMgr(CkMigrateMessage *m){}
  ~EnergyCommMgr(){}
  void initTemperCookie(ECookieMsg *msg){
    CkGetSectionInfo(temperCookie, msg);
  }
  void resumeFromTemper(ECookieMsg *msg){
    CkAbort("how did we get here?");
    //    EnergyGroup *eg         = UegroupProxy[thisInstance.proxyOffset].ckLocalBranch(); // find me the local copy
    //    eg->resumeFromTemper();
  }
private:
  CkSectionInfo temperCookie;
  UberCollection thisInstance;
};
#include "energyMessages.h"

#endif

