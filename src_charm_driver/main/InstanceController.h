#ifndef _INSTANCECONTROLLER_H_
#define _INSTANCECONTROLLER_H_
//#include "InstanceController.decl.h"
/****************************************************************************
* Handles launch synchronization for each instance.
* 
* Catch all for per instance synchronization of any kind.  
*
****************************************************************************/
class ICCookieMsg : public CkMcastBaseMsg, public CMessage_ICCookieMsg {
 public:
  int junk;
};
class InstanceController: public CBase_InstanceController {
 public:
  double Timer;
  InstanceController();
  ~InstanceController(){}
  InstanceController(CkMigrateMessage *m){}
  void doneInit(CkReductionMsg *msg);
  void initCookie(ICCookieMsg *msg);
  void printEnergyHart(CkReductionMsg *msg);
  void printEnergyEexc(CkReductionMsg *msg);
  void printFictEke(CkReductionMsg *msg);
  void allDoneCPForces(CkReductionMsg *m);
  void allDoneCPForcesAllKPoint(CkReductionMsg *m);
  void printEnergyEke(CkReductionMsg *m);
 private:
  int done_init;
  int numKpointforces;
  CkSectionInfo allKPcookie;
};

#endif
