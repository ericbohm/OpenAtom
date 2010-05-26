#ifndef _INSTANCECONTROLLER_H_
#define _INSTANCECONTROLLER_H_
//#include "InstanceController.decl.h"
/****************************************************************************
* Handles launch synchronization for each instance.
* 
* Catch all for per instance synchronization of any kind.  
*
****************************************************************************/
class InstanceController: public CBase_InstanceController {
 public:
  double Timer;
  InstanceController() {done_init=0;Timer=CmiWallTimer();}
  ~InstanceController(){}
  InstanceController(CkMigrateMessage *m){}
  void doneInit(CkReductionMsg *msg);
  void printEnergyHart(CkReductionMsg *msg);
  void printEnergyEexc(CkReductionMsg *msg);
  void printFictEke(CkReductionMsg *msg);
  void allDoneCPForces(CkReductionMsg *m);
  void printEnergyEke(CkReductionMsg *m);
  void cleanExit(CkReductionMsg *m);
  void cleanExitAll(CkReductionMsg *m);
 private:
  int done_init;
  
};

#endif
