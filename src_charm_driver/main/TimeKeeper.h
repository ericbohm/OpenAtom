//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file TimeKeeper.h
 *
 *  A place to collect substep times.  Each participating chare will register
 *  with the timekeeper on processor 0.  Giving its name and getting the next
 *  available timekeeper id.  Also handles performance counters and projections
 *  start/stop.
 *  
 *  Then at the beginning of its notional step every chare will contribute its
 *  start time to a minimum reduction.  At the end of its notional step each
 *  chare will contribute its end time to a maximum reduction.  The timekeeper
 *  will print out the delta time and name for each substep when it has the min
 *  and max.
 *
 *  HPM will instrument one step using BG/P UPC performance counters.  To use
 *  it, build libhpm.a in src_charm_driver/utilities/ copy it someplace in your
 *  lib search path and add -lhpm to link line Also add
 *  -L/bgsys/drivers/ppcfloor/runtime -lSPI.cna to link line.  DO NOT use
 *  /soft/apps/UPC/lib/libhpm.a it uses MPI and will cause you a lot of grief.
 */

#ifndef TimeKeeper_h
#define TimeKeeper_h

#include <vector>
#include <string>
#include "TopoManager.h"

#ifdef USE_HPM
extern "C" void HPM_Init(int);        
extern "C" void HPM_Start(char *label,int);
extern "C" void HPM_Stop(char *label,int);
extern "C" void HPM_Print(int,int);       
#endif

extern Config config;
extern int TimeKeeperID;
extern vector <string> TimeKeeperNames;

static int keeperRegister(string name)
{
  if(config.useTimeKeeper)
    {
      TimeKeeperNames.push_back(name);
      return(TimeKeeperID++);
    }
  else
    { //disable all timers
      return(-1);
    }
}

class TimeKeeper : public Group
{
 public:
  double *beginTime;
  int HPMCounter;
  int HPMEndCounter;
  int PRJCounter;
  int local_rank;
  TimeKeeper(){
      HPMCounter=-1;
    HPMEndCounter=0;
    PRJCounter=0;
#ifdef USE_HPM
    TopoManager *topoMgr = new TopoManager();
    int x,y,z;
    topoMgr->rankToCoordinates(CkMyPe(),x,y,z,local_rank);
    delete topoMgr;
#endif
  }

  void startTrace()
    {
#ifndef CMK_OPTIMIZE
      if(PRJCounter++==0)
	{
	  traceBegin();
	  //	  CkPrintf("[%d] trace begun \n",CkMyPe());
	}
#endif
    }

  void stopTrace()
    {
#ifndef CMK_OPTIMIZE
      if(--PRJCounter==0)
	{
	  traceEnd();
	  //	  CkPrintf("[%d] trace ended \n",CkMyPe());
	}
#endif
    }

#ifdef USE_HPM
  void initHPM(){
    if(HPMCounter==-1)
      {
	  if(CkMyPe()==0)
	      CkPrintf("[%d] HPM_Init\n",CkMyPe());
	  HPM_Init(local_rank);
	  HPMCounter=0;
      }
  }
  void startHPM(char *string){
      if(HPMCounter++==0)
      {
	  if(CkMyPe()==0)
	      CkPrintf("[%d] HPM_Start\n",CkMyPe());
	  HPM_Start(string, local_rank);
      }
    }
  void stopHPM(char *string)
    {
    if(--HPMCounter==0)
      {
	  if(CkMyPe()==0)
	      CkPrintf("[%d] HPM_Stop\n", CkMyPe());
	  HPM_Stop(string, local_rank);
	  HPMEndCounter=HPMCounter+1;
      }

    }

  void printHPM()
    {
    if(--HPMEndCounter==0)
      {
	  if(CkMyPe()==0)
	      CkPrintf("[%d] HPM_Print\n",CkMyPe());
	  HPM_Print(CkMyPe(), local_rank);
      }

    }
#endif
  // the timekeeper phase stuff is really only expected to work on pe 0
  void init ()
    {
      beginTime=new double[TimeKeeperID+1];
    }
  void collectStart(CkReductionMsg *msg)
    {
      int phase= msg->getUserFlag();
      CkAssert(phase<TimeKeeperID);
      CkAssert(phase>=0);
      beginTime[phase]  = *((double *)msg->getData());
      //      CkPrintf("Phase %s start %.10g\n",TimeKeeperNames[phase].c_str(), beginTime[phase]);
      delete msg;
    }

  // we assume sanity and do our printing when each max arrives
  void collectEnd(CkReductionMsg *msg)
    {
      int phase = msg->getUserFlag();
      CkAssert(phase<TimeKeeperID);
      double endTime  = *((double *)msg->getData());
      delete msg;
      //      CkPrintf("Phase %s start %.10g end %.10g duration %.10g\n",TimeKeeperNames[phase].c_str(),beginTime[phase],endTime, endTime-beginTime[phase]);
      CkPrintf("Phase %s duration %.10g\n",TimeKeeperNames[phase].c_str(), endTime-beginTime[phase]);
    }
};
#endif
