/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/
#ifndef TimeKeeper_h
#define TimeKeeper_h
#ifdef CMK_BLUEGENEP
void HPM_Init(void);
void HPM_Start(char *);
void HPM_Stop(char *);
void HPM_Print(void);
#endif
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file TimeKeeper.h
 *
 *  A place to collect substep times.  Each participating chare will
 *  register with the timekeeper on processor 0.  Giving its name and
 *  getting the next available timekeeper id.
 *
 *  Also handles performance counters and projections start/stop.
 *  
 *  Then at the beginning of its notional step every chare will
 *  contribute its start time to a minimum reduction.
 *
 *  At the end of its notional step each chare will contribute its end
 *  time to a maximum reduction.
 *
 *  The timekeeper will print out the delta time and name for each substep
 *  when it has the min and max.
 */

#include <vector>
#include <string>
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
  TimeKeeper(){
    HPMCounter=0;
    HPMEndCounter=0;
    PRJCounter=0;
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

#ifdef CMK_BLUEGENEP
  void initHPM(){
    if(HPMCounter==0)
      {
	HPM_Init();
      }
  }
  void startHPM(const char *string){
    if(HPMCounter++==0)
      {
	HPM_Start(string);
      }
    }
  void stopHPM(const char *string)
    {
    if(--HPMCounter==0)
      {
	HPM_Stop(string);
	HPMEndCounter=HPMCounter+1;
      }

    }

  void printHPM()
    {
    if(--HPMEndCounter==0)
      {
	HPM_Print();
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
