/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file TimeKeeper.h
 *
 *  A place to collect substep times.  Each participating chare will
 *  register with the timekeeper on processor 0.  Giving its name and
 *  getting the next available timekeeper id.
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
int keeperRegister(string name)
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

class TimeKeeper : public Chare
{
 public:
  double *beginTime;
  TimeKeeper(){
  }
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
