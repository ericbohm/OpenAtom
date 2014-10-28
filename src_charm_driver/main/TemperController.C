#include "charm++.h"
#include "utility/util.h"
#include "cpaimd.h"
#include "energyGroup.h"
#include "InstanceController.h"
#include "TemperController.h"


extern CProxy_InstanceController         instControllerProxy;

TemperController::TemperController()
{
  /* insert initialization here */  
  reportedIn=0;
  numTempers=config.UberKmax;
  numBeads=config.UberImax;
  CkPrintf("TemperController for Tempers %d and beads %d\n",numTempers,numBeads);
  temperEnergy=new EnergyStruct[numTempers];
  temperatures=new double[numTempers];
  for(int i=0;i<numTempers;i++)   temperatures[i]=.0; //placeholder
}

void TemperController::acceptData(int temper, EnergyStruct &energies)
{
  reportedIn++;
  //temperEnergy[temper]+=energies; // note += not defined on this
  //struct, but that is the notional operation

  /* we should be summing the energies here, as you will get one per
     bead, but which ones we actually care about summing is for Glenn
     to figure out. */
  if(reportedIn == numTempers*numBeads)
    acceptData();

}

void TemperController::acceptData()
{
  reportedIn=0;
  CkPrintf("Hey Glenn, this is where the Temper controller could do its job\n");
  /* insert temperature manipulation code here*/

  /* putting this in a for loop because we expect this to be personalized */
  UberCollection index;
  // average out the 
  for(int i=0; i < numTempers; i++)
  { 
    index.idxU.x=i;
    index.setPO();
    //      CkPrintf("sending temp to %d,%d,%d,%d\n",index.idxU.x, index.idxU.y, index.idxU.z, index.idxU.s);
    instControllerProxy[i].acceptNewTemperature(temperatures[i]);
  }
}

TemperController::~TemperController()
{
  delete [] temperEnergy;
  delete [] temperatures;
}
#include "temperController.def.h"
