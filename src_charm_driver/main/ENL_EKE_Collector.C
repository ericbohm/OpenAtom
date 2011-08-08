#include "charm++.h"
//#include "util.h"
#include "ENL_EKE_Collector.h"

// this stuff is in a C file because multiple includes of .def will
// give you too many implementations of this tiny thing for linktime joy
extern FILE *openScreenfWrite(const char *dirnameBase, const char *fname, int temper, int bead, bool beadfile);
ENL_EKE_Collector::ENL_EKE_Collector(int _numEnergyInputs)
{
  temperScreenFile=openScreenfWrite("TEMPER_OUT", "screen", thisIndex, 0, false);
  enlIteration=0; // having two of these is silly, but simple
  ekeIteration=0;
  energyExpected=_numEnergyInputs;
  ENL=EKE=0.0;
  countENL=countEKE=0;
}
void ENL_EKE_Collector::printENL()
{   
  
  fprintf(temperScreenFile,"[%d] ENL         = %5.8lf\n", enlIteration++, ENL/(double) energyExpected);   // tell the world
  ENL=0;
  countENL=0;
}
void ENL_EKE_Collector::printEKE()
{
  fprintf(temperScreenFile,"Iter [%d] EKE         = %5.8lf\n", ekeIteration++, EKE/(double) energyExpected);
  countEKE=0;
  EKE=0;
}

void ENL_EKE_Collector::acceptENL(double _enl)
{
  countENL++; 
  ENL+=_enl;
  if(countENL==energyExpected) printENL();
}
void ENL_EKE_Collector::acceptEKE(double _eke)
{
  countEKE++; 
  EKE+=_eke;
  if(countEKE==energyExpected) printEKE();
}


#include "ENL_EKE_Collector.def.h"
