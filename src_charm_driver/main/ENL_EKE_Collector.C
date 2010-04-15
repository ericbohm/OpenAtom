#include "charm++.h"
#include "ENL_EKE_Collector.h"

// this stuff is in a C file because multiple includes of .def will
// give you too many implementations of this tiny thing for linktime joy

ENL_EKE_Collector::ENL_EKE_Collector(int _numEnergyInputs)
{
  energyExpected=_numEnergyInputs;
  ENL=EKE=0.0;
  countENL=countEKE=0;
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
