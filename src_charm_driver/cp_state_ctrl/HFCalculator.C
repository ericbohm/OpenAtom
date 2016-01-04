#include "HFCalculator.h"
#include "main/cpaimd.h"
#include "charm++.h"

#include <iostream>
#include <fstream>
#include <cmath>

//============================================================================
extern CkVec <CProxy_CP_State_RealSpacePlane> UrealSpacePlaneProxy;
extern CkGroupID mCastGrpId;
extern Config config;
extern CPcharmParaInfo simReadOnly;
extern bool HartreeFockOn;

HFCalculator::HFCalculator() {
  CkPrintf("Memory usage on node: <%d> : <%.2lf> MB before construction\n",
    CkMyNode(), CmiMemoryUsage()/(1024.0*1024.0));
  mynode = CkMyNode();
  totalnodes = CkNumNodes();
  CkPrintf("HFCalculator created on: <node: %d>\n", mynode);
  totalStates = simReadOnly.nstates;
  mynumrows = totalStates / totalnodes;
  myrowstart = mynode * mynumrows;
  if ((mynode + 1) == totalnodes) {
    int remstates = totalStates % totalnodes;
    mynumrows += remstates;
  }
  mygridx = simReadOnly.sizeX;
  mygridy = simReadOnly.sizeY;
  mygridz = simReadOnly.sizeZ;
  oneStateSize = mygridx * mygridy * mygridz;
  oneState = new double[oneStateSize];
  myouterproduct = new double*[mynumrows];
  for(int i = 0 ; i < mynumrows ; i++) {
    myouterproduct[i] = new double[oneStateSize];
  }
#ifdef _HFCALCULATOR_VERBOSE_
  CkPrintf("totalnodes: <%d>, myrowstart: <%d>, mynumrows: <%d>\n",
    totalnodes, myrowstart, mynumrows);
  CkPrintf("Memory usage on node: <%d> : <%.2lf> MB after construction\n",
    CkMyNode(), CmiMemoryUsage()/(1024.0*1024.0));
  CkPrintf("HFCalculator: totalstates, gridxyz, oneStateSize (in MB): <%d,%d,%d,%d,%.2f>\n",
    totalStates, mygridx, mygridy, mygridz, (double) oneStateSize * sizeof(double)/(1024.0*1024.0));
#endif
}

void HFCalculator::fillOneState(HFInputMsg *msg) {
  int chunksize = mygridx * mygridy;
  CkAssert(chunksize == msg->inputSize);
  int offset = chunksize * (msg->chunknum);
  CmiMemcpy(&oneState[offset], msg->inputPsi, chunksize * sizeof(double));
  delete msg;
}

void HFCalculator::initializeOuterproduct() {
  for(int i = 0 ; i < mynumrows ; i++) {
    for(int j = 0 ; j < oneStateSize ; j++) {
      myouterproduct[i][j] = 0.0;
    }
  }
}

void HFCalculator::aggregateOuterproduct() {
  for(int i = 0 ; i < mynumrows ; i++) {
    double firstTerm = oneState[i];
    for(int j = 0 ; j < oneStateSize ; j++) {
      myouterproduct[i][j] += firstTerm * oneState[j];
    }
  }
}

#include "HartreeFock.def.h"
