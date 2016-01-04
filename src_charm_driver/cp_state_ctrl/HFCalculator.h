#ifndef _HFCALCULATOR_H_
#define _HFCALCULATOR_H_

#include "charm++.h"
#include "debug_flags.h"
#include "ckmulticast.h"
#include "HartreeFock.decl.h"

class HFInputMsg : public CMessage_HFInputMsg {
  public:
    int inputSize;
    int statenum;
    int chunknum;
    double *inputPsi;
};

class HFCalculator : public CBase_HFCalculator {
  public:
    HFCalculator_SDAG_CODE;
    HFCalculator();
    void initializeOuterproduct();
    void aggregateOuterproduct();
    void fillOneState(HFInputMsg *msg);
    int totalStates;
    int mygridx;
    int mygridy;
    int mygridz;
    int lcount;
    int scount;
    double *oneState;
    double **myouterproduct;
    int totalnodes;
    int mynode;
    int myrowstart;
    int mynumrows;
    int oneStateSize;
};
#endif
