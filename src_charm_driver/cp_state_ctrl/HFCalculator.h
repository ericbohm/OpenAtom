#ifndef _HFCALCULATOR_H_
#define _HFCALCULATOR_H_

#include "charm++.h"
#include "debug_flags.h"
#include "ckmulticast.h"
#include "ckcomplex.h"
#include "HartreeFock.decl.h"

class HFInputMsg : public CMessage_HFInputMsg {
  public:
    int inputSize;
    int statenum;
    int chunknum;
    double *inputPsi;
};

class HFTransposeMsg : public CMessage_HFTransposeMsg {
  public:
    int srcNode;
    int destNode;
    int tmsgsize;
    int startColumnNumber;
    int endColumnNumber;
    int totalColumns;
    int numRows;
    int startingRow;
    complex *partialColumn;
};

class HFCalculator : public CBase_HFCalculator {
  public:
    HFCalculator_SDAG_CODE;
    HFCalculator();
    void initializeOuterproduct();
    void aggregateOuterproduct();
    void fillOneState(HFInputMsg *msg);
    void doRowFFTs(complex **mymatrix, int mydir);
    int totalStates;
    int mygridx;
    int mygridy;
    int mygridz;
    int lcount;
    int scount;
    double *oneState;
    complex **myouterproduct;
    complex **transposedouterproduct;
    int totalnodes;
    int mynode;
    int myrowstart;
    int mynumrows;
    int mynumcols;
    int oneStateSize;
    int regularSize;
    int remSize;
    int extendedSize;
};
#endif
