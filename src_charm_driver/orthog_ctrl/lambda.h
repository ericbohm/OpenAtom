
/********************************************************************************
 * Lambda created by Eric Bohm 2006/7/23
 *
 *
 * 
********************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file lambda.h
 *
 */
//============================================================================

#ifndef _lambda_h_
#define _lambda_h_

#include "CLA_Matrix.h"
#include "ckPairCalculator.h"
extern Config config;
extern int nstates;
extern  PairCalcID pairCalcID1;
extern  PairCalcID pairCalcID2;

class Lambda : public CBase_Lambda{
 public:
  Lambda(){}
  Lambda(CkMigrateMessage *m){}
  Lambda(int lgrain, int ograin, int sgrain) : lambdaGrainSize(lgrain), orthoGrainSize(ograin), sGrainSize(sgrain) 
    {
      lPairCalcID2=pairCalcID2;
    }
  
  ~Lambda(){
  }
  // catch lambda for later non_minimization use
  void acceptAllLambda(CkReductionMsg *msg); 

  // catch lbresume reduction
  void lbresume(CkReductionMsg *msg); 

  // catch lambda for later non_minimization use
  void acceptSectionLambda(CkReductionMsg *msg); 

  void makeSections(int indexSize, int *indexZ);

  virtual void pup(PUP::er &p){
//    CBase_Lambda::pup(p);
    ArrayElement2D::pup(p);
    p | lambdaGrainSize;
    p | orthoGrainSize;
    p | sGrainSize;
    p | lPairCalcID2;
  }

 private:
  int lambdaGrainSize;
  int orthoGrainSize;
  int sGrainSize;
  PairCalcID lPairCalcID2;
};

/**
 * provide procnum mapping for Ortho
 */
class LambdaMap : public CkArrayMap {
  public:
  LambdaMap(int NN,int _nLambda, int _stride):N(NN), nLambda(_nLambda), stride(_stride)
      {
	offset=0;
	if(nLambda<CkNumPes())
	  offset=1;  //skip proc 0
      }
    virtual int procNum(int arrayHdl, const CkArrayIndex &iIndex){
      int *index=(int *) iIndex.data();
      return (stride*(N * index[0] + index[1]) + offset) % CkNumPes();
    }
  private:
    int N;
    int nLambda;
    int offset;
    int stride;
};

#endif // #ifndef _lambda_h_
