//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file ortho.h
 *
 */
//============================================================================

#ifndef _ortho_h_
#define _ortho_h_

#include "CLA_Matrix.h"
#include "ckPairCalculator.h"
extern Config config;
extern int nstates;
extern  PairCalcID pairCalcID1;
extern  PairCalcID pairCalcID2;
#define INVSQR_TOLERANCE	1.0e-15
#define INVSQR_MAX_ITER		10

class Ortho : public CBase_Ortho{
 public:
  Ortho(){}
  Ortho(CkMigrateMessage *m){}
  ~Ortho(){
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] tmp_arr;
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	delete [] ortho;
	delete [] orthoT;
	delete [] wallTimeArr;
      }
  }
  // catch lambda for later non_minimization use
  void acceptAllLambda(CkReductionMsg *msg); 

  // catch lbresume reduction
  void lbresume(CkReductionMsg *msg); 

  // get our copy of the pcproxy
  void setPCproxy(CProxySection_PairCalculator inproxy){pcProxy=inproxy;}

  // catch lambda for later non_minimization use
  void acceptSectionLambda(CkReductionMsg *msg); 

  //		void S_to_T(CkReductionMsg *); replaced
  void resume();
  void multiplyForGamma(double *orthoT, double *lambda, double *gamma,int n);

  /* start step 1 on proxy 1, the callback will be for step 2
   * S2 = 3 * I - T * S1
   * currently A has T, B has S1, need to construct 3*I in C
   */
  void do_iteration(void);

  void ResumeFromSync() {
    /*
    if(thisIndex.x <= thisIndex.y)
      {

	CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pairCalcID1.mCastGrpId).ckLocalBranch();               
	mcastGrp->resetSection(pcProxy);
	setGredProxy(&pcProxy, pairCalcID1.mCastGrpId,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)),true,CkCallback(CkIndex_Ortho::lbresume(NULL),thisProxy));
	if(thisIndex.x!=thisIndex.y)
	  thisProxy(thisIndex.y,thisIndex.x).setPCproxy(pcProxy);	  

      }

      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pairCalcID2.mCastGrpId).ckLocalBranch();               
//      mcastGrp->resetSection(pcLambdaProxy);

//      resume(); // should be handled by ortho reduction
*/
  }
  void makeSections(int indexSize, int *indexZ);

  void print_results(void){
    char outname[80];
    snprintf(outname,80,"tmatrix%d_%d.out",thisIndex.x,thisIndex.y);
    FILE *outfile = fopen(outname, "w");
	for(int i=0; i<m; i++){
	  for(int j=0; j<n; j++){
           fprintf(outfile, "%d %d %.10g \n",i+thisIndex.x*n+1,j+thisIndex.y*n+1, A[i*n+j]);
	  }
	}
	fclose(outfile);
  }

  void collect_results(void);

  virtual void pup(PUP::er &p){
//    CBase_Ortho::pup(p);
      ArrayElement2D::pup(p);
    p | m;
    p | n;
    p | step;
    p | iterations;
    p | num_ready;
    p | got_start;
    p | pcProxy;
    p | pcLambdaProxy;
    p | numGlobalIter;
    p | lbcaught;
    if(p.isUnpacking() && thisIndex.x==0 &&thisIndex.y==0)
      { 
	ortho = new double[nstates * nstates];
	orthoT = new double[nstates * nstates];
	wallTimeArr = new double[config.maxIter];
      }
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	p(ortho,nstates*nstates);
	p(orthoT,nstates*nstates);
	p(wallTimeArr,config.maxIter);
      }
    if(p.isUnpacking()){
      A = new double[m * n];
      B = new double[m * n];
      C = new double[m * n];
      tmp_arr = new double[m * n];
    }
    p(A, m * n);
    p(B, m * n);
    p(C, m * n);
    p(tmp_arr, m * n);
  }

  void collect_error(CkReductionMsg *msg);

  void start_calc(CkReductionMsg *msg);

  /* called from each CLA_Matrix array (3 per multiplication, 3 mults) */
  void all_ready(){
    num_ready++;
    if(num_ready == 9)
      thisProxy.ready();
  }

  void ready(){ // startup initialization synchronization
    num_ready = 1;
    if(got_start)
      do_iteration();
  }

  Ortho(int m, int n, CLA_Matrix_interface matA1,
   CLA_Matrix_interface matB1, CLA_Matrix_interface matC1,
   CLA_Matrix_interface matA2, CLA_Matrix_interface matB2,
   CLA_Matrix_interface matC2, CLA_Matrix_interface matA3,
   CLA_Matrix_interface matB3, CLA_Matrix_interface matC3);

  void tolerance_check(void);
  void step_2();
  void step_3();
  void gamma_done();

  static inline void step_2_cb(void *obj){
    ((Ortho*) obj)->step_2();
  }

  static inline void step_3_cb(void *obj){
    ((Ortho*) obj)->step_3();
  }

  static inline void tol_cb(void *obj){
    ((Ortho*) obj)->tolerance_check();
  }

  static inline void gamma_done_cb(void *obj){
    ((Ortho*) obj)->gamma_done();
  }

 private:
  double *orthoT; // only used on [0,0]
  double *ortho; //only used on [0,0]
  double *wallTimeArr;//only used on [0,0]
  int numGlobalIter; // global leanCP iterations

  // used in each element
  int iterations; //local inv_sq iterations
  CProxySection_PairCalculator pcProxy;
  CProxySection_PairCalculator pcLambdaProxy;
  int num_ready;
  bool got_start;
  int lbcaught;

  double *A, *B, *C, *tmp_arr;
  int step;
  /* Note, for now m and n are always equal. When we move to chunks not all
   * having the same grain size, these will not be the same and there will
   * probably be some debugging to do.
   */
  int m, n;
  CLA_Matrix_interface matA1, matB1, matC1, matA2, matB2, matC2, matA3,
   matB3, matC3;

};

/**
 * provide procnum mapping for Ortho
 */
class OrthoMap : public CkArrayMap {
  public:
    OrthoMap(int NN):N(NN){}
    virtual int procNum(int arrayHdl, const CkArrayIndex &idx){
      CkArrayIndex2D idx2d = *(CkArrayIndex2D *) &idx;
      return (N * idx2d.index[0] + idx2d.index[1]) % CkNumPes();
    }
  private:
    int N;
};

#endif // #ifndef _ortho_h_
