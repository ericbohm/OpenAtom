#ifndef _ortho_h_
#define _ortho_h_

#include "matmul.h"
#include "ckPairCalculator.h"
extern Config config;
extern CProxy_matmul matmulProxy1;
extern CProxy_matmul matmulProxy2;
extern CProxy_matmul matmulProxy3;
extern int Ortho_UE_step2;
extern int Ortho_UE_step3;
extern int Ortho_UE_error;
extern bool Ortho_use_local_cb;
#define INVSQR_TOLERANCE	1.0e-15
#define INVSQR_MAX_ITER		10

class Ortho : public CBase_Ortho{
 public:
  Ortho(){}
  Ortho(CkMigrateMessage *m){}
  Ortho(int s_grain);
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

  // get our copy of the pcproxy
  void setPCproxy(CProxySection_PairCalculator inproxy){pcProxy=inproxy;}

  // catch lambda for later non_minimization use
  void acceptSectionLambda(CkReductionMsg *msg); 

  // redundant
  void S_to_T(CkReductionMsg *msg); 

  /* start step 1 on proxy 1, the callback will be for step 2
   * S2 = 3 * I - T * S1
   * currently A has T, B has S1, need to construct 3*I in C
   */

  //		void S_to_T(CkReductionMsg *); replaced
  void resume();
  void multiplyForGamma(double *orthoT, double *lambda, double *gamma,int n);

  void do_iteration(void);
  /* S1 = 0.5 * S2 * S3 (on proxy 2)
   * currently A has T, B has S1, C has S2
   */
  void makeSections(int indexSize, int *indexZ);
  void step_2();
  void step_2_go(); //for local callback
  /* T = 0.5 * S2 * S3 (S3 = T before) (on proxy 3)
   * currently A has T, B has S1 (old), C has S2, tmp_arr has new S1
   */
  void step_3();
  void step_3_go(); //for local callback
  /* calculate error and reset pointers (for step 1)
   * current: T -> A, S1 -> tmp_arr, S2 -> C, S3 -> B
   */
  void error_step();

  void print_results(void){
    char outname[80];
    snprintf(outname,80,"tmatrix%d_%d.out",thisIndex.x,thisIndex.y);
    FILE *outfile = fopen(outname, "w");
	for(int i=0; i<k; i++){
	  for(int j=0; j<k; j++){
           fprintf(outfile, "%d %d %.10g \n",i+thisIndex.x*k+1,j+thisIndex.y*k+1, A[i*k+j]);
	  }
	}
	fclose(outfile);

  }

  void collect_results(void);

  virtual void pup(PUP::er &p){
    CBase_Ortho::pup(p);
    p | n;
    p | k;
    p | k2;
    p | step;
    p | iterations;
    p | num_ready;
    p | got_start;
    p | pcProxy;
    p | pcLambdaProxy;
    p | numGlobalIter;
    if(p.isUnpacking() && thisIndex.x==0 &&thisIndex.y==0)
      { 
	ortho = new double[n * n];
	orthoT = new double[n * n];
	wallTimeArr = new double[config.maxIter];
      }
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	p(ortho,n*n);
	p(orthoT,n*n);
	p(wallTimeArr,config.maxIter);
      }
    if(p.isUnpacking()){
      A = new double[k2];
      B = new double[k2];
      C = new double[k2];
      tmp_arr = new double[k2];
      if(step == 1)
        if(Ortho_use_local_cb)
	  matmulProxy1(thisIndex.x, thisIndex.y).ckLocal()->reset(C, Ortho::step_2_cb, (void*) this);
        else
	  matmulProxy1(thisIndex.x, thisIndex.y).ckLocal()->reset(C);
      else if(step == 2)
        if(Ortho_use_local_cb)
	  matmulProxy2(thisIndex.x, thisIndex.y).ckLocal()->reset(tmp_arr, Ortho::step_3_cb, (void*) this);
        else
	  matmulProxy2(thisIndex.x, thisIndex.y).ckLocal()->reset(tmp_arr);
      else if(step == 3)
        if(Ortho_use_local_cb)
	  matmulProxy3(thisIndex.x, thisIndex.y).ckLocal()->reset(A, Ortho::step_3_cb, (void*) this);
        else
	  matmulProxy3(thisIndex.x, thisIndex.y).ckLocal()->reset(A);
    }
    p(A, k2);
    p(B, k2);
    p(C, k2);
    p(tmp_arr, k2);
  }

  void collect_error(CkReductionMsg *msg);

  void start_calc(CkReductionMsg *msg);

  void ready(){ // startup initialization synchronization
    num_ready++;
    if(got_start)
      do_iteration();
  }

  /* below functions used when local callback is being used */
  static inline void step_2_cb(void *obj){
	double start_t = CmiWallTimer();
	((Ortho*) obj)->step_2();
	traceUserBracketEvent(Ortho_UE_step2, start_t, CmiWallTimer());
  }
  static inline void step_3_cb(void *obj){
	double start_t = CmiWallTimer();
	((Ortho*) obj)->step_3();
	traceUserBracketEvent(Ortho_UE_step3, start_t, CmiWallTimer());
  }
  static inline void error_cb(void *obj){
	double start_t = CmiWallTimer();
	((Ortho*) obj)->error_step();
	traceUserBracketEvent(Ortho_UE_error, start_t, CmiWallTimer());
  }

 private:
  double *orthoT; // only used on [0,0]
  double *ortho; //only used on [0,0]
  double *wallTimeArr;//only used on [0,0]
  int numGlobalIter; // global leanCP iterations
  // used in each element
  double *A, *B, *C, *tmp_arr;
  int n, k, k2;
  int step;
  int iterations; //local inv_sq iterations
  CProxySection_PairCalculator pcProxy;
  CProxySection_PairCalculator pcLambdaProxy;
  int num_ready;
  bool got_start;
};



#endif // #ifndef _ortho_h_
