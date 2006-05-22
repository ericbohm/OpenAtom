
/********************************************************************************/
/* Ortho is currently decomposed by sGrainSize.  Which makes sense for
 * the old statewise only decomposition of the paircalculator.  For
 * the new 4d decomposition it may make sense to decompose the
 * paircalculator more along the numPoints axis than the states axis.
 * Which may result in comparitively fewer Ortho objects.  Thereby
 * further motivating the need to chop up Ortho by some unit smaller
 * than grainsize.  So as to avoid Amdaling on the orthonormalization
 * process.
 *
 * If we allow orthograin to be entirely distinct from sGrainSize we
 * have arbitrary overlap situations between orthograins and
 * scalcgrains.  Supporting that requires a rather complicated (and
 * therefore bug prone) contiguousReducer reduction/multicast
 * stitcher/splitter implementation to reduce data from scalc->ortho
 * and multicast data from ortho->scalc.  
 * 
 * We don't want to do that if we don't have to.  If we restrict
 * orthograin to be a factor of sGrainsize then we have no section
 * overlap issues.  Thereby leaving us with ortho sections that need a
 * simple tiling split of the sgrain sections.  Mirrored by a
 * stitching of the submatrix inputs for the backward path.  
 *
 * This can be accomplished manually within the current codebase with
 * some waste in data replication and computation replication to
 * handle the splitting/stiching operations.  
 *
 * A more efficient implementation would adopt the multicast manager
 * group model of building a tree of participants for these
 * operations.  The reduction side from the PC would be broken up into
 * multiple reductions, one for each orthograin within the sgrain.
 * With a separate contribution for each orthograin.  The multicast
 * requires us to stitch together the input matrices into one per
 * sgrain section.  This might be accomplished in two stages, one in
 * which the stitching is done, and a second in which the stitched
 * sgrainsize matrices are multicast.  The alternative is to just
 * multicast the orthograin submatrices where needed and have each
 * scalc do its strided copying stitching.  As stitching is not
 * computationally intensive, this may be the simplest and fastest
 * solution.  The second approach allows you to simply use the
 * reductions and multicasts as mirror uses of the tree.  Where each
 * little ortho can run once it gets its input, while the scalcs would
 * have to assemble their inputs from multiple multicasts.  
 *
 * The assembling approach allows for the later introduction of
 * streaming computations, where the scalc does its multiply using
 * each orthograin submatrix as it arrives.  Considering how messy the
 * backward path is already, introducing that streaming will require
 * significant correctness testing.
 *
 * Implementation details for this require that each ortho object
 * participate in a section which has a section multicast client
 * directed to the sGrainSize PC section.  The converse PC sGrainSize
 * elements will have an array of section cookies, one for each of the
 * subsections for all orthograin elements within the sGrain.  The
 * forward path of the PC will contribute its orthograin tile (via a
 * strided contribute) which will end up at the correct ortho object.
 * 
 * Note: these PC sections must include all 4th dim blocks.  
 *
 * 
********************************************************************************/

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

  void maxCheck(CkReductionMsg *msg);

  void resumeV(CkReductionMsg *msg);

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
    snprintf(outname,80,"tmatrix_t:%d_%d_%d.out",numGlobalIter,thisIndex.x,thisIndex.y);
    FILE *outfile = fopen(outname, "w");
	for(int i=0; i<m; i++){
	  for(int j=0; j<n; j++){
           fprintf(outfile, "%d %d %.12g \n",i+thisIndex.x*n+1,j+thisIndex.y*n+1, A[i*n+j]);
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
    p | multiproxy;
    p | pcProxy;
    p | pcLambdaProxy;
    p | numGlobalIter;
    p | lbcaught;
    p | toleranceCheckOrthoT;
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

/**
 * OrthoT tolerance check util return max value
 */
  inline double array_diag_max(int sizem, int sizen, double *array){
      double absval, max_ret;
      if(thisIndex.x!=thisIndex.y){ //not diagonal
          max_ret=fabs(array[0]);          
	  for(int i=0;i<sizem;i++){
          for(int j=0;j<sizen;j++){
 	    absval=fabs(array[i*sizen+j]);
	    max_ret = (max_ret>absval) ? max_ret : absval;
	  }}//endfor
      }else{ //on diagonal 
          absval=fabs(fabs(array[0]-2.0));
          max_ret = absval;
	  for(int i=0;i<sizem;i++){
	  for(int j=0;j<sizen;j++){
	    absval=fabs(array[i*sizen+j]);
	    if(i!=j){
 	       max_ret = (max_ret>absval) ? max_ret : absval;
	     }else{
	       absval=fabs(absval-2.0);
	       max_ret = (max_ret>absval) ? max_ret : absval;
	     }//endif
	  }}//endfor
      }//endif
      return max_ret;
  }//end routine

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
  CProxySection_Ortho multiproxy;
  CProxySection_PairCalculator pcProxy;
  CProxySection_PairCalculator pcLambdaProxy;
  int num_ready;
  bool got_start;
  int lbcaught;

  bool toleranceCheckOrthoT; //trigger tolerance failure PsiV conditions
  

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
