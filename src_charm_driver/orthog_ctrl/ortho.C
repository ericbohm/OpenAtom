//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//   The Ortho object is basically used to accept the rows of the "S" matrix,
//   use the orthoTransform subroutine to obtain the "T" matrix, and send the
//   rows of the "T" matrix back to all the GSpace states.
// 
//   The S_to_T(...) method is invoked when the S-Calculators send their 
//   contributions either through the S_Reducer group or otherwise.
//   Once the conversion from S to T is done, the preparation for the 
//   broadcast back to the S-Calculators is done.
// 
//   This method also takes of performing the load-balancing step.
// 
// 2D version of ortho is permit better parallelization of the diagonal
// and to facilitate pair calculator communication.
//
// Each pair calculator will perform its multiply and reduce the
// result within its "column" AKA quadrant via the Ortho array.
// Resulting in sub matrices on Ortho which could be pasted together
// on ortho[0,0] for the complete matrix.
//============================================================================

#include "charm++.h"
#include "../../include/debug_flags.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"
#include "ortho.h"
#include <unistd.h>
#include "cksparsecontiguousreducer.h"
#include "../../src_mathlib/mathlib.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cporthog.h"

//============================================================================

extern Config config;
extern int nstates;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern  PairCalcID pairCalcID1;
extern  PairCalcID pairCalcID2;
extern CProxy_AtomsGrp atomsGrpProxy;

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
Ortho::Ortho(int s_grain){
    this->n = nstates;
    this->k = s_grain;
    this->k2 = k * k;
    this->A = new double[k2];
    this->B = new double[k2];
    this->C = new double[k2];
    this->tmp_arr = new double[k2];
    num_ready = 0;
    usesAtSync=CmiTrue;
    setMigratable(false);
    got_start = false;
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	ortho = new double[nstates * nstates];
	orthoT = new double[nstates * nstates];
	wallTimeArr = new double[config.maxIter];
      }
    else
      {
	ortho=NULL;
	orthoT=NULL;
	wallTimeArr=NULL;
      }
    numGlobalIter = 0;
    
  }


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::makeSections(int indexSize, int *indexZ){
    int s1=thisIndex.x*k;
    int s2=thisIndex.y*k;
    if(thisIndex.x <= thisIndex.y) //we get the reduction
    {

	pcProxy = initOneRedSect(indexSize, indexZ, 1, &pairCalcID1,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)), s1, s2, 0);
	if(thisIndex.x!=thisIndex.y)
	    thisProxy(thisIndex.y,thisIndex.x).setPCproxy(pcProxy);
    }
    pcLambdaProxy = initOneRedSect(indexSize, indexZ, 1, &pairCalcID2, CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)) , s1, s2, 0);

}

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::do_iteration(void){
    step = 1;
    memset(C, 0, sizeof(double) * k2);
    if(thisIndex.x == thisIndex.y)
      for(int i = 0; i < k; i++)
	C[i * k + i] = 3;
    if(Ortho_use_local_cb)
	matmulProxy1(thisIndex.x, thisIndex.y).ckLocal()->start_mat_mul(A, B, C, -1, 1, Ortho::step_2_cb, (void*) this);
    else
      matmulProxy1(thisIndex.x, thisIndex.y).ckLocal()->start_mat_mul(A, B, C, -1, 1, CkCallback(CkIndex_Ortho::step_2(), thisProxy(thisIndex.x, thisIndex.y)));
  }

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_2(){
    step = 2;
    if(Ortho_use_local_cb)
      matmulProxy2(thisIndex.x, thisIndex.y).ckLocal()->start_mat_mul(B, C, tmp_arr, 0.5, 0, Ortho::step_3_cb, (void*) this);
    else
      matmulProxy2(thisIndex.x, thisIndex.y).ckLocal()->start_mat_mul(B, C, tmp_arr, 0.5, 0, CkCallback(CkIndex_Ortho::step_3(), thisProxy(thisIndex.x, thisIndex.y)));
  }

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_3(){
    step = 3;
    CmiMemcpy(B, A, sizeof(double) * k2);
    if(Ortho_use_local_cb)
      matmulProxy3(thisIndex.x, thisIndex.y).ckLocal()->start_mat_mul(C, B, A, 0.5, 0, Ortho::error_cb, (void*) this);
    else
      matmulProxy3(thisIndex.x, thisIndex.y).ckLocal()->start_mat_mul(C, B, A, 0.5, 0, CkCallback(CkIndex_Ortho::error_step(), thisProxy(thisIndex.x, thisIndex.y)));
  }

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::error_step(){
    step = 4;
    double ret = 0;
    for(int i = 0; i < k2; i++){
      double tmp = B[i] - A[i];
      ret += tmp * tmp;
    }
    double *tmp = B;
    B = tmp_arr;
    tmp_arr = tmp;
    contribute(sizeof(double), &ret, CkReduction::sum_double);
}

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_error(CkReductionMsg *msg) {
    CmiAssert(thisIndex.x == 0 && thisIndex.y == 0);
    //			end_t = CmiWallTimer();
    double error = *((double *) msg->getData());
    delete msg;
    iterations++;
    error = sqrt(error / (n * n));
    //			CkPrintf("%d\t%f\t%g\n", iterations, end_t - start_t, error);
    if(error > INVSQR_TOLERANCE && iterations < INVSQR_MAX_ITER){
      //				start_t = CmiWallTimer();
      thisProxy.do_iteration();
    }
    else{
      thisProxy.collect_results();
    }
  }

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::start_calc(CkReductionMsg *msg){
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	atomsGrpProxy.StartRealspaceForces(); 
	CkPrintf("------------------------------------------------------\n");
	CkPrintf("Iteration %d done\n", numGlobalIter+1);
	CkPrintf("======================================================\n\n");
	CkPrintf("======================================================\n");
      }
    got_start = true;
    int chunksize=k;
    double *S = (double*) msg->getData();
    step = 0;
    iterations = 0;
    if(thisIndex.x < thisIndex.y)
      { //we have a spare to copy to  
	//make a copy
	//	CkPrintf("[%d,%d] sending copy to [%d,%d]\n",thisIndex.x, thisIndex.y,thisIndex.y,thisIndex.x);
	CkReductionMsg *omsg=CkReductionMsg::buildNew(msg->getSize(),msg->getData());
	// transpose it
	double *dest= (double*) omsg->getData();;
	double tmp;
	for(int i = 0; i < chunksize; i++)
	  for(int j = i + 1; j < chunksize; j++){
	    tmp = dest[i * chunksize + j];
	    dest[i * chunksize + j] = dest[j*chunksize + i];
	    dest[j * chunksize + i] = tmp;
	  }
	thisProxy(thisIndex.y,thisIndex.x).start_calc(omsg);

      }
    else if(thisIndex.x > thisIndex.y)
      { //we are a spare 
	//	CkPrintf("[%d,%d] has its copy\n",thisIndex.x, thisIndex.y);
      }

#ifdef _CP_DEBUG_SMAT_
    char fname[80];
    snprintf(fname,80,"smatrix.out_ortho_%d_%d",thisIndex.x,thisIndex.y);
    FILE *outfile = fopen(fname, "w");
    for(int i=0; i<chunksize; i++){
      for(int j=0; j<chunksize; j++){
	fprintf(outfile, "[%d %d] %10.9f \n", i + chunksize*thisIndex.x+1, j+chunksize*thisIndex.y+1, S[i*chunksize+j]);
      }
    }
    fclose(outfile);
#endif
    for(int i = 0; i < k2; i++)
      B[i] = S[i] / 2;
    delete msg;
    memset(A, 0, sizeof(double) * k2);
    step = 0;
    iterations = 0;
    /* see if we have a non-zero part of I or T (=3I) */
    if(thisIndex.x == thisIndex.y)
      for(int i = 0; i < k; i++)
	A[i * k + i] = 1;
    if(num_ready == 3)
      do_iteration();
  }




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_results(void)
  {
#ifdef _CP_DEBUG_TMAT_
    print_results();
#endif
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	wallTimeArr[numGlobalIter] = CkWallTimer();
	if (numGlobalIter+1 == config.maxIter) {
	  ckout << "wall times from within ortho" << endl;
	  int t;
	  for (t = 1; t < config.maxIter; t++)
	    ckout << wallTimeArr[t] - wallTimeArr[t-1] << endl;
	  ckout << CkWallTimer() - wallTimeArr[t-1] << endl;
	  CkPrintf("======================================================\n");
	  CkExit();
	}
	else
	{
	    if(numGlobalIter>0)
		ckout <<"Iteration time : " << wallTimeArr[numGlobalIter] - wallTimeArr[numGlobalIter-1] << endl;
	}
      }
    numGlobalIter++;
    if (numGlobalIter < config.maxIter) 
      {
	if ((config.lbgspace || config.lbpaircalc) &&(numGlobalIter== FIRST_BALANCE_STEP || numGlobalIter % LOAD_BALANCE_STEP == 0)) {
	  CkPrintf("[%d %d] ortho calling atsync with paircalc %d gspace %d iter %d\n",thisIndex.x, thisIndex.y,config.lbpaircalc, config.lbgspace, numGlobalIter);
	  if(thisIndex.x==0 && thisIndex.y==0)
	    {
	      gSpacePlaneProxy.isAtSync(numGlobalIter);
	    }
	  AtSync();
	}
	else
	  resume();
      }
  }




void Ortho::resume() 

//============================================================================
    {//begin routine
//============================================================================
	/*
	 * Consider a set of rows (of thickness config.sGrainSize) starting at s.
	 * Send this set of rows to S_Calculator (s, s). That will use some part
	 * and send the appropriate parts to the correct S_Calculator guys.
	 */
//#if !CONVERSE_VERSION_ELAN || CMK_BLUEGENE_CHARM
      if(thisIndex.y<=thisIndex.x)
	finishPairCalcSection(k2, A, pcProxy);
      //else no one to talk to 
//#endif
	
	// reset values for the next round
        //         memset(ortho, 0, sizeof(double)*nstates*nstates);
	//         memset(orthoT, 0, sizeof(double)*nstates*nstates);

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


void 
Ortho::acceptAllLambda(CkReductionMsg *msg) {
    delete msg;
    CkAbort("do not call acceptAllLambda");

}



void 
Ortho::acceptSectionLambda(CkReductionMsg *msg) {
    
  double *lambda = (double *)msg->getData();
  int lambdaCount = msg->getSize()/sizeof(double);
#ifdef _CP_DEBUG_LMAT_
  char lmstring[80];
  snprintf(lmstring,80,"lmatrix_%d_%d.out",thisIndex.x,thisIndex.y);
  FILE *fp = fopen(lmstring,"w");
  for(int i=0; i<lambdaCount; i++) {
      fprintf(fp,"[%d] %.12g\n",i, lambda[i]);
    }
  fclose(fp);
#endif
  // revise this to do a matmul replacing multiplyforgamma
  if(!scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt)
    {
	CkAbort("the gamma code is broken\n");
      // transform Tlambda = T*lambda: store in lambda
      double *scr = new double [nstates*nstates];
      CmiMemcpy(scr,lambda, nstates*nstates*sizeof(double));
      multiplyForGamma(orthoT, scr, lambda, nstates);
      delete [] scr;
      finishPairCalcSection2(lambdaCount, lambda, orthoT, pcLambdaProxy);
      // finish pair calc
      CkPrintf("[%d,%d] finishing\n",thisIndex.x, thisIndex.y);
    }
  else
    {
     // finish pair calc
	// look up asymmetric backward path make sure we send right part
	finishPairCalcSection(lambdaCount, lambda, pcLambdaProxy);
      //CkPrintf("[%d,%d] finishing asymm\n",thisIndex.x, thisIndex.y);
    }
    
  delete msg;
  //----------------------------------------------------------------------------
}// end routine
//==============================================================================


void Ortho::multiplyForGamma(double *orthoT, double *lambda, double *gamma, int n)
{
  double alpha=1.0;
  double beta=0.0;
  int n_in=n;
  int m_in=n;
  int k_in=n;
  char transform='N';
  char transformT='T';
  DGEMM(&transform, &transform, &n_in, &m_in, &k_in, &alpha, orthoT, &n_in, lambda, &k_in, &beta, gamma, &n_in);      

}
