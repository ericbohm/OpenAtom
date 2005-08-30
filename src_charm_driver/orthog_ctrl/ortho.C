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
Ortho::Ortho(int s_grain, CProxySection_PairCalculator pcProxy, CProxySection_PairCalculator pcLambdaProxy){
    this->n = nstates;
    this->k = s_grain;
    this->k2 = k * k;
    this->A = new double[k2];
    this->B = new double[k2];
    this->C = new double[k2];
    this->pcProxy = pcProxy;
    this->pcLambdaProxy = pcLambdaProxy;
    this->tmp_arr = new double[k2];
    num_ready = 0;
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
    memcpy(B, A, sizeof(double) * k2);
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
      if(thisIndex.x==0 && thisIndex.y==0){
	FILE *outfile = fopen("smatrix.out_ortho_0_0", "w");
	for(int i=0; i<chunksize; i++){
	  for(int j=0; j<chunksize; j++){
	    fprintf(outfile, "[%d %d] %10.9f \n", i, j, S[i*chunksize+j]);
	  }
	}
	fclose(outfile);
      }
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
// this won't ever be called now
void Ortho::S_to_T(CkReductionMsg *msg) 

//============================================================================
    { //begin routine 
//============================================================================
    CkDataSegHeader r;
    double zero=0;

#ifdef CONVERSE_VERSION_ELAN
    double *data = (double *)msg->getData();
#else
    double *data = (double *)decompressMsg(msg, r, zero);
#endif
    int i,j;
    for (i = 0; i < nstates; i++)
        for (j = 0; j < nstates; j++)
	  if (j >= i){
	    ortho[i * nstates + j] = data[(i*nstates+j)];
	  }
	  else {
	    ortho[i * nstates + j] = data[(j*nstates+i)];
	  }

#ifdef _CP_DEBUG_SMAT_
	FILE *outfile = fopen("smatrix.out_ortho", "w");
	for(int i=0; i<nstates; i++){
	  for(int j=0; j<nstates; j++){
	    fprintf(outfile, "[%d %d] %10.9f \n", i, j, ortho[i*nstates+j]);
	  }
	}
	fclose(outfile);
#endif
	CkPrintf("------------------------------------------------------\n");
	CkPrintf("Iteration %d done\n", numGlobalIter+1);
	CkPrintf("======================================================\n\n");
	CkPrintf("======================================================\n");

    
/*
#if CONVERSE_VERSION_ELAN && ! CMK_BLUEGENE_CHARM
	sReducerProxy[0].prepareBroadcast();
#endif
*/
	
	// convert the "S" matrix to "T"  
   CPORTHOG::CP_orthoTransform(orthoT, ortho , nstates);
   //   CP_orthoTransform(scProxy.ckLocalBranch()->cpcharmParaInfo, orthoT, ortho , nstates * nstates);
	// S and T matrix are symmetric and real
#ifdef _CP_DEBUG_TMAT_
	FILE *outfile2 = fopen("tmatrix.out", "w");
	for(int i=0; i<nstates; i++){
	  for(int j=0; j<nstates; j++){
           fprintf(outfile2, "[%d %d] %10.9f \n",i+1,j+1, orthoT[i*nstates+j]);
	  }
	}
	fclose(outfile2);
#endif
    wallTimeArr[numGlobalIter] = CkWallTimer();
    numGlobalIter++;
    if (numGlobalIter == config.maxIter) {
		ckout << "wall times from within ortho" << endl;
		int t;
		for (t = 1; t < config.maxIter; t++)
			ckout << wallTimeArr[t] - wallTimeArr[t-1] << endl;
		ckout << CkWallTimer() - wallTimeArr[t-1] << endl;
  	        CkPrintf("======================================================\n");
		CkExit();
    }

    
    if ((config.lbgspace || config.lbpaircalc) &&(numGlobalIter== FIRST_BALANCE_STEP || numGlobalIter % LOAD_BALANCE_STEP == 0)) {
      CkPrintf("ortho calling atsync with paircalc %d gspace %d\n",config.lbpaircalc, config.lbgspace);
      //      if(config.lbpaircalc)
      //	{
	  isAtSyncPairCalc(&pairCalcID1);
	  isAtSyncPairCalc(&pairCalcID2);
	  //	}
	  //      if(config.lbgspace)
        gSpacePlaneProxy.isAtSync(numGlobalIter);
    }
    else
        resume();
    
    delete msg;
#ifndef CONVERSE_VERSION_ELAN
    delete [] data;
#endif
//----------------------------------------------------------------------------
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_results(void)
  {
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
	  CkPrintf("ortho calling atsync with paircalc %d gspace %d\n",config.lbpaircalc, config.lbgspace);
	  //      if(config.lbpaircalc)
	  //	{
	  if(thisIndex.x==0 && thisIndex.y==0)
	    {
	      isAtSyncPairCalc(&pairCalcID1);
	      isAtSyncPairCalc(&pairCalcID2);
	      gSpacePlaneProxy.isAtSync(numGlobalIter);
	    }
	  AtSync();
	  //	}
	  //      if(config.lbgspace)
	}
	else
	  resume();
      }

#ifdef _CP_DEBUG_TMAT_
    print_results();
#endif
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
    
  CkDataSegHeader r;
  double zero=0;
#ifdef CONVERSE_VERSION_ELAN
  double *lambda = (double *)msg->getData();
#else
  double *lambda = (double *)decompressMsg(msg, r, zero);
#endif



#ifdef _CP_DEBUG_LMAT_
  FILE *fp = fopen("lmatrix.out","w");
  for(int i=0; i<nstates; i++) {
    for(int j=0; j<nstates; j++) {
      fprintf(fp,"[%d %d] %.12g\n",i+1, j+1, lambda[i*nstates+j]);
    }
  }
  fclose(fp);
#endif

  /* This is done before createPaircalc in cpaimd.C now
    CkArrayIndex2D myindex(thisIndex.x, thisIndex.y);
  CkCallback cb(CkIndex_CP_State_GSpacePlane::acceptLambda(NULL), 
		myindex, thisProxy.ckGetArrayID());
  */
  // move load balancing and iteration control stuff here


  if(!scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt)
    {
      // transform Tlambda = T*lambda: store in lambda
      double *scr = new double [nstates*nstates];
      memcpy(scr,lambda, nstates*nstates*sizeof(double));
      multiplyForGamma(orthoT, scr, lambda, nstates);
      delete [] scr;
      finishPairCalc2(&pairCalcID2, nstates*nstates, lambda, orthoT);
      // finish pair calc
      CkPrintf("[%d,%d] finishing\n",thisIndex.x, thisIndex.y);
     
    }
  else
    {
     // finish pair calc
      finishPairCalc(&pairCalcID2, nstates*nstates, lambda);
      CkPrintf("[%d,%d] finishing asymm\n",thisIndex.x, thisIndex.y);
    }
    
  delete msg;

#ifndef CONVERSE_VERSION_ELAN
  delete [] lambda;
#endif

  //----------------------------------------------------------------------------
}// end routine
//==============================================================================

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
      memcpy(scr,lambda, nstates*nstates*sizeof(double));
      multiplyForGamma(orthoT, scr, lambda, nstates);
      delete [] scr;
      finishPairCalc2(&pairCalcID2, nstates*nstates, lambda, orthoT);
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
