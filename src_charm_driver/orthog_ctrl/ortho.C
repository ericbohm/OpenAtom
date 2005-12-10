//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file ortho.C
 * The Ortho object is basically used to accept the rows of the "S" matrix,
 * use the orthoTransform subroutine to obtain the "T" matrix, and send the
 * rows of the "T" matrix back to all the GSpace states.
 * 
 * The S_to_T(...) method is invoked when the S-Calculators send their 
 * contributions either through the S_Reducer group or otherwise.
 * Once the conversion from S to T is done, the preparation for the 
 * broadcast back to the S-Calculators is done.
 * 
 * This method also takes of performing the load-balancing step.
 * 
 * 2D version of ortho is permit better parallelization of the diagonal
 * and to facilitate pair calculator communication.
 *
 * Each pair calculator will perform its multiply and reduce the
 * result within its "column" AKA quadrant via the Ortho array.
 * Resulting in sub matrices on Ortho which could be pasted together
 * on ortho[0,0] for the complete matrix.
 *
 * For the dynamics (non minization) case we have an additional
 * multiplication of gamma= lambda*orthoT.  And we pass both gamma and
 * orthoT to the PC instead.
 *
 * For dynamics we also have a tolerance check for orthoT.  The
 * tolerance check is a min and max reduction on OrthoT.  If it falls
 * out of the tolerance range we have to correct the velocities.  We
 * notify the PC of this via a flag in the finishsection call.  It
 * will then be handled by gspace.  We will keep the tolerance flag
 * and set orthoT to the unit AKA identity matrix for the gamma
 * calculation.
 *
 * We will perform that tolerance check every config.toleranceInterval
 * steps.
 *
 */
//============================================================================

#include "charm++.h"
#include "../../include/debug_flags.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
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

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_error(CkReductionMsg *msg) {
    CmiAssert(thisIndex.x == 0 && thisIndex.y == 0);
    //			end_t = CmiWallTimer();
    double error = *((double *) msg->getData());
    delete msg;
    iterations++;
    error = sqrt(error / (nstates * nstates));
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
//	atomsGrpProxy.StartRealspaceForces(); 
        if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==1){
  	  CkPrintf("------------------------------------------------------\n");
	  CkPrintf("Iteration %d done\n", numGlobalIter+1);
	  CkPrintf("======================================================\n\n");
	  CkPrintf("======================================================\n");
	}else{
          if(numGlobalIter>0){
 	   CkPrintf("======================================================\n\n");
	   CkPrintf("======================================================\n");
	  }//endif
	  CkPrintf("Beginning Iteration %d \n", numGlobalIter);
  	  CkPrintf("------------------------------------------------------\n");
	}//endif
      }
    got_start = true;
    int chunksize = m;
    double *S = (double*) msg->getData();
    step = 0;
    iterations = 0;
    if(thisIndex.x < thisIndex.y)
      { //we have a spare to copy to  
	//make a copy
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
    for(int i = 0; i < m * n; i++)
      B[i] = S[i] / 2;
    memset(A, 0, sizeof(double) * m * n);
    step = 0;
    iterations = 0;
    /* see if we have a non-zero part of I or T (=3I) */
    if(thisIndex.x == thisIndex.y)
      for(int i = 0; i < m; i++)
	A[i * m + i] = 1;
    if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt==0 && numGlobalIter % config.toleranceInterval==0 && numGlobalIter!=0) 
      { // do tolerance check on smat, do_iteration will be called by reduction root
 	CkPrintf("doing tolerance check on SMAT \n");
	double max =array_diag_max(m,n, S);
	contribute(sizeof(double),&max, CkReduction::max_double, CkCallback(CkIndex_Ortho::maxCheck(NULL),CkArrayIndex2D(0,0),thisProxy.ckGetArrayID()));
      }
    else if(num_ready == 1)
      do_iteration();
    delete msg;
//============================================================================
  }
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_results(void){
//============================================================================
// Output Timings and debug information then increment iteration counter

    int cp_min_opt  = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
    int iprintout;
    if(cp_min_opt==1){iprintout=config.maxIter-1;}else{iprintout=config.maxIter;}
    int itime       = numGlobalIter;
    if(config.maxIter>=30){itime=1; wallTimeArr[0]=wallTimeArr[1];}

#ifdef _CP_DEBUG_TMAT_
    print_results();
#endif

    if(thisIndex.x==0 && thisIndex.y==0){
  	wallTimeArr[itime] = CkWallTimer();
	if (numGlobalIter == iprintout && config.maxIter<30) {
          CkPrintf("----------------------------\n");
	  ckout << "wall times from within ortho" << endl;
	  for (int t = 1; t < iprintout; t++){
	    ckout << wallTimeArr[t] - wallTimeArr[t-1] << endl;
	  }//endfor
	  ckout << wallTimeArr[itime] - wallTimeArr[itime-1] << endl;
	  CkPrintf("------------------------------\n");
        }else{
	  if(numGlobalIter>0){
	    ckout <<"Iteration time (ORTHO) : " << 
               wallTimeArr[itime] - wallTimeArr[itime-1] << endl;
	  }//endif
        }//endif
    }//endif

    numGlobalIter++;

//=======================================================================
// Load balance controller

    if (numGlobalIter <= config.maxIter+1){

      if ((config.lbgspace || config.lbpaircalc) &&
          (numGlobalIter== FIRST_BALANCE_STEP||(numGlobalIter % LOAD_BALANCE_STEP) == 0)){
           CkPrintf("[%d %d] ortho calling atsync with paircalc %d gspace %d iter %d\n",
               thisIndex.x, thisIndex.y,config.lbpaircalc, config.lbgspace, numGlobalIter);
   	   if(thisIndex.x==0 && thisIndex.y==0){
             gSpacePlaneProxy.isAtSync(numGlobalIter);
	   }//endif
       }else{
 	  resume();
       }//endif

    }//endif : we are still going

//----------------------------------------------------------------------
  }//end routine
//=======================================================================


//============================================================================
/**
 * Resume execution by finishing the backward path of the Psi calculator
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::resume(){
//============================================================================
      int actionType=0; //normal
      if(toleranceCheckOrthoT)
	{
	  actionType=1;  //keepOrtho
	}
      if(scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt!=1){
	// copy orthoT for use in gamma computation
	if(orthoT==NULL) //allocate if null
	  { orthoT = new double[m * n];}
	memcpy(orthoT,A,m*n*sizeof(double));
      }
      if(thisIndex.y<=thisIndex.x)
	finishPairCalcSection(m * n, A, pcProxy, actionType);

//----------------------------------------------------------------------------
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::maxCheck(CkReductionMsg *msg){
//============================================================================

  double tolMax=fabs(((double *) msg->getData())[0]);
  delete msg;

  CkPrintf("SMAT tol    = %g\n", tolMax);
  if(tolMax < scProxy.ckLocalBranch()->cpcharmParaInfo->tol_norb){
      toleranceCheckOrthoT=false;
      thisProxy.do_iteration();
  }else{
      // smat is outside of the tolerance range  need new PsiV
      toleranceCheckOrthoT=true;
      CkPrintf("recalculating PsiV due to tolerance failure \n");
      gSpacePlaneProxy.requirePsiV();  //gspace will trigger our resume
  }//endif

//============================================================================
   }//end routine
//============================================================================

/**
 * Resume execution with the vpsi tolerance correction flag on
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::resumeV(CkReductionMsg *msg){ // gspace tolerance check entry
  delete msg;
  toleranceCheckOrthoT=true;
  do_iteration();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::acceptAllLambda(CkReductionMsg *msg) {
    delete msg;
    CkAbort("do not call acceptAllLambda");
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::lbresume(CkReductionMsg *msg) {
//============================================================================
    delete msg;
    lbcaught++;
    if(thisIndex.x ==0 && thisIndex.y==0)
      CkPrintf("O [%d %d] caught lb %d",thisIndex.x, thisIndex.y, lbcaught);
    if(lbcaught==1) //gspace is all done lambda reduction reset
	gSpacePlaneProxy.syncpsi();
    if(lbcaught==2) //gspace is all done lambda and psi reduction resets
      setGredProxy(&pcLambdaProxy, pairCalcID2.mCastGrpId,  CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)),true,CkCallback(CkIndex_Ortho::lbresume(NULL),thisProxy));      
    if(lbcaught==3)
      {
	if(thisIndex.x <= thisIndex.y) //lambda is done
	  {
	    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pairCalcID1.mCastGrpId).ckLocalBranch();               
	    mcastGrp->resetSection(pcProxy);
	    setGredProxy(&pcProxy, pairCalcID1.mCastGrpId,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)),true,CkCallback(CkIndex_Ortho::lbresume(NULL),thisProxy));
	    if(thisIndex.x!=thisIndex.y)
	      thisProxy(thisIndex.y,thisIndex.x).setPCproxy(pcProxy);	  
	  }
      }
    if(lbcaught==4) //everyone is done
    {
	CkPrintf("O [%d %d] resumes\n",thisIndex.x,thisIndex.y);
	resume();
	lbcaught=0;
    }
//============================================================================
  }//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::acceptSectionLambda(CkReductionMsg *msg) {
//============================================================================

  double *lambda = (double *)msg->getData();
  int lambdaCount = msg->getSize()/sizeof(double);

#ifdef _CP_DEBUG_LMAT_
  char lmstring[80];
  snprintf(lmstring,80,"lmatrix_%d_%d.out",thisIndex.x,thisIndex.y);
  FILE *fp = fopen(lmstring,"w");
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fprintf(fp, "[%d %d] %10.9f \n", i + n*thisIndex.x+1, j+n*thisIndex.y+1, lambda[i*n+j]);
    }
  }
#endif

  // revise this to do a matmul replacing multiplyforgamma
  if(!scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt){

    if(toleranceCheckOrthoT)
      {// replace orthoT with the identity matrix
	if(thisIndex.x!=thisIndex.y)
	  { //not on diagonal
	    bzero(orthoT,m*n*sizeof(double));
	  }
	else 
	  { // on diagonal
	    for(int i=0;i<m;i++)
	      for(int j=0;j<n;j++)
		{
		  if(i!=j)
		    orthoT[i*n+j]=0.0;	      
		  else
		    orthoT[i*n+j]=1.0;	      
		}
	  }
	toleranceCheckOrthoT=false;
      }
    matA1.multiply(1, 0, orthoT, Ortho::gamma_done_cb, (void*) this,
		   thisIndex.x, thisIndex.y);
    matB1.multiply(1, 0, lambda, Ortho::gamma_done_cb, (void*) this,
		   thisIndex.x, thisIndex.y);
    matC1.multiply(1, 0, B, Ortho::gamma_done_cb, (void*) this,
		   thisIndex.x, thisIndex.y);
    //completed gamma will call finishPairCalcSection2
  }
  else
    {
      // finish pair calc
      finishPairCalcSection(lambdaCount, lambda, pcLambdaProxy,0);
#ifdef _CP_DEBUG_ORTHO_VERBOSE_
      if(thisIndex.x==0 && thisIndex.y==0)
	CkPrintf("[%d,%d] finishing asymm\n",thisIndex.x, thisIndex.y);
#endif
    }
  delete msg;

//----------------------------------------------------------------------------
  }// end routine
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::makeSections(int indexSize, int *indexZ){
    int s1=thisIndex.x*m;
    int s2=thisIndex.y*n;
    if(thisIndex.x <= thisIndex.y) //we get the reduction
    {

	pcProxy = initOneRedSect(indexSize, indexZ, 1, &pairCalcID1,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)), s1, s2, 0);
	if(thisIndex.x!=thisIndex.y)
	    thisProxy(thisIndex.y,thisIndex.x).setPCproxy(pcProxy);
    }
    pcLambdaProxy = initOneRedSect(indexSize, indexZ, 1, &pairCalcID2, CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)) , s1, s2, 0);

//----------------------------------------------------------------------------
  }// end routine
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::gamma_done(){
//============================================================================
  //  CkPrintf("[%d %d] sending ortho %g %g %g %g gamma %g %g %g %g\n",thisIndex.x, thisIndex.y,orthoT[0],orthoT[1],orthoT[m*n-2],orthoT[m*n-1],B[0],B[1],B[m*n-2],B[m*n-1]);
#ifdef _CP_DEBUG_GMAT_
    char fname[80];
    snprintf(fname,80,"gamma.out_ortho_%d_%d",thisIndex.x,thisIndex.y);
    FILE *outfile = fopen(fname, "w");
    for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
	fprintf(outfile, "[%d %d] %10.9f \n", i + n*thisIndex.x+1, j+n*thisIndex.y+1, B[i*n+j]);
      }
    }
    fclose(outfile);
#endif

  finishPairCalcSection2(m * n, B, orthoT, pcLambdaProxy,0);

//----------------------------------------------------------------------------
  }// end routine
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::multiplyForGamma(double *orthoT, double *lambda, double *gamma, int n){
//============================================================================

  double alpha=1.0;
  double beta=0.0;
  int n_in=n;
  int m_in=n;
  int k_in=n;
  char transform='N';
  char transformT='T';
  DGEMM(&transform, &transform, &n_in, &m_in, &k_in, &alpha, orthoT, &n_in, lambda, 
        &k_in, &beta, gamma, &n_in);      

//-----------------------------------------------------------------------------
  }// end routine
//==============================================================================


//============================================================================
//New functions necessary for using CLA_Matrix
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
Ortho::Ortho(int m, int n, CLA_Matrix_interface matA1,
 CLA_Matrix_interface matB1, CLA_Matrix_interface matC1,
 CLA_Matrix_interface matA2, CLA_Matrix_interface matB2,
 CLA_Matrix_interface matC2, CLA_Matrix_interface matA3,
 CLA_Matrix_interface matB3, CLA_Matrix_interface matC3){
//============================================================================
/* do basic initialization */

  this->m = m;
  this->n = n;
  this->matA1 = matA1; this->matB1 = matB1; this->matC1 = matC1;
  this->matA2 = matA2; this->matB2 = matB2; this->matC2 = matC2;
  this->matA3 = matA3; this->matB3 = matB3; this->matC3 = matC3;
  A = new double[m * n];
  B = new double[m * n];
  C = new double[m * n];
  tmp_arr = new double[m * n];
  step = 0;
  lbcaught=0;
  num_ready = 0;
  usesAtSync=CmiTrue;
  setMigratable(false);
  got_start = false;
  toleranceCheckOrthoT=false;
  ortho=NULL;
  orthoT=NULL;
  wallTimeArr=NULL;
  if(thisIndex.x==0 && thisIndex.y==0 && config.maxIter<30){
    wallTimeArr = new double[config.maxIter+2];
  }else{
    wallTimeArr = new double[30];
    wallTimeArr[0]=0.0;
    wallTimeArr[1]=0.0;
  }//endif
  numGlobalIter = 0;


//============================================================================
   }//end routine
//============================================================================


//============================================================================
/* start step 1 on proxy 1, the callback will be for step 2
 * S2 = 3 * I - T * S1
 * currently A has T, B has S1, need to construct 3*I in C
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::do_iteration(void){
  step = 1;
  memset(C, 0, m * n * sizeof(double));
  if(thisIndex.x == thisIndex.y){
    for(int i = 0; i < n; i++)
      C[i * n + i] = 3;
  }
  matA1.multiply(-1, 1, A, Ortho::step_2_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matB1.multiply(-1, 1, B, Ortho::step_2_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matC1.multiply(-1, 1, C, Ortho::step_2_cb, (void*) this,
   thisIndex.x, thisIndex.y);
}
//============================================================================


//============================================================================
/* S1 = 0.5 * S3 * S2 (on proxy 2)
 * currently A has T, B has S1, C has S2
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_2(void){
  step = 2;
  matA2.multiply(0.5, 0, B, Ortho::step_3_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matB2.multiply(0.5, 0, C, Ortho::step_3_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matC2.multiply(0.5, 0, tmp_arr, Ortho::step_3_cb, (void*) this,
   thisIndex.x, thisIndex.y);
}
//============================================================================


//============================================================================
/* T = 0.5 * S2 * S3 (do S3 = T before) (on proxy 3)
 * currently A has T, B has S1 (old), C has S2, tmp_arr has new S1
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_3(){
  step = 3;
  memcpy(B, A, m * n * sizeof(double));
  matA3.multiply(0.5, 0, C, Ortho::tol_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matB3.multiply(0.5, 0, B, Ortho::tol_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matC3.multiply(0.5, 0, A, Ortho::tol_cb, (void*) this,
   thisIndex.x, thisIndex.y);
}
//============================================================================


//============================================================================
/* calculate error and reset pointers (for step 1)
 * current: T -> A, S1 -> tmp_arr, S2 -> C, S3 -> B
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::tolerance_check(){
  step = 4;
  double ret = 0;
  for(int i = 0; i < m * n; i++){
    double tmp = B[i] - A[i];
    ret += tmp * tmp;
  }
  double *tmp = B;
  B = tmp_arr;
  tmp_arr = tmp;
  contribute(sizeof(double), &ret, CkReduction::sum_double);
}
//============================================================================


