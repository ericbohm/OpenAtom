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
 * tolerance check is a min reduction on OrthoT.  If it falls out of
 * the tolerance range we have to correct the velocities.  We notify
 * the PC of this via a flag in the finishsection call.  It will then
 * be handled by gspace.  We will keep the tolerance flag and set
 * orthoT to the unit AKA identity matrix for the gamma calculation.
 *
 * We will perform that tolerance check every config.toleranceInterval
 * steps.
 *
 * This complicates interaction with PC a bit since we can no longer
 * just multiply our index by the sGrainSize to determine it.
 *
 */
//============================================================================

#include "charm++.h"
#include "../../include/debug_flags.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"
#include "orthoHelper.h"
#include "ortho.h"
#include <unistd.h>
#include "../../src_mathlib/mathlib.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cporthog.h"
#include "../../src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#define PRINTF CkPrintf

//============================================================================

extern Config config;
extern int nstates;
extern CProxy_TimeKeeper              TimeKeeperProxy;
extern CProxy_main                    mainProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern  PairCalcID pairCalcID1;
extern  PairCalcID pairCalcID2;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_CP_Rho_GHartExt rhoGHartExtProxy;
extern CProxy_OrthoHelper orthoHelperProxy;
extern ComlibInstanceHandle orthoInstance;
extern CProxy_Ortho orthoProxy;
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_error(CkReductionMsg *msg) {
    CmiAssert(thisIndex.x == 0 && thisIndex.y == 0);
    //			end_t = CmiWallTimer();
    double error = *((double *) msg->getData());
    delete msg;
    error = sqrt(error / (nstates * nstates));
    //			CkPrintf("%d\t%f\t%g\n", iterations, end_t - start_t, error);
    if(error > invsqr_tolerance && iterations < invsqr_max_iter){
      //				start_t = CmiWallTimer();
      if(config.useOrthoSection)
	{
	  orthoMtrigger *tmsg= new orthoMtrigger;
	  multiproxy.do_iteration(tmsg);
	}
      else
	thisProxy.do_iteration();
    }
    else{
      if(config.useOrthoSection)
	{
	  orthoMtrigger *tmsg= new orthoMtrigger;
	  multiproxy.collect_results(tmsg);
	}
      else
	thisProxy.collect_results();
    }
  }

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::start_calc(CkReductionMsg *msg){
  int cp_min_opt = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
  int gen_wave   = scProxy.ckLocalBranch()->cpcharmParaInfo->gen_wave;
#ifdef _CP_SUBSTEP_TIMING_
  if(timeKeep>0)
    {
      double ostart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),TimeKeeperProxy);
      contribute(sizeof(double),&ostart,CkReduction::min_double, cb , timeKeep);
    }
#endif

  if(thisIndex.x==0 && thisIndex.y==0)
    {
      if(cp_min_opt==1){
	PRINT_LINE_DASH;
        int iii = numGlobalIter;
        if(gen_wave==0){iii+=1;}
	CkPrintf("Iteration %d done\n",iii);
	PRINT_LINE_STAR; CkPrintf("\n");
	PRINT_LINE_STAR; 
      }else{
	if(numGlobalIter>0){
	  PRINT_LINE_STAR; CkPrintf("\n");
	  PRINT_LINE_STAR;
	}//endif
        if(numGlobalIter<config.maxIter){
  	  CkPrintf("Beginning Iteration %d \n", numGlobalIter);
	}else{
  	  CkPrintf("Completing Iteration %d \n", numGlobalIter-1);
	}//endif
	PRINT_LINE_DASH;
      }//endif
    }
  got_start = true;
  double *S = (double*) msg->getData();
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->getSize()/sizeof(double) ;i++)
    {
      CkAssert(isnan(((double*) msg->getData())[i])==0);
    }
#endif

  step = 0;
  iterations = 0;
  int s1=thisIndex.x * config.orthoGrainSize;
  int s2=thisIndex.y * config.orthoGrainSize;
  int maxpcstateindex=(config.nstates/config.sGrainSize-1)*config.sGrainSize;
  if(config.orthoGrainSize!=config.sGrainSize)
    {
      // do something clever
      s1=s1/config.sGrainSize*config.sGrainSize;
      s2=s2/config.sGrainSize*config.sGrainSize;
      s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
      s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

    }
  if(s1 < s2)   
    {
      //we get the reduction and //we have a spare to copy to  
      //make a copy
      CkReductionMsg *omsg=CkReductionMsg::buildNew(msg->getSize(),msg->getData());
      // transpose it
      /*
      double *dest= (double*) omsg->getData();;
      double tmp;
      for(int i = 0; i < m; i++)
	for(int j = i + 1; j < n; j++){
	  tmp = dest[i * n + j];
	  dest[i * n + j] = dest[j*m + i];
	  dest[j * m + i] = tmp;
	}
      */
      // simple out of place scheme

      double *dest= (double*) omsg->getData();;
      double tmp;
      for(int i = 0; i < m; i++)
	for(int j = 0; j < n; j++){
	  dest[j * m + i] = S[i *n + j];
	}
      thisProxy(thisIndex.y,thisIndex.x).start_calc(omsg);

    }
  else if(s2 < s1)
    { //we get a transposed copy be happy
      
    }
  else if((s1==s2) && (thisIndex.x > thisIndex.y))
    { //we are a spare, got our matrix direct from Scalc 

    }

#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrixDouble("smat",(double *)S, m, n, numGlobalIter, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     

#endif

#ifdef _CP_ORTHO_DEBUG_COMPARE_SMAT_
  if(savedsmat==NULL)
    { // load it
      savedsmat= new double[m*n];
      loadMatrixDouble("smat",(double *)savedsmat, m, n,numGlobalIter, thisIndex.x*config.orthoGrainSize,thisIndex.y*config.orthoGrainSize,0,false);     
    }
  for(int i=0;i<m*n;i++)
    {
      if(fabs(S[i]-savedsmat[i])>0.0001)
	{
	  fprintf(stderr, "O [%d,%d] %d element ortho %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, S[i], savedsmat[i]);
	}

      CkAssert(fabs(S[i]-savedsmat[i])<0.0001);
      CkAssert(fabs(S[i]-savedsmat[i])<0.0001);
    }
#endif
  for(int i = 0; i < m * n; i++){
    B[i] = S[i] / 2.0;
  }
  memset(A, 0, sizeof(double) * m * n);
  step = 0;
  iterations = 0;
  /* see if we have a non-zero part of I or T (=3I) */
  if(thisIndex.x == thisIndex.y){
    for(int i = 0; i < m; i++){
      A[i * m + i] = 1;
    }//endfor
  }//endif
  // do tolerance check on smat, do_iteration will be called by reduction root
  if(cp_min_opt==0 && (numGlobalIter % config.toleranceInterval)==0 && numGlobalIter>1){
    if(thisIndex.x==0 && thisIndex.y==0){
      CkPrintf("doing tolerance check on SMAT \n");
    }//endif
    double max =array_diag_max(m,n,S);
    contribute(sizeof(double),&max, CkReduction::max_double, 
	       CkCallback(CkIndex_Ortho::maxCheck(NULL),CkArrayIndex2D(0,0),
			  thisProxy.ckGetArrayID()));
  }else{
    if(num_ready == 1){do_iteration();}
  }//endif
  delete msg;

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_results(void){
//============================================================================
// Output Timings and debug information then increment iteration counter

    int cp_min_opt  = scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt;
    int gen_wave    = scProxy.ckLocalBranch()->cpcharmParaInfo->gen_wave;

    int iprintout   = config.maxIter;
    if(cp_min_opt==1 && gen_wave==0){iprintout-=1;}

    int itime       = numGlobalIter;
    if(config.maxIter>=30){itime=1; wallTimeArr[0]=wallTimeArr[1];}

    if(thisIndex.x==0 && thisIndex.y==0){
  	wallTimeArr[itime] = CkWallTimer();
	if (numGlobalIter == iprintout && config.maxIter<30) {
	      PRINT_LINE_DASH;
	      CkPrintf("Wall Times from within Ortho\n\n");
	      for (int t = 1; t < iprintout; t++){
		ckout << wallTimeArr[t] - wallTimeArr[t-1] << endl;
	      }//endfor
	  if(itime>0)
	    {
	      CkPrintf("%g\n", wallTimeArr[itime] - wallTimeArr[itime-1]);
	      PRINT_LINE_DASH; CkPrintf("\n");
	    }
        }else{
	  if(numGlobalIter>0){
	    CkPrintf("Iteration time (ORTHO) : %g\n", 
               wallTimeArr[itime] - wallTimeArr[itime-1]);
	    CkPrintf("Ortho S->T iterations: %d\n",iterations); 
	  }//endif
        }//endif
    }//endif

#ifdef _CP_DEBUG_TMAT_
    print_results();
#endif
    numGlobalIter++;
//=======================================================================
// Load balance controller


    if (numGlobalIter <= config.maxIter+1){

      if ((config.lbgspace || config.lbpaircalc ||config.lbdensity) &&
          (numGlobalIter== FIRST_BALANCE_STEP||(numGlobalIter % LOAD_BALANCE_STEP) == 0)){
           CkPrintf("[%d %d] ortho calling atsync with paircalc %d gspace %d iter %d\n",
               thisIndex.x, thisIndex.y,config.lbpaircalc, config.lbgspace, numGlobalIter);
	   AtSync();
   	   if(thisIndex.x==0 && thisIndex.y==0){
             gSpacePlaneProxy.isAtSync(numGlobalIter);
	     rhoRealProxy.isAtSync(numGlobalIter);
	     rhoGProxy.isAtSync(numGlobalIter);
	     rhoGHartExtProxy.isAtSync(numGlobalIter);
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
	//	CkPrintf("O [%d %d] making copy of orthoT m %d n %d\n",thisIndex.x,thisIndex.y,m,n);
	if(orthoT==NULL) //allocate if null
	  { orthoT = new double[m * n];}
	CmiMemcpy(orthoT,A,m*n*sizeof(double));
      }
    int s1=thisIndex.x * config.orthoGrainSize;
    int s2=thisIndex.y * config.orthoGrainSize;
    int maxpcstateindex=(config.nstates/config.sGrainSize-1)*config.sGrainSize;

    if(config.orthoGrainSize!=config.sGrainSize)
      {
	// do something clever
	s1=s1/config.sGrainSize*config.sGrainSize;
	s2=s2/config.sGrainSize*config.sGrainSize;
	s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
	s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
      }
    //    if(thisIndex.y <= thisIndex.x)   //we have the answer scalc wants
    //    if((s2 < s1) || ((s2==s1)&&()))   //we have the answer scalc wants
#ifdef _CP_ORTHO_DUMP_TMAT_
    dumpMatrixDouble("tmat",(double *)A, m, n,numGlobalIter,thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
#endif

#ifdef _CP_ORTHO_DEBUG_COMPARE_TMAT_
  if(savedtmat==NULL)
    { // load it
      savedtmat= new double[m*n];
      loadMatrixDouble("tmat",(double *)savedtmat, m, n, numGlobalIter, thisIndex.x * orthoGrainSize,thisIndex.y*orthoGrainSize,0,false);     
    }
  for(int i=0;i<m*n;i++)
    {
      if(fabs(A[i]-savedtmat[i])>0.0001)
	{
	  fprintf(stderr, "O [%d,%d] %d element ortho %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, A[i], savedtmat[i]);
	}

      CkAssert(fabs(A[i]-savedtmat[i])<0.0001);
      CkAssert(fabs(A[i]-savedtmat[i])<0.0001);
    }
#endif

    if(s1 == s2)   //we have the answer scalc wants
      finishPairCalcSection(m * n, A, &oPairCalcID1, thisIndex.x, thisIndex.y, actionType,  0);
    else if(thisIndex.y < thisIndex.x)   //we have the answer scalc wants
      finishPairCalcSection(m * n, A, &oPairCalcID1, thisIndex.y, thisIndex.x, actionType, 0);
    else if(thisIndex.y > thisIndex.x && config.phantomSym)
      {
	transpose(A,m,n);
	/*	double *dest= (double*) A;
	double tmp;
	for(int i = 0; i < m; i++)
	  for(int j = i + 1; j < n; j++){
	    tmp = dest[i * n + j];
	    dest[i * n + j] = dest[j * m + i];
	    dest[j * m + i] = tmp;
	  }
	*/
	// we have a transposed copy of what scalc wants
	finishPairCalcSection(m * n, A, &oPairCalcID1, thisIndex.x, thisIndex.y, actionType, 0);
#ifdef _CP_ORTHO_DUMP_TMAT_
	dumpMatrixDouble("tmatT",(double *)A, m, n,numGlobalIter,thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     

#endif
      }


#ifdef _CP_SUBSTEP_TIMING_
    if(timeKeep>0)
      {
	double oend=CmiWallTimer();
	CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),TimeKeeperProxy);
	contribute(sizeof(double),&oend,CkReduction::max_double, cb , timeKeep);
      }
#endif
    /* this should be triggered by psi after sym fw path is done */
  if(!scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt && !config.PCCollectTiles)
    {
      sendOrthoTtoAsymm();
    }
//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::sendOrthoTtoAsymm(){
//============================================================================
  int actionType=0; //normal
  sendMatrix(m * n, orthoT, &oPairCalcID2, thisIndex.x, thisIndex.y,0,oPairCalcID2.priority-1 );

}
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
/**
 * this is a very paranoid heavily barriered resumption to avoid
 *  nasty migration race conditions.
 */
void Ortho::lbresume(CkReductionMsg *msg) {
//============================================================================

    delete msg;
    lbcaught++;
    int lambdas=1;
    if(!scProxy.ckLocalBranch()->cpcharmParaInfo->cp_min_opt)
      lambdas=2;
    if(thisIndex.x ==0 && thisIndex.y==0)
      CkPrintf("O [%d %d] caught lb %d",thisIndex.x, thisIndex.y, lbcaught);
    if(lbcaught==lambdas) //gspace is all done lambda reduction reset
	gSpacePlaneProxy.syncpsi();
    if(lbcaught==lambdas+1) //gspace is all done lambda and psi reduction resets
      {
	CkAbort("must fix ortho proxy reset!\n");
	setGredProxy(&oPairCalcID2.proxyAsym, oPairCalcID2.mCastGrpId[0],  CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)),true,CkCallback(CkIndex_Ortho::lbresume(NULL),thisProxy),thisIndex.x, thisIndex.y); 
      }
    if(lbcaught==lambdas+2)
      {
	CkAbort("must fix ortho proxy reset!\n");
	if(thisIndex.x <= thisIndex.y) //lambda is done
	  {
	    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(oPairCalcID1.orthomCastGrpId).ckLocalBranch();               
	    mcastGrp->resetSection(oPairCalcID1.proxySym);
	    setGredProxy(&oPairCalcID1.proxySym, oPairCalcID1.orthomCastGrpId,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)),true,CkCallback(CkIndex_Ortho::lbresume(NULL),thisProxy), thisIndex.x, thisIndex.y);
	    if(thisIndex.x!=thisIndex.y)
	      thisProxy(thisIndex.y,thisIndex.x).setPCproxy(oPairCalcID1.proxySym);	  
	  }
      }
    if(lbcaught==lambdas+3) //everyone is done
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

#ifdef _CP_ORTHO_DUMP_LMAT_
    dumpMatrixDouble("lmat",lambda, m, n, numGlobalIter, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     

#endif
#ifdef _CP_ORTHO_DEBUG_COMPARE_LMAT_
  if(savedlmat==NULL)
    { // load it
      savedlmat= new double[m*n];
      loadMatrixDouble("lmat",(double *)savedlmat, m, n,numGlobalIter, thisIndex.x*config.orthoGrainSize,thisIndex.y*config.orthoGrainSize,0,false);     
    }
  for(int i=0;i<m*n;i++)
    {
      if(fabs(lambda[i]-savedlmat[i])>0.0001)
	{
	  fprintf(stderr, "O [%d,%d] %d element ortho %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, lambda[i], savedlmat[i]);
	}

      CkAssert(fabs(lambda[i]-savedlmat[i])<0.0001);
      CkAssert(fabs(lambda[i]-savedlmat[i])<0.0001);
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
    if(ortho==NULL)
      ortho= new double[m*n];
    CmiMemcpy(ortho,orthoT,m*n*sizeof(double));
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

#ifdef _NAN_CHECK_
      CkAssert(lambdaCount==m*n);
      for(int i=0; i<m*n; i++)
	CkAssert(finite(lambda[i]));
#endif

      // finish pair calc
      finishPairCalcSection(lambdaCount, lambda, &oPairCalcID2, thisIndex.x, thisIndex.y, 0, oPairCalcID2.priority+1);
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

  int s1=thisIndex.x * config.orthoGrainSize;
  int s2=thisIndex.y * config.orthoGrainSize;
  int maxpcstateindex=(config.nstates/config.sGrainSize-1)*config.sGrainSize;
  if(config.orthoGrainSize!=config.sGrainSize)
    {
      // do something clever
      s1=s1/config.sGrainSize*config.sGrainSize;
      s2=s2/config.sGrainSize*config.sGrainSize;
      s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
      s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
    }
  
  // m and n are orthograinsize which must be <=config.sGrainSize
  // thisIndex.x and thisIndex.y range from 0 to nstates/config.orthoGrainSize

  if(s1 <= s2)   //we get the reduction
    {  
      initOneRedSect(indexSize, indexZ, config.numChunksSym, &oPairCalcID1,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)), CkCallback(CkIndex_main::doneInit(NULL),mainProxy), s1, s2, thisIndex.x, thisIndex.y, config.orthoGrainSize, false,false,false);
      if(config.phantomSym)
	{
	  initOneRedSect(indexSize, indexZ, config.numChunksSym, &oPairCalcID1,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)), CkCallback(CkIndex_main::doneInit(NULL),mainProxy), s1, s2, thisIndex.x, thisIndex.y, config.orthoGrainSize, true, true, config.useOrthoDirect);
	}
      else
	{
	  if(s1!=s2)
	    thisProxy(thisIndex.y,thisIndex.x).setPCproxy(oPairCalcID1.proxySym);
	  initOneRedSect(indexSize, indexZ, config.numChunksSym, &oPairCalcID1,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)), CkCallback(CkIndex_main::doneInit(NULL),mainProxy), s1, s2, thisIndex.x, thisIndex.y, config.orthoGrainSize, false, true, config.useOrthoDirect);
	}
    }
  else if(config.phantomSym)
    {  // we are not phantoms
      initOneRedSect(indexSize, indexZ, config.numChunksSym, &oPairCalcID1,  CkCallback(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y)), CkCallback(CkIndex_main::doneInit(NULL),mainProxy),s1, s2, thisIndex.x, thisIndex.y, config.orthoGrainSize, false, true,config.useOrthoDirect);
    }
  if(config.lambdaGrainSize==config.orthoGrainSize)
    { //no point in having a different chare if you won't have more of them
      // in the != case this will happen in the lambda chare
      initOneRedSect(indexSize, indexZ, config.numChunksAsym, &oPairCalcID2, CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)), CkCallback(CkIndex_main::doneInit(NULL),mainProxy) ,s1, s2, thisIndex.x, thisIndex.y, config.orthoGrainSize, false, false, false);
      //everybody sends in lambda
      initOneRedSect(indexSize, indexZ, config.numChunksAsym, &oPairCalcID2, CkCallback(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)) , CkCallback(CkIndex_main::doneInit(NULL),mainProxy), s1, s2, thisIndex.x, thisIndex.y, config.orthoGrainSize, false, true, config.useOrthoDirect);
    }

//----------------------------------------------------------------------------
  }// end routine
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::gamma_done(){
//============================================================================
//  CkPrintf("[%d %d] sending ortho %g %g %g %g gamma %g %g %g
//%g\n",thisIndex.x,
//thisIndex.y,orthoT[0],orthoT[1],orthoT[m*n-2],orthoT[m*n-1],B[0],B[1],B[m*n-2],B[m*n-1]);

#ifdef _CP_ORTHO_DUMP_GMAT_
    dumpMatrixDouble("gmat",(double *)B, m, n, numGlobalIter, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     

#endif

    //CkAbort("HEY I ONLY WANT ONE DYNAMICS STEP");

#ifdef _CP_ORTHO_DEBUG_COMPARE_GMAT_
  if(savedgmat==NULL)
    { // load it
      savedgmat= new double[m*n];
      loadMatrixDouble("gmat",(double *)savedgmat, m, n,numGlobalIter,thisIndex.x*config.orthoGrainSize,thisIndex.y*config.orthoGrainSize,0,false);     
    }
  for(int i=0;i<m*n;i++)
    {
      if(fabs(S[i]-savedgmat[i])>0.0001)
	{
	  fprintf(stderr, "O [%d,%d] %d element ortho %.10g not %.10g\n",thisIndex.x, thisIndex.y,i, S[i], savedgmat[i]);
	}

      CkAssert(fabs(S[i]-savedgmat[i])<0.0001);
      CkAssert(fabs(S[i]-savedgmat[i])<0.0001);
    }
#endif
  if(config.PCCollectTiles)
    {
      // if not streaming we might as well just lockstep these
      finishPairCalcSection2(m * n, B, ortho, &oPairCalcID2,thisIndex.x, thisIndex.y, 0,  oPairCalcID2.priority+1);
    }
  else // orthoT was already sent ahead for processing
    {
      finishPairCalcSection(m * n, B, &oPairCalcID2,thisIndex.x, thisIndex.y, 0,  oPairCalcID2.priority+1);
    }

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
 CLA_Matrix_interface matB3, CLA_Matrix_interface matC3, int timekeep){
//============================================================================
/* do basic initialization */
  parallelStep2=config.useOrthoHelpers;
  
  invsqr_max_iter=config.invsqr_max_iter;
  invsqr_tolerance=config.invsqr_tolerance;
  if(invsqr_tolerance==0)
    invsqr_tolerance=INVSQR_TOLERANCE;
  if(invsqr_max_iter==0)
    invsqr_tolerance=INVSQR_MAX_ITER;
  this->matA1 = matA1; this->matB1 = matB1; this->matC1 = matC1;
  this->matA2 = matA2; this->matB2 = matB2; this->matC2 = matC2;
  this->matA3 = matA3; this->matB3 = matB3; this->matC3 = matC3;
  timeKeep=timekeep;
  int borderOrtho= config.nstates / config.orthoGrainSize-1;
  int remOrtho = config.nstates%config.orthoGrainSize;
  if(thisIndex.x==borderOrtho)
    this->m = m + remOrtho;
  else
    this->m = m;
  if(thisIndex.y==borderOrtho)
    this->n = n + remOrtho;
  else
    this->n = n;
  A = new double[this->m * this->n];
  B = new double[this->m * this->n];
  C = new double[this->m * this->n];
  tmp_arr = new double[this->m * this->n];
  step2done=false;
  step3done=false;
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
  oPairCalcID1=pairCalcID1;
  oPairCalcID2=pairCalcID2;
  if( (thisIndex.x==0 && thisIndex.y==0) && (config.useOrthoSection || config.useOrthoSectionRed))
    {
      int numOrtho=config.nstates/config.orthoGrainSize;
      multiproxy = 
      CProxySection_Ortho::ckNew(thisProxy.ckGetArrayID(),  
				 0, numOrtho-1,1,
				 0, numOrtho-1, 1);

      if(config.useOrthoSectionRed)
	{
	  CProxySection_Ortho rproxy =   multiproxy;
	  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(oPairCalcID1.orthoRedGrpId).ckLocalBranch();               
	  rproxy.ckSectionDelegate(mcastGrp);
	  initCookieMsg *redMsg=new initCookieMsg;
	  rproxy.orthoCookieinit(redMsg);

	}
      if(config.useOrthoSection)
	{
	  if( config.useCommlib && config.useOrthoDirect)
	    {
	      ComlibAssociateProxy(&orthoInstance,multiproxy);	  
	    }
	  else
	    {
	      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(oPairCalcID1.orthomCastGrpId).ckLocalBranch();               
	      multiproxy.ckSectionDelegate(mcastGrp);
	    }
	}
  }
  if(thisIndex.x==0 && thisIndex.y==0 && config.maxIter<30){
    wallTimeArr = new double[config.maxIter+2];
  }else{
    wallTimeArr = new double[30];
  }//endif
  wallTimeArr[0]=0.0;
  wallTimeArr[1]=0.0;

  numGlobalIter = 0;
#ifdef _CP_ORTHO_DEBUG_COMPARE_GMAT_
  savedgmat=NULL;
#endif
#ifdef _CP_ORTHO_DEBUG_COMPARE_SMAT_
  savedsmat=NULL;
#endif
#ifdef _CP_ORTHO_DEBUG_COMPARE_LMAT_
  savedlmat=NULL;
#endif

#ifdef _CP_ORTHO_DEBUG_COMPARE_TMAT_
  savedtmat=NULL;
#endif

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
#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrixDouble("step1:A:",(double *)A, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
    dumpMatrixDouble("step1:B:",(double *)B, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
    dumpMatrixDouble("step1:C:",(double *)C, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
#endif

  matA1.multiply(-1, 1, A, Ortho::step_2_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  CmiNetworkProgress();
  matB1.multiply(-1, 1, B, Ortho::step_2_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  CmiNetworkProgress();
  matC1.multiply(-1, 1, C, Ortho::step_2_cb, (void*) this,
   thisIndex.x, thisIndex.y);
}
//============================================================================


//============================================================================
/* S1 = 0.5 * S3 * S2 (on proxy 2)
 * currently A has T, B has S1, C has S2
 * Multiply tmp_arr = B*C
 * tmp_arr not used in step3, therefore no data dependence
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_2(void){

  if(config.useOrthoHelpers)
    {
      step_2_send();
      step_3();
    }
  else
    {
#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrixDouble("step2:A:",(double *)B, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
    dumpMatrixDouble("step2:B:",(double *)C, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
    //    dumpMatrixDouble("step2:C:",(double *)tmp_arr, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
#endif

      step = 2;
      matA2.multiply(0.5, 0, B, Ortho::step_3_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      matB2.multiply(0.5, 0, C, Ortho::step_3_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      matC2.multiply(0.5, 0, tmp_arr, Ortho::step_3_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
    }
}
//============================================================================


//============================================================================
/** S1 = 0.5 * S3 * S2 (on proxy 2)
 * currently A has T, B has S1, C has S2
 * Multiply tmp_arr = B*C
 * tmp_arr not used in step3, therefore no data dependence
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_2_send(void){
  step = 2;
  // copy our data to the helper
  OrthoHelperMsg *omsg= new (m*n, m*n, 0) OrthoHelperMsg;
  omsg->init(m*n, B,C,0.5, 0.5, 0.5);
  orthoHelperProxy(thisIndex.x,thisIndex.y).recvAB(omsg);
  // will come back in recvStep2
}
//============================================================================

//============================================================================
/** result of 0.5 * S3 * S2 arrives from helper
 *
 */
void Ortho::recvStep2(double *step2result, int size){
  // copy our data into the tmp_arr
  
    CmiMemcpy(tmp_arr, step2result, m * n * sizeof(double));
    step2done=true;
    if(step3done) //end of iteration check
      { 
	tolerance_check();
      }
}

//============================================================================
/* T = 0.5 * S2 * S3 (do S3 = T before) (on proxy 3)
 * currently A has T, B has S1 (old), C has S2, tmp_arr has new S1
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::step_3(){
  step = 3;
  CmiMemcpy(B, A, m * n * sizeof(double));
#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrixDouble("step3:A:",(double *)C, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
    dumpMatrixDouble("step3:B:",(double *)B, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
    dumpMatrixDouble("step3:C:",(double *)A, m, n, iterations, thisIndex.x * config.orthoGrainSize, thisIndex.y * config.orthoGrainSize, 0, false);     
#endif

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
  step2done=false;
  step3done=false;
   
#ifdef _NAN_CHECK_ 
  for(int i = 0; i < m * n; i++)
  {
    CkAssert(finite(A[i]));
    CkAssert(finite(B[i]));
  }
#endif

  double ret = 0;
  for(int i = 0; i < m * n; i++){
    double tmp = B[i] - A[i];
    ret += tmp * tmp;
  }
  double *tmp = B;
  B = tmp_arr;
  tmp_arr = tmp;
  if(config.useOrthoSectionRed)
    {
      CkCallback mycb(CkIndex_Ortho::collect_error(NULL), thisProxy(0, 0));
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(oPairCalcID1.orthomCastGrpId).ckLocalBranch();               
      mcastGrp->contribute(sizeof(double),  &ret, CkReduction::sum_double, orthoCookie, mycb);
    }
  else
    contribute(sizeof(double), &ret, CkReduction::sum_double);
  iterations++;
}


//============================================================================



void Ortho::pup(PUP::er &p){
//    CBase_Ortho::pup(p);
    ArrayElement2D::pup(p);
    p | m;
    p | n;
    p | step;
    p | iterations;
    p | num_ready;
    p | got_start;
    p | multiproxy;
    p | numGlobalIter;
    p | lbcaught;
    p | orthoCookie;
    p | toleranceCheckOrthoT;
    p | oPairCalcID1;    
    p | oPairCalcID2;
    p | step2done;
    p | step3done;
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


void OrthoHelper::sendMatrix()
{
  if(trigger!=NULL)
    delete trigger;
  orthoProxy(thisIndex.x, thisIndex.y).recvStep2(C, m*n);
}


// highly questionable that this is safe
void Ortho::setPCproxy(CProxySection_PairCalculator inproxy)
{
  oPairCalcID1.proxySym=inproxy;
}
