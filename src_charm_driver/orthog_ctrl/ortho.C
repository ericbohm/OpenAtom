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

#include "ortho.h"
#include "orthoHelper.h"
//#include "gSpaceDriver.decl.h"
#include "timeKeeper.decl.h"

#include "charm++.h"
#include "utility/matrix2file.h"
#include "utility/util.h"
//#include "main/groups.h"

//#include "fft_slab_ctrl/fftCacheSlab.h"
#include <unistd.h>

#ifdef TRACE_MEMORY
    #include <trace-projections.h>
#endif

#include "src_mathlib/mathlib.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cporthog.h"
#include "src_piny_physics_v1.0/include/class_defs/piny_constants.h"
#define PRINTF CkPrintf
//============================================================================

extern Config config;
extern CProxy_TimeKeeper              TimeKeeperProxy;
extern ComlibInstanceHandle orthoInstance;
//============================================================================


Ortho::~Ortho()
{
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




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::collect_error(CkReductionMsg *msg) {
    #ifdef VERBOSE_ORTHO
        CkPrintf("[%d,%d] Ortho::collect_error \n", thisIndex.x, thisIndex.y);
    #endif
    CmiAssert(thisIndex.x == 0 && thisIndex.y == 0);
    //			end_t = CmiWallTimer();
    double error = *((double *) msg->getData());
    error = sqrt(error / (cfg.numStates * cfg.numStates));
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
      if(iterations>=config.invsqr_max_iter)
	{
	  CkPrintf("Ortho reached max_iter %d with residual of %.10g, which is greater than invsqr_tolerance %.10g!\nEither increase invsqr_max_iter or invsqr_tolerance.\n  If this is not the first step in a GenWave try lowering your timestep\n",config.invsqr_max_iter, error, config.invsqr_tolerance);
	  CkAbort("invsqr_tolerance not met within invsqr_max_iter\n");
	}


      if(config.useOrthoSection)
	{
	  orthoMtrigger *tmsg= new orthoMtrigger;
	  multiproxy.collect_results(tmsg);
	}
      else
	thisProxy.collect_results();
    }
    delete msg;
    
  }

//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::start_calc(CkReductionMsg *msg){
    #ifdef VERBOSE_ORTHO
        CkPrintf("[%d,%d] Ortho::start_calc \n", thisIndex.x, thisIndex.y);
    #endif
#ifdef _CP_SUBSTEP_TIMING_
  if(timeKeep>0)
    {
      double ostart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&ostart,CkReduction::min_double, cb , timeKeep);
    }
#endif

  if(thisIndex.x==0 && thisIndex.y==0)
    {
      if(!cfg.isDynamics){
	PRINT_LINE_DASH;
        int iii = numGlobalIter;
        if(!cfg.isGenWave){iii+=1;}
	CkPrintf("{%d} Iteration %d done\n",cfg.instanceIndex,iii);
	PRINT_LINE_STAR; CkPrintf("\n");
	PRINT_LINE_STAR;
      }else{
	if(numGlobalIter>0){
	  PRINT_LINE_STAR; CkPrintf("\n");
	  PRINT_LINE_STAR;
	}//endif
        if(numGlobalIter<config.maxIter){
  	  CkPrintf("{%d} Beginning Iteration %d \n",cfg.instanceIndex, numGlobalIter);
	}else{
  	  CkPrintf("{%d} Completing Iteration %d \n", cfg.instanceIndex,numGlobalIter-1);
	}//endif
	PRINT_LINE_DASH;
      }//endif
    }
  got_start = true;
  internalType *S = (internalType*) msg->getData();
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->getSize()/sizeof(internalType) ;i++)
      CkAssert( isfinite(S[i]) );
#endif

  step = 0;
  iterations = 0;
  CkIndex2D pc = symmSectionMgr.computePCStateIndices();
  if(pc.x < pc.y)   
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

      internalType *dest= (internalType*) omsg->getData();;
      internalType tmp;
      for(int i = 0; i < m; i++)
          for(int j = 0; j < n; j++)
              dest[j * m + i] = S[i *n + j];

      thisProxy(thisIndex.y,thisIndex.x).start_calc(omsg);

    }
  else if(pc.y < pc.x)
    { //we get a transposed copy be happy
      
    }
  else if((pc.x==pc.y) && (thisIndex.x > thisIndex.y))
    { //we are a spare, got our matrix direct from Scalc 

    }

#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrix("smat",(double *)S, m, n, numGlobalIter, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);

#endif

#ifdef _CP_ORTHO_DEBUG_COMPARE_SMAT_
  if(savedsmat==NULL)
    { // load it
      savedsmat= new double[m*n];
      loadMatrix("smat",(double *)savedsmat, m, n,numGlobalIter, thisIndex.x*cfg.grainSize,thisIndex.y*cfg.grainSize,0,false);
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
  memset(A, 0, sizeof(internalType) * m * n);
  step = 0;
  iterations = 0;
  /* see if we have a non-zero part of I or T (=3I) */
  if(thisIndex.x == thisIndex.y){
    for(int i = 0; i < m; i++){
      A[i * m + i] = 1;
    }//endfor
  }//endif
  // do tolerance check on smat, do_iteration will be called by reduction root
  if(cfg.isDynamics && (numGlobalIter % config.toleranceInterval)==0 && numGlobalIter>1){
    if(thisIndex.x==0 && thisIndex.y==0){
      CkPrintf("{%d} doing tolerance check on SMAT \n",cfg.instanceIndex);
    }//endif
    double max =array_diag_max(m,n,S);
    contribute(sizeof(double),&max, CkReduction::max_double, 
	       CkCallback(CkIndex_Ortho::maxCheck(NULL),CkArrayIndex2D(0,0),
			  thisProxy.ckGetArrayID()));
  }else{
    if(num_ready == 10){do_iteration();}
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

    int iprintout   = config.maxIter;
    if(!cfg.isDynamics && !cfg.isGenWave){iprintout-=1;}

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
	    CkPrintf("{%d}Iteration time (ORTHO) : %g\n", 
		     cfg.instanceIndex,
               wallTimeArr[itime] - wallTimeArr[itime-1]);
	    CkPrintf("{%d} Ortho S->T iterations: %d\n",cfg.instanceIndex,iterations);
	  }//endif
        }//endif
    }//endif

#ifdef _CP_DEBUG_TMAT_
    print_results();
#endif
#ifdef TRACE_MEMORY
    traceMemoryUsage();
#endif
    numGlobalIter++;
//=======================================================================
// Load balance controller


    if (numGlobalIter <= config.maxIter+1)
 	  resume();
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
      if(cfg.isDynamics){
	// copy orthoT for use in gamma computation
	//	CkPrintf("O [%d %d] making copy of orthoT m %d n %d\n",thisIndex.x,thisIndex.y,m,n);
	if(orthoT==NULL) //allocate if null
	  { orthoT = new internalType[m * n];}
	CmiMemcpy(orthoT,A,m*n*sizeof(internalType));
      }
    CkIndex2D pc = symmSectionMgr.computePCStateIndices();
    //    if(thisIndex.y <= thisIndex.x)   //we have the answer scalc wants
    //    if((pc.y < pc.x) || ((pc.y==pc.x)&&()))   //we have the answer scalc wants
#ifdef _CP_ORTHO_DUMP_TMAT_
    dumpMatrix("tmat",(double *)A, m, n,numGlobalIter,thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
#endif

#ifdef _CP_ORTHO_DEBUG_COMPARE_TMAT_
  if(savedtmat==NULL)
    { // load it
      savedtmat= new double[m*n];
      loadMatrix("tmat",(double *)savedtmat, m, n, numGlobalIter, thisIndex.x * cfg.grainSize, thisIndex.y*cfg.grainSize,0,false);
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

    /**
     * orthos talking to non-phantoms in the symm pc instance will have to perform an extra transpose of their data.
     * Hence, they do this only when they have to. We avoid this when phantoms are turned off by rigging the pc
     * sections to make the mirror orthos deliver the data to our non-phantoms. 
     * Refer PCSectionManager::setupArraySection for more info.
     */
    if(pc.x == pc.y)   //we have the answer scalc wants
      symmSectionMgr.sendResults(m*n, A, 0, thisIndex.x, thisIndex.y, actionType, 0);
    else if(thisIndex.y < thisIndex.x)   //we have the answer scalc wants
      symmSectionMgr.sendResults(m*n, A, 0, thisIndex.y, thisIndex.x, actionType, 0);
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
    symmSectionMgr.sendResults(m*n, A, 0, thisIndex.x, thisIndex.y, actionType, 0);
#ifdef _CP_ORTHO_DUMP_TMAT_
	dumpMatrix("tmatT",(double *)A, m, n,numGlobalIter,thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);

#endif
      }


#ifdef _CP_SUBSTEP_TIMING_
    if(timeKeep>0)
      {
	double oend=CmiWallTimer();
	CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
	contribute(sizeof(double),&oend,CkReduction::max_double, cb , timeKeep);
      }
#endif
    /* this should be triggered by psi after sym fw path is done */
  if(cfg.isDynamics && !config.PCCollectTiles)
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
  //  CkPrintf("[%d,%d] sending orthoT\n",thisIndex.x,thisIndex.y);
  asymmSectionMgr.sendMatrix(m*n, orthoT, 0, thisIndex.x, thisIndex.y, actionType, asymmSectionMgr.msgPriority-1);

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::maxCheck(CkReductionMsg *msg){
//============================================================================

  double tolMax=fabs(((double *) msg->getData())[0]);
  delete msg;

  CkPrintf("{%d} SMAT tol    = %g\n", cfg.instanceIndex,tolMax);
  if(tolMax < cfg.maxTolerance){
      toleranceCheckOrthoT=false;
      thisProxy.do_iteration();
  }else{
      // smat is outside of the tolerance range  need new PsiV
      toleranceCheckOrthoT=true;
      CkPrintf("recalculating PsiV due to tolerance failure \n");
      // Use the callback to trigger tolerance failure events.
      cfg.uponToleranceFailure.send();
      // Simply suspend all work. We'll be resumed (by GSpaceDriver) when tolerance updates are done.
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
void Ortho::acceptSectionLambda(CkReductionMsg *msg) {
//============================================================================

  internalType *lambda = (internalType*)msg->getData();
  int lambdaCount = msg->getSize()/sizeof(internalType);

#ifdef _CP_ORTHO_DUMP_LMAT_
    dumpMatrix("lmat",lambda, m, n, numGlobalIter, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);

#endif
#ifdef _CP_ORTHO_DEBUG_COMPARE_LMAT_
  if(savedlmat==NULL)
    { // load it
      savedlmat= new double[m*n];
      loadMatrix("lmat",(double *)savedlmat, m, n,numGlobalIter, thisIndex.x*cfg.grainSize,thisIndex.y*cfg.grainSize,0,false);
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
  if(cfg.isDynamics){

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
      ortho= new internalType[m*n];
    CmiMemcpy(ortho,orthoT,m*n*sizeof(internalType));
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
          CkAssert( isfinite(lambda[i]) );
#endif

      // finish pair calc
      asymmSectionMgr.sendResults(lambdaCount, lambda, 0, thisIndex.x, thisIndex.y, 0, asymmSectionMgr.msgPriority+1);
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
void Ortho::makeSections(const pc::pcConfig &cfgSymmPC, const pc::pcConfig &cfgAsymmPC, CkArrayID symAID, CkArrayID asymAID)
{
    /** For runs using a large numPE, Orthos chares are typically mapped onto a small fraction of the cores
     * However, array broadcasts in charm++ involve all PEs (due to some legacy quirk present because of any-time 
     * migration support). Hence, to avoid this anti-scaling overhead, we delegate array collectives to comlib / CkMulticast
     * by building a section that includes all chare elements! :)
     *
     * The (0,0) Ortho chare sets up these sections and delegates them
     */
    if( (thisIndex.x==0 && thisIndex.y==0) && (config.useOrthoSection || config.useOrthoSectionRed))
    {
        /// Create an array section that includes the whole Ortho chare array
        int numOrtho = cfg.numStates / cfg.grainSize;
        multiproxy = CProxySection_Ortho::ckNew(thisProxy.ckGetArrayID(), 0, numOrtho-1,1, 0, numOrtho-1, 1);
        /// If reductions are being delegated to a comm library
        if(config.useOrthoSectionRed)
        {
            CProxySection_Ortho rproxy =   multiproxy;
            CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(oRedGID).ckLocalBranch();
            CkAssert(mcastGrp != NULL);
            rproxy.ckSectionDelegate(mcastGrp);
            initCookieMsg *redMsg=new initCookieMsg;
            /// Ask the rest of the section (the whole array) to init their CkSectionInfo cookies that identify the mgr etc
            rproxy.orthoCookieinit(redMsg);
        }
        /// If multicasts are being delegated to a comm library
        if(config.useOrthoSection)
        {
            /// Use the appropriate library for multicasts
            if( config.useCommlib && config.useOrthoDirect)
            {
                #ifdef USE_COMLIB
                    ComlibAssociateProxy(orthoInstance,multiproxy);	  
                #endif
            }
            else
            {
                CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(oMCastGID).ckLocalBranch();
                CkAssert(mcastGrp != NULL);
                multiproxy.ckSectionDelegate(mcastGrp);
            }
        }
    }

    // Initialize the paircalc section managers with data from the paircalc config objects
    symmSectionMgr.init (thisIndex, cfgSymmPC , symAID , oMCastGID, oRedGID);
    asymmSectionMgr.init(thisIndex, cfgAsymmPC, asymAID, oMCastGID, oRedGID);
    
    /// Symmetric PC sections should trigger S -> T computations in Ortho via this method
    CkCallback orthoCB(CkIndex_Ortho::start_calc(NULL), thisProxy(thisIndex.x, thisIndex.y));
    /// Asymmetric sections should simply drop off lambda at this method
    CkCallback orthoLambdaCB(CkIndex_Ortho::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y));

    /// Setup the symmetric instance paircalc array section for communication with the symm PC chares
    symmSectionMgr.setupArraySection(orthoCB,config.phantomSym,config.useOrthoDirect);
    /// Setup the asymmetric instance paircalc array section for gather/scatter of lambda data from/to the asymm PC chares
    asymmSectionMgr.setupArraySection(orthoLambdaCB, config.phantomSym, config.useOrthoDirect);
}


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Ortho::gamma_done(){
//============================================================================
//  CkPrintf("[%d %d] sending ortho %g %g %g %g gamma %g %g %g
//%g\n",thisIndex.x,
//thisIndex.y,orthoT[0],orthoT[1],orthoT[m*n-2],orthoT[m*n-1],B[0],B[1],B[m*n-2],B[m*n-1]);

#ifdef _CP_ORTHO_DUMP_GMAT_
    dumpMatrix("gmat",(double *)B, m, n, numGlobalIter, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);

#endif

    //CkAbort("HEY I ONLY WANT ONE DYNAMICS STEP");

#ifdef _CP_ORTHO_DEBUG_COMPARE_GMAT_
  if(savedgmat==NULL)
    { // load it
      savedgmat= new double[m*n];
      loadMatrix("gmat",(double *)savedgmat, m, n,numGlobalIter,thisIndex.x*cfg.grainSize,thisIndex.y*cfg.grainSize,0,false);
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
      asymmSectionMgr.sendResults(m*n, B, ortho, thisIndex.x, thisIndex.y, 0, asymmSectionMgr.msgPriority+1);
    }
  else // orthoT was already sent ahead for processing
    {
      asymmSectionMgr.sendResults(m*n, B, 0, thisIndex.x, thisIndex.y, 0, asymmSectionMgr.msgPriority+1);
    }

//----------------------------------------------------------------------------
  }// end routine
//==============================================================================




//============================================================================
//New functions necessary for using CLA_Matrix
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

Ortho::Ortho(int _m, int _n, CLA_Matrix_interface _matA1,
 CLA_Matrix_interface _matB1, CLA_Matrix_interface _matC1,
 CLA_Matrix_interface _matA2, CLA_Matrix_interface _matB2,
 CLA_Matrix_interface _matC2, CLA_Matrix_interface _matA3,
 CLA_Matrix_interface _matB3, CLA_Matrix_interface _matC3,
 orthoConfig &_cfg,
 CkArrayID _step2Helper,
 int timekeep, CkGroupID _oMCastGID, CkGroupID _oRedGID) : 
    cfg(_cfg),
    oMCastGID(_oMCastGID), oRedGID(_oRedGID),
    step2Helper(_step2Helper)
{

/* do basic initialization */
  parallelStep2=config.useOrthoHelpers;
  
  invsqr_max_iter=config.invsqr_max_iter;
  invsqr_tolerance=config.invsqr_tolerance;
  if(invsqr_tolerance==0)
    invsqr_tolerance=INVSQR_TOLERANCE;
  if(invsqr_max_iter==0)
    invsqr_max_iter=INVSQR_MAX_ITER;
  this->matA1 = _matA1; this->matB1 = _matB1; this->matC1 = _matC1;
  this->matA2 = _matA2; this->matB2 = _matB2; this->matC2 = _matC2;
  this->matA3 = _matA3; this->matB3 = _matB3; this->matC3 = _matC3;
  timeKeep=timekeep;
  int borderOrtho = cfg.numStates / cfg.grainSize - 1;
  int remOrtho = cfg.numStates%cfg.grainSize;
  if(thisIndex.x==borderOrtho)
    this->m = _m + remOrtho;
  else
    this->m = _m;
  if(thisIndex.y==borderOrtho)
    this->n = _n + remOrtho;
  else
    this->n = _n;
  A = new internalType[this->m * this->n];
  B = new internalType[this->m * this->n];
  C = new internalType[this->m * this->n];
  tmp_arr = new internalType[this->m * this->n];
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
}//end routine




/** start step 1 on proxy 1, the callback will be for step 2
 * S2 = 3 * I - T * S1
 * currently A has T, B has S1, need to construct 3*I in C
 */
void Ortho::do_iteration(void){
    #ifdef VERBOSE_ORTHO
        CkPrintf("[%d,%d] Ortho::do_iteration \n", thisIndex.x, thisIndex.y);
    #endif
  step = 1;
  memset(C, 0, m * n * sizeof(internalType));
  if(thisIndex.x == thisIndex.y){
    for(int i = 0; i < n; i++)
      C[i * n + i] = 3;
  }
#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrix("step1:A:",(double *)A, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
    dumpMatrix("step1:B:",(double *)B, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
    dumpMatrix("step1:C:",(double *)C, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
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




/* S1 = 0.5 * S3 * S2 (on proxy 2)
 * currently A has T, B has S1, C has S2
 * Multiply tmp_arr = B*C
 * tmp_arr not used in step3, therefore no data dependence
 */
void Ortho::step_2(void){

    step = 2;
  if(config.useOrthoHelpers)
    {
    // Send our data to the helper and await results which will arrive in recvStep2
    OrthoHelperMsg *omsg= new (m*n, m*n, 0) OrthoHelperMsg;
    omsg->init(m*n, B,C,0.5, 0.5, 0.5);
    step2Helper(thisIndex.x,thisIndex.y).recvAB(omsg);
      step_3();
    }
  else
    {
#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrix("step2:A:",(double *)B, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
    dumpMatrix("step2:B:",(double *)C, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
    //    dumpMatrix("step2:C:",(double *)tmp_arr, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
#endif

      matA2.multiply(0.5, 0, B, Ortho::step_3_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      matB2.multiply(0.5, 0, C, Ortho::step_3_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
      matC2.multiply(0.5, 0, tmp_arr, Ortho::step_3_cb, (void*) this,
		     thisIndex.x, thisIndex.y);
    }
}



/**
 * T = 0.5 * S2 * S3 (do S3 = T before) (on proxy 3)
 * currently A has T, B has S1 (old), C has S2, tmp_arr has new S1
 */
void Ortho::step_3(){
  step = 3;
  CmiMemcpy(B, A, m * n * sizeof(internalType));
#ifdef _CP_ORTHO_DUMP_SMAT_
    dumpMatrix("step3:A:",(double *)C, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
    dumpMatrix("step3:B:",(double *)B, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
    dumpMatrix("step3:C:",(double *)A, m, n, iterations, thisIndex.x * cfg.grainSize, thisIndex.y * cfg.grainSize, 0, false);
#endif

  matA3.multiply(0.5, 0, C, Ortho::tol_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matB3.multiply(0.5, 0, B, Ortho::tol_cb, (void*) this,
   thisIndex.x, thisIndex.y);
  matC3.multiply(0.5, 0, A, Ortho::tol_cb, (void*) this,
   thisIndex.x, thisIndex.y);
}
//============================================================================


/** calculate error and reset pointers (for step 1)
 * current: T -> A, S1 -> tmp_arr, S2 -> C, S3 -> B
 */
void Ortho::tolerance_check(){
  step = 4;
  step2done=false;
  step3done=false;
   
#ifdef _NAN_CHECK_ 
  for(int i = 0; i < m * n; i++)
  {
      CkAssert( isfinite(A[i]) );
      CkAssert( isfinite(B[i]) );
  }
#endif

  internalType ret = 0;
  for(int i = 0; i < m * n; i++){
    internalType tmp = B[i] - A[i];
    ret += tmp * tmp;
  }
  internalType *tmp = B;
  B = tmp_arr;
  tmp_arr = tmp;
  if(config.useOrthoSectionRed)
    {
      CkCallback mycb(CkIndex_Ortho::collect_error(NULL), thisProxy(0, 0));
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(symmSectionMgr.orthomCastGrpID).ckLocalBranch();
      mcastGrp->contribute(sizeof(internalType),  &ret, CkReduction::sum_double, orthoCookie, mycb);
    }
  else
    contribute(sizeof(internalType), &ret, CkReduction::sum_double);
  iterations++;
}




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
    p | step2done;
    p | step3done;
    if(p.isUnpacking() && thisIndex.x==0 &&thisIndex.y==0)
      { 
	ortho = new internalType[cfg.numStates * cfg.numStates];
	orthoT = new internalType[cfg.numStates * cfg.numStates];
	wallTimeArr = new double[config.maxIter];
      }
    if(thisIndex.x==0 && thisIndex.y==0)
      {
	PUParray(p,ortho,cfg.numStates*cfg.numStates);
	PUParray(p,orthoT,cfg.numStates*cfg.numStates);
	p(wallTimeArr,config.maxIter);
      }
    if(p.isUnpacking()){
      A = new internalType[m * n];
      B = new internalType[m * n];
      C = new internalType[m * n];
      tmp_arr = new internalType[m * n];
    }
    PUParray(p,A, m * n);
    PUParray(p,B, m * n);
    PUParray(p,C, m * n);
    PUParray(p,tmp_arr, m * n);
  }

#include "orthoMap.h"
#include "ortho.def.h"

