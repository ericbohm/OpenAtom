//**************************************************************************
 /** \file pairCalculator.C						   *
 * This is a matrix multiply library with extra frills to communicate the  *
 * results back to gspace or the calling ortho char as directed by the     *
 * callback.                                                               *
 *                                                                         *
 * The pairCalculator handles initialization and creation of the           *
 * ckPairCalculator chare arrays, their reduction group, their multicast   *
 * manager, and their section proxies                                      *
 *                                                                         *
 * Expected usage begings with createPairCalculator(*,PairCalcID *,*)      *
 * the PairCalcID contains meta information about the calculator.          *
 * In particular, the various section proxies and array ids necessary to   *
 * handle the expected communication modalities between a parent array and *
 * the ckPairCalculator array.                                             *
 *                                                                         *
 * Folloup usage goes through:                                             *
 *  startPairCalcLeft(PairCalcID, datasize, data *, index1, index2)        *
 *                                                                         *
 * The result is returned by the callback set in the create routine        * 
 * The backward path is trigered by:                                       *
 *                                                                         *
 * finishPairCalc(PairCalcID, datasize, data *)                            *
 *  Its result is returned via the end entry point which was also set      *
 *   during creation                                                       *
 *                                                                         *
 * The results of the backward path are returned in a set of section       *
 * reductions.  The reduction sums a matrix of doubles across a            *
 * section of the paircalculator to the scatter client.  The scatter       *
 * client then sends slices of the result to the appropriate chares in     *
 * gspace.                                                                 *
 *                                                                         */
//**************************************************************************
 
#include "ckPairCalculator.h"
#include "pairCalculator.h"
extern ComlibInstanceHandle mcastInstanceCP;


void createPairCalculator(bool sym, int s, int grainSize, int numZ, int* z, 
			  CkCallback cb,  PairCalcID* pcid, int cb_ep, 
			  int cb_ep_tol, 
			  CkArrayID cb_aid, int comlib_flag, CkGroupID *mapid,
			  int flag_dp, bool conserveMemory, bool lbpaircalc, 
			  int priority, CkGroupID mCastGrpId) {

  traceRegisterUserEvent("calcpairDGEMM", 210);
  traceRegisterUserEvent("calcpairContrib", 220);
  traceRegisterUserEvent("multiplyResultDGEMM1", 230);
  traceRegisterUserEvent("multiplyResultDGEMM2", 240);
  traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

  //CkPrintf("create pair calculator %d, %d\n", s, grainSize);

  /*
   *CProxy_PairCalcReducer pairCalcReducerProxy = CProxy_PairCalcReducer::ckNew();
   *
   * CkCallback rcb = CkCallback(CkIndex_PairCalcReducer::__idx_startMachineReduction_void,
   *                            pairCalcReducerProxy.ckGetGroupID());
   *pairCalcReducerProxy.ckSetReductionClient(&rcb); not using this anymore
   */

  // FIXME: nuke this block size and unused 4th dimension BS
  int blkSize = 1;
  
  CkArrayOptions options;
  CProxy_PairCalculator pairCalculatorProxy;
  redtypes cpreduce=section;

#ifdef CONVERSE_VERSION_ELAN
  bool machreduce=(s/grainSize * numZ* blkSize>=CkNumNodes()) ? true: false;
#else
  bool machreduce=false;
#endif
  
  if(machreduce)
    cpreduce=machine;
  if(!mapid) {
    pairCalculatorProxy = CProxy_PairCalculator::ckNew();
  }
  else {
    options.setMap(*mapid);
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(sym, grainSize, s, blkSize,  cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc, cpreduce,options);
  }

  int proc = 0;

  pcid->Init(pairCalculatorProxy.ckGetArrayID(), grainSize, blkSize, s, sym, comlib_flag, flag_dp, conserveMemory, lbpaircalc, mCastGrpId, priority);

  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  if(sym)// cheap hack to only do this once 
    mcastInstanceCP=ComlibRegister(multistrat);


  if(sym)
    for(int numX = 0; numX < numZ; numX += blkSize){
      for (int s1 = 0; s1 < s; s1 += grainSize) {
	for (int s2 = s1; s2 < s; s2 += grainSize) {
	  for (int c = 0; c < blkSize; c++) {
	    if(mapid) {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
	      pairCalculatorProxy(z[numX],s1,s2,c).
		insert(sym, grainSize, s, blkSize,  cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc, cpreduce );
	    }
	    else
	      {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, blkSize, cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc, cpreduce, proc);
		proc++;
		if (proc >= CkNumPes()) proc = 0;
	      }
	  }
	}
      }
    }
  else
    {
      for(int numX = 0; numX < numZ; numX += blkSize){
	for (int s1 = 0; s1 < s; s1 += grainSize) {
	  for (int s2 = 0; s2 < s; s2 += grainSize) {
	    for (int c = 0; c < blkSize; c++) {
	      if(mapid) {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, blkSize, cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc,  cpreduce );
	      }
	      else{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, blkSize,  cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc,   cpreduce, proc);
		proc++;
		if (proc >= CkNumPes()) proc = 0;
	      }
	    }
	  }
	}          
      }
    }
  pairCalculatorProxy.doneInserting();
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("    Finished init {grain=%d, sym=%d, blk=%d, Z=%d, S=%d}\n", grainSize, sym, blkSize, numZ, s);
#endif
}




//! initialize  plane and row wise section reduction for lambda->gspace
/**
 * The makeOneResultSection functions all have the same mission.  Make one
 * section at a time using only the relevant processors instead of making them
 * all at once. We have each gspaceplane chare initialize its own section.
 * Each section will have S/grainsize members.  Such that PC(w,*,y,*)
 * contribute to GSP(y,w).  Symmetric case will additionally have
 * PC(w,x,y!=x,*) contributing to GSP(x,w) to fill out the total S/grainsize
 * contributions in each section.
 *
 * Then return the section proxy.
 */
CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  int ecount=0;
  int offset=state%pcid->GrainSize;
  int s2=state/pcid->GrainSize*pcid->GrainSize;
  int S=pcid->S;
  int GrainSize=pcid->GrainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  
									       plane, plane, 1,
									       0, S-GrainSize, GrainSize,
									       s2, s2, 1,
									       0, 0, 1);
  sectProxy.ckSectionDelegate(mcastGrp);
  //initialize proxy
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * initialize  plane and row wise section reduction for psi->gspace
 */
CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  int offset=state%pcid->GrainSize;
  int s2=state/pcid->GrainSize*pcid->GrainSize; //row
  int S=pcid->S;
  int GrainSize=pcid->GrainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  
									       plane, plane, 1,
									       0, s2, GrainSize,
									       s2, s2, 1,
									       0, 0, 1);
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}


/**
 * initialize  plane and column wise section reduction for psi->gspace
 */
CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  int offset=state%pcid->GrainSize;
  int s1=state/pcid->GrainSize*pcid->GrainSize; //column
  int S=pcid->S;
  int GrainSize=pcid->GrainSize;
  CkAssert(s1+GrainSize<S);
  CProxySection_PairCalculator sectProxy = 
      CProxySection_PairCalculator::ckNew(pcid->Aid,  
					  plane, plane, 1,
					  s1, s1, 1,
					  s1+GrainSize, S-GrainSize, GrainSize,
					  0, 0, 1);
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * initialize the planewise section reduction for Ortho
 */
CProxySection_PairCalculator initOneRedSect(int numZ, int* z, int blkSize,  PairCalcID* pcid, CkCallback cb, int s1, int s2, int c)
{
  int ecount=0;
  CkArrayIndexMax *elems= new CkArrayIndexMax[numZ/blkSize];
  for(int numX = 0; numX < numZ; numX += blkSize){
    CkArrayIndex4D idx4d(z[numX],s1,s2,c);
    elems[ecount++]=idx4d;
  }

  // now that we have the section, make the proxy and do delegation
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  elems, ecount); 
  delete [] elems;

  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  sectProxy.ckSectionDelegate(mcastGrp);

  // send the message to initialize it with the callback and groupid
  setGredProxy(&sectProxy, pcid->mCastGrpId, cb, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * send the multcast message to initialize the ortho section tree and set the cookie
 */
void setGredProxy(CProxySection_PairCalculator *sectProxy, CkGroupID mCastGrpId, CkCallback cb, bool lbsync, CkCallback synccb)
{
  initGRedMsg *gredMsg=new initGRedMsg;
  gredMsg->cb=cb;
  gredMsg->mCastGrpId=mCastGrpId;
  gredMsg->lbsync=lbsync;
  gredMsg->synccb=synccb;
  sectProxy->initGRed(gredMsg);
}

/**
 * send the multcast message to initialize the section tree and set the cookie
 */
void setResultProxy(CProxySection_PairCalculator *sectProxy, int state, int GrainSize, CkGroupID mCastGrpId, bool lbsync, CkCallback synccb)
{
    int offset=state%GrainSize;
    int dest=state/GrainSize*GrainSize; //row or column
    initResultMsg *redMsg=new initResultMsg;
    redMsg->mCastGrpId=mCastGrpId;
    redMsg->dest=dest;
    redMsg->offset=offset;
    redMsg->lbsync=lbsync;
    redMsg->synccb=synccb;
    sectProxy->initResultSection(redMsg);
}

// Deposit data and start calculation
void startPairCalcLeft(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ, bool psiV){
#ifdef _PAIRCALC_NO_MULTI_
  startPairCalcLeftSlow(pcid, n, ptr, myS, myZ);
#else
  int symmetric = pcid->Symmetric;
  bool flag_dp = pcid->isDoublePacked;
  if(!(pcid->existsLproxy||pcid->existsLNotFromproxy)){
    makeLeftTree(pcid,myS,myZ);
  }
  //use proxy to send
#ifdef _PAIRCALC_DEBUG_PARANOID_
  double re;
  double im;
  for(int i=0;i<n;i++)
    {
      re=ptr[i].re;
      im=ptr[i].im;
      if(fabs(re)>0.0)
	CkAssert(fabs(re)>1.0e-300);
      if(fabs(im)>0.0)
	CkAssert(fabs(im)>1.0e-300);
    }
#endif
  if(pcid->existsLproxy)
    {
      calculatePairsMsg *msgfromrow=new (n, 8* sizeof(int)) calculatePairsMsg;
      CkSetQueueing(msgfromrow, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msgfromrow) = pcid->priority;    
      msgfromrow->init(n, myS, true, flag_dp, ptr, psiV);
      pcid->proxyLFrom.acceptPairData(msgfromrow);
    }
  if(pcid->existsLNotFromproxy)
    { //symmetric
      calculatePairsMsg *msg= new ( n,8*sizeof(int) ) calculatePairsMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = pcid->priority;    
      msg->init(n, myS, false, flag_dp, ptr, psiV);   
      pcid->proxyLNotFrom.acceptPairData(msg);
    }


  if(pcid->useComlib && _PC_COMMLIB_MULTI_) {
  }  
#endif
}


// create multicast proxies
void makeLeftTree(PairCalcID* pcid, int myS, int myZ){
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);
  int s1, s2, x, c;
  int grainSize = pcid->GrainSize;
  int blkSize =  pcid->BlkSize;
  int S = pcid->S;
  int symmetric = pcid->Symmetric;
  bool flag_dp = pcid->isDoublePacked;

  bool conserveMemory = pcid->conserveMemory;
  x = myZ;
  s1 = (myS/grainSize) * grainSize;
  if(!(pcid->existsLproxy||pcid->existsLNotFromproxy)){
    int numElems;
    //create multicast proxy array section list 
    if(symmetric){
      CkArrayIndexMax *elems= new CkArrayIndexMax[blkSize*S/grainSize];
      CkArrayIndexMax *elemsfromrow= new CkArrayIndexMax[blkSize*S/grainSize];
      // 1 proxy for left and 1 for right
      int erowcount=0;
      int ecount=0;
      CkArrayIndex4D idx(x,0,0,0);
      for (c = 0; c < blkSize; c++)
	for(s2 = 0; s2 < S; s2 += grainSize){
	  if(s1 <= s2)
	    {
	      idx.index[1]=s1;
	      idx.index[2]=s2;
	      idx.index[3]=c;
	      elemsfromrow[erowcount++]=idx;
	    }
	  else // swap s1 : s2 and toggle fromRow
	    {
	      idx.index[1]=s2;
	      idx.index[2]=s1;
	      idx.index[3]=c;
	      elems[ecount++]=idx;
	    }
	}
      if(ecount)
	{

	  pcid->proxyLNotFrom = CProxySection_PairCalculator::ckNew(pairCalculatorID, elems, ecount); 
	  pcid->existsLNotFromproxy=true;	  
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
	  if(pcid->useComlib && _PC_COMMLIB_MULTI_)
	    {
	      ComlibAssociateProxy(&mcastInstanceCP,pcid->proxyLNotFrom);
	    }
	  else
	    {

	      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
	      pcid->proxyLNotFrom.ckSectionDelegate(mcastGrp);
	    }
#endif
	}
      if(erowcount)
	{
	  pcid->proxyLFrom  = CProxySection_PairCalculator::ckNew(pairCalculatorID, elemsfromrow, erowcount); 
	  pcid->existsLproxy=true;	  
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
	  if(pcid->useComlib && _PC_COMMLIB_MULTI_)
	    {
	      ComlibAssociateProxy(&mcastInstanceCP, pcid->proxyLFrom);
	    }
	  else
	    {
	      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch(); 
	      pcid->proxyLFrom.ckSectionDelegate(mcastGrp);
	      // MultiCastMgr makes its own copy
	      delete [] elemsfromrow;
	      delete [] elems;

	    }
#endif
	}
    }
    else { //just make left here, right will be taken care of in startRight
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("initializing multicast proxy in %d %d \n",x,s1);
#endif
      pcid->proxyLFrom = CProxySection_PairCalculator::ckNew(pcid->Aid,  
							     x, x, 1,
							     s1, s1, 1,
							     0, S-grainSize, grainSize,
							     0, 0, 1);
      pcid->existsLproxy=true;      

#ifndef _PAIRCALC_DO_NOT_DELEGATE_

      if(pcid->useComlib && _PC_COMMLIB_MULTI_ )
	{
	  ComlibAssociateProxy(&mcastInstanceCP,pcid->proxyLFrom);
	}
      else
	{
	  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch(); 
	  pcid->proxyLFrom.ckSectionDelegate(mcastGrp);
	  // MultiCastMgr makes its own copy
	}
#endif

    }

  }

}

/**
 * synchronize for migration
 */
void isAtSyncPairCalc(PairCalcID* pcid){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     lbsync symm=%d\n", pcid->Symmetric);
#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);
  pairCalculatorProxy.lbsync();
}

void startPairCalcRight(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ){
#ifdef _PAIRCALC_NO_MULTI_
  startPairCalcRightSlow(pcid, n, ptr, myS, myZ);
#else
  bool flag_dp = pcid->isDoublePacked;
  if(!pcid->existsRproxy)
    {
      makeRightTree(pcid,myS,myZ);
    }
  if(pcid->existsRproxy)
    {
#ifdef _DEBUG_PAIRCALC_PARANOID_
      double re;
      double im;
      for(int i=0;i<n;i++)
	{
	  re=ptr[i].re;
	  im=ptr[i].im;
	  if(fabs(re)>0.0)
	    CkAssert(fabs(re)>1.0e-300);
	  if(fabs(im)>0.0)
	    CkAssert(fabs(im)>1.0e-300);
	}
#endif
      calculatePairsMsg *msg= new ( n,8*sizeof(int) ) calculatePairsMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = pcid->priority;    
      msg->init(n,myS,false,flag_dp,ptr,false);
      pcid->proxyRNotFrom.acceptPairData(msg);
    }
  else
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("Warning! No Right proxy ! \n");
#endif
    }
#endif //_NO_MULTI
}

void makeRightTree(PairCalcID* pcid, int myS, int myZ){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Right symm=%d\n", pcid->Symmetric);
#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);

  int s1, s2, x, c;
  int grainSize = pcid->GrainSize;
  int blkSize =  pcid->BlkSize;
  int S = pcid->S;
  bool symmetric = pcid->Symmetric;
  bool flag_dp = pcid->isDoublePacked;

  CkAssert(symmetric == false);
  
  x = myZ;
  s2 = (myS/grainSize) * grainSize;
  //create multicast proxy list 
  if(!pcid->existsRproxy)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("initializing multicast proxy in %d %d \n",x,s2);
#endif
      pcid->proxyRNotFrom = 
	  CProxySection_PairCalculator::ckNew(pcid->Aid,  
					      x, x, 1,
					      0, S-grainSize, grainSize,
					      s2, s2, 1,
					      0, 0, 1);
      pcid->existsRproxy=true;      
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
      if(pcid->useComlib && _PC_COMMLIB_MULTI_)
      {
	  ComlibAssociateProxy(&mcastInstanceCP,pcid->proxyRNotFrom);
      }
      else
      {
	  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch(); 
	  pcid->proxyRNotFrom.ckSectionDelegate(mcastGrp);
      }
#endif
    }
}




void finishPairCalcSection(int n, double *ptr, CProxySection_PairCalculator sectionProxy, int actionType) {
  finishPairCalcSection2(n, ptr, NULL, sectionProxy, actionType);
}


/* This version uses a section multicast to only send the part of the matrix needed by each section */
void finishPairCalcSection2(int n, double *ptr1, double *ptr2, CProxySection_PairCalculator sectionProxy, int actionType) {
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Finish Mcast 2\n");
#endif

  if(ptr2==NULL){
    multiplyResultMsg *omsg=new ( n,0,0 ) multiplyResultMsg;
    omsg->init1(n, ptr1, actionType);
    sectionProxy.multiplyResult(omsg);
  }
  else {
    multiplyResultMsg *omsg=new ( n,n,0 ) multiplyResultMsg;
    omsg->init(n, n, ptr1, ptr2, actionType);
    sectionProxy.multiplyResult(omsg);
  }
}




