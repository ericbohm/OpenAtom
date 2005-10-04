
#include "ckPairCalculator.h"
#include "pairCalculator.h"
extern ComlibInstanceHandle mcastInstanceCP;
/***************************************************************************
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
 *                                                                         *
 ***************************************************************************/

void createPairCalculator(bool sym, int s, int grainSize, int numZ, int* z, int op1, 
                          FuncType f1, int op2, FuncType f2, CkCallback cb, 
                          PairCalcID* pcid, int cb_ep, CkArrayID cb_aid, 
                          int comlib_flag, CkGroupID *mapid, int flag_dp, bool conserveMemory, bool lbpaircalc, int priority, CkGroupID mCastGrpId) {

  traceRegisterUserEvent("calcpairDGEMM", 210);
  traceRegisterUserEvent("calcpairContrib", 220);
  traceRegisterUserEvent("multiplyResultDGEMM1", 230);
  traceRegisterUserEvent("multiplyResultDGEMM2", 240);
  traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

  //CkPrintf("create pair calculator %d, %d\n", s, grainSize);

  CProxy_PairCalcReducer pairCalcReducerProxy = CProxy_PairCalcReducer::ckNew();

  CkCallback rcb = CkCallback(CkIndex_PairCalcReducer::__idx_startMachineReduction_void,
                              pairCalcReducerProxy.ckGetGroupID());
  
  pairCalcReducerProxy.ckSetReductionClient(&rcb); 

  // FIXME: to choose grain size and block size instead of hard-coding
  int blkSize = 1;
  
  CkArrayOptions options;
  CProxy_PairCalculator pairCalculatorProxy;
  redtypes cpreduce=sparsecontiguous;
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
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(sym, grainSize, s, blkSize, op1, f1, op2, f2, cb, pairCalcReducerProxy.ckGetGroupID(), cb_aid, cb_ep, conserveMemory, lbpaircalc, cpreduce,options);
  }

  pairCalculatorProxy.ckSetReductionClient(&cb);  

  int proc = 0;
  /*
    CharmStrategy * pstrat = new PipeBroadcastStrategy(USE_HYPERCUBE, pairCalculatorProxy.ckGetArrayID());
  */
  CharmStrategy *bstrat = new BroadcastStrategy(USE_HYPERCUBE);

  ComlibInstanceHandle bcastInstance = CkGetComlibInstance();
  bcastInstance.setStrategy(bstrat);
  CharmStrategy *ringstrat = new RingMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  pcid->Init(pairCalculatorProxy.ckGetArrayID(), pairCalcReducerProxy.ckGetGroupID(), grainSize, blkSize, s, sym, comlib_flag, bcastInstance, flag_dp, conserveMemory, lbpaircalc, mCastGrpId, priority);
  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  if(sym)// cheap hack to only do this once 
    mcastInstanceCP=ComlibRegister(multistrat);

  pcid->rinst=ComlibRegister(ringstrat);

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
		insert(sym, grainSize, s, blkSize, op1, f1, op2, f2, cb, pairCalcReducerProxy.ckGetGroupID(), cb_aid, cb_ep, conserveMemory, lbpaircalc, cpreduce );
	    }
	    else
	      {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, blkSize, op1, f1, op2, f2, cb, pairCalcReducerProxy.ckGetGroupID(), cb_aid, cb_ep, conserveMemory, lbpaircalc, cpreduce, proc);
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
		  insert(sym, grainSize, s, blkSize, op1, f1, op2, f2, cb, pairCalcReducerProxy.ckGetGroupID(), cb_aid, cb_ep, conserveMemory, lbpaircalc,  cpreduce );
	      }
	      else{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, blkSize, op1, f1, op2, f2, cb, pairCalcReducerProxy.ckGetGroupID(), cb_aid, cb_ep, conserveMemory, lbpaircalc,   cpreduce, proc);
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



/**
 * This is the one at a time version.  Instead of making them all at once,
 * we have each gspaceplane chare initialize its own section.
 * Each section will have S/grainsize members.  Such that PC(w,*,y,*)
 * contribute to GSP(y,w).  Symmetric case will additionally have
 * PC(w,x,y!=x,*) contributing to GSP(x,w) to fill out the total
 * S/grainsize contributions in each section.
 *
 * Return the section proxy
 */
CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  int ecount=0;
  int offset=state%pcid->GrainSize;
  int s2=state/pcid->GrainSize*pcid->GrainSize;
  CkArrayIndexMax *elems= new CkArrayIndexMax[pcid->S/pcid->GrainSize];
  for (int s1 = 0; s1 < pcid->S; s1 += pcid->GrainSize) {
    for (int c = 0; c < pcid->BlkSize; c++) {
      CkArrayIndex4D idx4d(plane,s1,s2,c);
      elems[ecount++]=idx4d;
#ifdef _PAIRCALC_DEBUG_SPROXY_
      CkPrintf("[%d %d] proxy to %d %d %d %d sym %d\n",state, plane, plane, s1, s2,c, pcid->Symmetric);  
#endif
    }
  }
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  elems, ecount);
  sectProxy.ckSectionDelegate(mcastGrp);
  //initialize proxy
  initResultMsg *redMsg=new initResultMsg;
  redMsg->mCastGrpId=pcid->mCastGrpId;
  //need some way to distinguish between row and col cookies
  //during contribute  use the base index
  redMsg->dest=s2;
  redMsg->offset=offset;
  sectProxy.initResultSection(redMsg);
  delete [] elems;
  return sectProxy;
}


CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  int ecount=0;
  int offset=state%pcid->GrainSize;
  int s2=state/pcid->GrainSize*pcid->GrainSize;
  CkArrayIndexMax *elems= new CkArrayIndexMax[pcid->S/pcid->GrainSize];
  for (int s1 = 0; s1 <=s2; s1 += pcid->GrainSize) {
    for (int c = 0; c < pcid->BlkSize; c++) {
	CkArrayIndex4D idx4d(plane,s1,s2,c);
	elems[ecount++]=idx4d;
#ifdef _PAIRCALC_DEBUG_SPROXY_
	  CkPrintf("[%d %d] proxy to %d %d %d %d sym %d S %d\n", state, plane, plane, s1, s2,c, pcid->Symmetric,pcid->S);  
#endif
    }
  }
  if(ecount)
    {
      CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid, elems, ecount);
      sectProxy.ckSectionDelegate(mcastGrp);
      //initialize proxy
      initResultMsg *redMsg=new initResultMsg;
      redMsg->mCastGrpId=pcid->mCastGrpId;
      //need some way to distinguish between row and col cookies
      //during contribute  use the base index
      redMsg->dest=s2;
      redMsg->offset=offset;
      sectProxy.initResultSection(redMsg);
      delete [] elems;
      return sectProxy;
    }
  else
    {
      CkPrintf("[%d %d] no elements for sym1?\n",state,plane);
    }
}

CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  int ecount=0;
  int offset=state%pcid->GrainSize;
  int s1=state/pcid->GrainSize*pcid->GrainSize; //column
  CkArrayIndexMax *elems= new CkArrayIndexMax[pcid->S/pcid->GrainSize];
  for (int s2 = s1+pcid->GrainSize; s2 <  pcid->S; s2 += pcid->GrainSize) {
    for (int c = 0; c < pcid->BlkSize; c++) {
	CkArrayIndex4D idx4d(plane,s1,s2,c);
	elems[ecount++]=idx4d;
#ifdef _PAIRCALC_DEBUG_SPROXY_
	CkPrintf("[%d %d] other proxy to %d %d %d %d sym %d S%d\n",state, plane, plane, s1, s2,c, pcid->Symmetric,pcid->S);  
#endif
    }
  }
  if(ecount)
    {
      CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid, elems, ecount);
      sectProxy.ckSectionDelegate(mcastGrp);
      //initialize proxy
      initResultMsg *redMsg=new initResultMsg;
      redMsg->mCastGrpId=pcid->mCastGrpId;
      //need some way to distinguish between row and col cookies
      //during contribute  use the base index
      redMsg->dest=s1;
      redMsg->offset=offset;
      sectProxy.initResultSection(redMsg);
      delete [] elems;
      return sectProxy;
    }
  else
    {
      CkPrintf("[%d %d] no elements for sym2?\n",state,plane);
    }
}

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

  /* Don't use Commlib with a section reduction proxy.  */

  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();       
  sectProxy.ckSectionDelegate(mcastGrp);
  // send the message to initialize it with the callback and groupid
  initGRedMsg *gredMsg=new initGRedMsg;
  gredMsg->cb=cb;
  gredMsg->mCastGrpId=pcid->mCastGrpId;
  sectProxy.initGRed(gredMsg);
  return sectProxy;
}


// Deposit data and start calculation
void startPairCalcLeft(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ){
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
      msgfromrow->init(n, myS, true, flag_dp, ptr);
      pcid->proxyLFrom.calculatePairs_gemm(msgfromrow);
    }
  if(pcid->existsLNotFromproxy)
    { //symmetric
      calculatePairsMsg *msg= new ( n,8*sizeof(int) ) calculatePairsMsg;
      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(msg) = pcid->priority;    
      msg->init(n, myS, false, flag_dp, ptr);   
      pcid->proxyLNotFrom.calculatePairs_gemm(msg);
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
      int erowcount=0;
      CkArrayIndexMax *elemsfromrow= new CkArrayIndexMax[blkSize*S/grainSize];
      CkArrayIndex4D idx(x,s1,0,0);
      for (c = 0; c < blkSize; c++)
	for(s2 = 0; s2 < S; s2 += grainSize){
	  idx.index[2]=s2;
	  idx.index[3]=c;
	  elemsfromrow[erowcount++]=idx;
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("Add Section ID %d: %d,%d,%d,%d \n",erowcount, idx.index[0],idx.index[1],idx.index[2],idx.index[3]);
#endif
	}
      CkAssert(erowcount>0);
#ifdef _PAIRCALC_DEBUG_
      for(int count = 0; count < erowcount; count ++) {
	CkArrayIndex4D idx4d;

	idx4d.nInts = elemsfromrow[count].nInts;
	idx4d.data()[0] = elemsfromrow[count].data()[0];
	idx4d.data()[1] = elemsfromrow[count].data()[1];
	
	CkPrintf("DEBUG ID %d: %d,%d,%d,%d \n",count, idx4d.index[0], idx4d.index[1], idx4d.index[2], idx4d.index[3]);          
      }
#endif      
      pcid->proxyLFrom = CProxySection_PairCalculator::ckNew(pairCalculatorID, elemsfromrow, erowcount); 
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
	  delete [] elemsfromrow;
	}
#endif

    }

  }

}


void isAtSyncPairCalc(PairCalcID* pcid){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     lbsync symm=%d\n", pcid->Symmetric);
#endif
  //nuke the register will rebuild in ResumeFromSync
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
      msg->init(n,myS,false,flag_dp,ptr);
      pcid->proxyRNotFrom.calculatePairs_gemm(msg);
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
      int ecount=0;
      CkArrayIndexMax *elems= new CkArrayIndexMax[blkSize*S/grainSize];
      CkArrayIndex4D idx(x,0,s2,0);
      for (c = 0; c < blkSize; c++)
	for(s1 = 0; s1 < S; s1 += grainSize){

	  idx.index[1]=s1;
	  idx.index[3]=c;
	  elems[ecount++] = idx;
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("Add Section ID %d: %d,%d,%d,%d \n",ecount, idx.index[0],idx.index[1],idx.index[2],idx.index[3]);
#endif
	}

#ifdef _PAIRCALC_DEBUG_      
      for(int count = 0; count < ecount; count ++) {
	CkArrayIndex4D idx4d;
	idx4d.nInts = elems[count].nInts;
	idx4d.data()[0] = elems[count].data()[0];
	idx4d.data()[1] = elems[count].data()[1];
	
	CkPrintf("DEBUG ID %d: %d,%d,%d,%d \n",count, idx4d.index[0], idx4d.index[1], idx4d.index[2], idx4d.index[3]);          
      }
#endif
      if(ecount)
	{
	  pcid->proxyRNotFrom = CProxySection_PairCalculator::ckNew(pairCalculatorID, elems, ecount); 
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
	      // MultiCastMgr makes its own copy
	      delete [] elems;
	    }
#endif
	}

    }
}

void finishPairCalc(PairCalcID* pcid, int n, double *ptr) {
  finishPairCalc2(pcid, n, ptr, NULL);
}


void finishPairCalc2(PairCalcID* pcid, int n, double *ptr1, double *ptr2) {
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Finish 2\n");
#endif

#ifdef _PAIRCALC_SLOW_FAT_SIMPLE_CAST_
  CkArrayID pairCalcID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalcID);
#endif

  CkGroupID pairCalcReducerID = (CkArrayID)pcid->Gid; 
  CProxy_PairCalcReducer pairCalcReducerProxy(pairCalcReducerID);

  ComlibInstanceHandle bcastInstance = pcid->cinst;
  /* if{
     #ifdef _PAIRCALC_SLOW_FAT_SIMPLE_CAST_
     ComlibDelegateProxy(&pairCalculatorProxy);
     #else
     ComlibDelegateProxy(&pairCalcReducerProxy);
     #endif
     }
  */
#ifdef _PAIRCALC_SLOW_FAT_SIMPLE_CAST_
  /* 
     Just broadcast directly to the paircalculators. We expect this to
     perform badly, this is mostly a debugging comparison block.
  */ 
  if(ptr2==NULL){
    multiplyResultMsg *omsg=new ( n,0 ) multiplyResultMsg;
    omsg->init1(n, ptr1);
    pairCalculatorProxy.multiplyResult(omsg);
  }
  else {
    multiplyResultMsg2 *omsg=new ( n,n,0 ) multiplyResultMsg2;
    omsg->init(n,n, ptr1, ptr2);
    pairCalculatorProxy.multiplyResult(omsg);
  }
#else
  if(ptr2==NULL){
    /*
      entireResultMsg *omsg=new ( n, 0 ) entireResultMsg;
      omsg->init(n, ptr1, pcid->Symmetric);
      pairCalcReducerProxy.broadcastEntireResult(omsg);
    */
    pairCalcReducerProxy.broadcastEntireResult(n, ptr1, pcid->Symmetric);
  }
  else {

    /*entireResultMsg2 *omsg=new ( n, n, 0 ) entireResultMsg2;
      omsg->init(n, ptr1, ptr2, pcid->Symmetric);

      pairCalcReducerProxy.broadcastEntireResult(omsg);
    */
    pairCalcReducerProxy.broadcastEntireResult(n, ptr1, ptr2, pcid->Symmetric);
  }
#endif

}


void finishPairCalcSection(int n, double *ptr, CProxySection_PairCalculator sectionProxy) {
  finishPairCalcSection2(n, ptr, NULL, sectionProxy);
}


/* This version uses a section multicast to only send the part of the matrix needed by each section */
void finishPairCalcSection2(int n, double *ptr1, double *ptr2, CProxySection_PairCalculator sectionProxy) {
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Finish Mcast 2\n");
#endif

  if(ptr2==NULL){
    multiplyResultMsg *omsg=new ( n,0,0 ) multiplyResultMsg;
    omsg->init1(n, ptr1);
    sectionProxy.multiplyResult(omsg);
  }
  else {
    multiplyResultMsg *omsg=new ( n,n,0 ) multiplyResultMsg;
    omsg->init(n, n, ptr1, ptr2);
    sectionProxy.multiplyResult(omsg);
  }
}


void startPairCalcLeftAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ){
    
}

void startPairCalcRightAndFinish(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ){

}



/* These are the classic no multicast version for comparison and debugging */

void startPairCalcLeftSlow(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ){
#ifdef _DEBUG_
  CkPrintf("     Calc Left ptr %d\n", ptr);

#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);

  int s1, s2, x, c;
  int grainSize = pcid->GrainSize;
  int blkSize =  pcid->BlkSize;
  int S = pcid->S;
  int symmetric = pcid->Symmetric;
  bool flag_dp = pcid->isDoublePacked;

  x = myZ;
  s1 = (myS/grainSize) * grainSize;
  if(symmetric){
    for (c = 0; c < blkSize; c++)
      for(s2 = 0; s2 < S; s2 += grainSize){
	if(s1 <= s2)
 	  {
	    calculatePairsMsg *msg=new ( n,0 ) calculatePairsMsg;
	    msg->init(n, myS, true, flag_dp, ptr);
	    pairCalculatorProxy(x, s1, s2, c).calculatePairs_gemm(msg);
	  }
	else
	  {
	    calculatePairsMsg *msg=new ( n,0 ) calculatePairsMsg;
	    msg->init(n, myS, false, flag_dp, ptr);
	    pairCalculatorProxy(x, s2, s1, c).calculatePairs_gemm(msg);
	  }
      }
  }
  else {
    for (c = 0; c < blkSize; c++)
      for(s2 = 0; s2 < S; s2 += grainSize){
	calculatePairsMsg *msg=new ( n,0 ) calculatePairsMsg;
	msg->init(n, myS, true, flag_dp, ptr);
	pairCalculatorProxy(x, s1, s2, c).calculatePairs_gemm(msg);
      }
  }
}

void startPairCalcRightSlow(PairCalcID* pcid, int n, complex* ptr, int myS, int myZ){
#ifdef _DEBUG_
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
  for (c = 0; c < blkSize; c++)
    for(s1 = 0; s1 < S; s1 += grainSize){
      calculatePairsMsg *msg=new ( n,0 ) calculatePairsMsg;
      msg->init(n, myS, false, flag_dp, ptr);
      pairCalculatorProxy(x, s1, s2, c).calculatePairs_gemm(msg);
    }
}
