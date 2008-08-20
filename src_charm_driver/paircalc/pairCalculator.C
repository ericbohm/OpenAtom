//**************************************************************************
 /** \file pairCalculator.C						   
 * This is a matrix multiply library with extra frills to communicate the  
 * results back to gspace or the calling ortho char as directed by the     
 * callback.                                                               
 *                                                                         
 * The pairCalculator handles initialization and creation of the           
 * ckPairCalculator chare arrays, their reduction group, their multicast   
 * manager, and their section proxies                                      
 *                                                                         
 * Expected usage begings with createPairCalculator(*,PairCalcID *,*)      
 * the PairCalcID contains meta information about the calculator.          
 * In particular, the various section proxies and array ids necessary to   
 * handle the expected communication modalities between a parent array and 
 * the ckPairCalculator array.                                             
 *                                                                         
 * After the main::main proc 0 phase the array sections used to populate    
 * and return data from the paircalculator are called.                     
 * The forward path section reduction (to and from ortho) is initialized   
 * via initOneRedSect()  the backward path is initialized via the          
 * appropriate makeOneResultSection_X() call.  In each case the call       
 * should be made by each GSP or Ortho object.  That way each one has its  
 * own proxy and the section tree will only include relevant processors.   
 *                                                                         
 * Followup usage goes through:                                            
 *  startPairCalcLeft(PairCalcID, datasize, data *, index1, index2)        
 * and in the asymmetric case                                              
 *  startPairCalcRight(PairCalcID, datasize, data *, index1, index2)       
 *                                                                         
 * The result is returned by the callback set in the create routine         
 * The backward path is triggered by:                                      
 *                                                                         
 * finishPairCalcSection(PairCalcID, datasize, data *)                     
 *                                                                         
 *  Its result is returned via the callback entry point passed in 
 *   during creation                                                       
 *                                                                         
 * The results of the backward path are returned in a set of section
 * reductions.  The reduction contributes its slice of its matrix of
 * doubles with the offset=thisIndex.z.  The client then returns the
 * sum of each slice to the GSP element that created the section with
 * the offset so it can be copied into the correct place in the points
 * array.
 *
 * The chunk decomposition changes the setup substantially.  In chunk
 * decomposition we send a piece of the nonzero points of gspace to
 * the paircalculators.  It is not a multicast.  Each GSP[P,S] will
 * send its ith chunk to PC[P,0,0,i].  Once nstate chunks arrive at a
 * PC the multiply can proceed.
 *
 * In the hybrid case this becomes a multicast of chunks to the
 * appropriate state decomposition destination as before.
 *
 */
//*************************************************************************
#include "ckPairCalculator.h"
#include "pairCalculator.h"
#include <algorithm>
extern ComlibInstanceHandle mcastInstanceCP;
extern ComlibInstanceHandle mcastInstanceACP;
extern ComlibInstanceHandle gAsymInstance;
extern ComlibInstanceHandle gSymInstance;


void createPairCalculator(bool sym, int s, int grainSize, int numZ, int* z, 
			  CkCallback cb,  PairCalcID* pcid, int cb_ep, 
			  int cb_ep_tol, int rdma_ep,
			  CkArrayID cb_aid, int comlib_flag, CkGroupID *mapid,
			  int flag_dp, bool conserveMemory, bool lbpaircalc, 
			  int priority, CkVec <CkGroupID> mCastGrpId,
			  CkGroupID orthomCastGrpId, CkGroupID orthoRedGrpId, 
			  int numChunks, int orthoGrainSize, bool collectTiles,
			  bool streamBWout, bool delayBWSend, int streamFW, 
			  bool useDirectSend, bool gSpaceSum, int gpriority, 
			  bool phantomSym, bool useBWBarrier, 
			  int gemmSplitFWk, int gemmSplitFWm, int gemmSplitBW,
			  bool expectOrthoT, int instance)
			   {

  traceRegisterUserEvent("calcpairDGEMM", 210);
  traceRegisterUserEvent("calcpairContrib", 220);
  traceRegisterUserEvent("multiplyResultDGEMM1", 230);
  traceRegisterUserEvent("multiplyResultDGEMM2", 240);
  traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

  CkArrayOptions options;
  CProxy_PairCalculator pairCalculatorProxy;
  redtypes cpreduce=section;
  
#ifdef CONVERSE_VERSION_ELAN
  bool machreduce=(s/grainSize * numZ* numChunks>=CkNumNodes()) ? true: false;
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
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, gSpaceSum, gpriority, phantomSym, useBWBarrier, 
						       gemmSplitFWk, gemmSplitFWm, 
						       gemmSplitBW,expectOrthoT, instance,
						       options);
  }

  int proc = 0;

  pcid->Init(pairCalculatorProxy.ckGetArrayID(), grainSize, numChunks, s, sym, comlib_flag, flag_dp, conserveMemory, lbpaircalc,  priority, useDirectSend);
  pcid->orthomCastGrpId=orthomCastGrpId;
  pcid->orthoRedGrpId=orthoRedGrpId;
  pcid->cproxy=pairCalculatorProxy;
  pcid->mCastGrpId=mCastGrpId;
#ifdef CMK_BLUEGENEL
  //  CharmStrategy *multistrat = new RectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
#else
  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
#endif

  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  if(sym)
    mcastInstanceCP=ComlibRegister(multistrat);
  else
    mcastInstanceACP=ComlibRegister(multistrat);
  CkAssert(mapid);
  if(sym)
    for(int numX = 0; numX < numZ; numX ++){
      for (int s1 = 0; s1 <= maxpcstateindex; s1 += grainSize) {
	int s2start=(phantomSym) ? 0 : s1;
	for (int s2 = s2start; s2 <= maxpcstateindex; s2 += grainSize) {
	  for (int c = 0; c < numChunks; c++) {
	    if(mapid) {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
	      pairCalculatorProxy(z[numX],s1,s2,c).
		insert(sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, gSpaceSum, gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance );
	    }
	    else
	      {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, numChunks, cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, gSpaceSum, gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance, proc);
		proc++;
		if (proc >= CkNumPes()) proc = 0;
	      }
	  }
	}
      }
    }
  else
    {
      for(int numX = 0; numX < numZ; numX ++){
	for (int s1 = 0; s1 <= maxpcstateindex; s1 += grainSize) {
	  for (int s2 = 0; s2 <= maxpcstateindex; s2 += grainSize) {
	    for (int c = 0; c < numChunks; c++) {
	      if(mapid) {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, numChunks, cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc,  cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, gSpaceSum,  gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance);
	      }
	      else{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc,   cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, gSpaceSum, gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance,proc);
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
  CkPrintf("    Finished init {grain=%d, sym=%d, blk=%d, Z=%d, S=%d}\n", grainSize, sym, numChunks, numZ, s);
#endif
}




//! initialize  plane and row wise section reduction for lambda->gspace
/**
 * The makeOneResultSection functions all have the same mission.  Make one
 * section at a time using only the relevant processors instead of making them
 * all at once. We have each gspaceplane chare initialize its own section.
 * Each section will have S/grainsize members.  Such that PC(w,*,y,*)
 * contribute to GSP(y,w).  
 *
 * Symmetric and asymm dynamics case will additionally have
 * PC(w,x,y!=x,*) contributing to GSP(x,w) to fill out the total
 * S/grainsize contributions in each section.
 *
 * Then return the section proxy.
 */
CProxySection_PairCalculator makeOneResultSection_asym(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[plane]).ckLocalBranch();       
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  int s2=state/pcid->GrainSize*pcid->GrainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
  int nstates=pcid->nstates;
  int GrainSize=pcid->GrainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  
									       plane, plane, 1,
									       0, maxpcstateindex, GrainSize,
									       s2, s2, 1,
									       chunk, chunk,1);
  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart=state%GrainSize;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list_max( sid._elems, sid._nElems, newListStart);
  CkAssert(order);
  sectProxy.ckSectionDelegate(mcastGrp);
  //initialize proxy
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId[plane], false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * initialize  plane and column wise section reduction for lambda->gspace
 */
CProxySection_PairCalculator makeOneResultSection_asym_column(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[plane]).ckLocalBranch();       
  int GrainSize=pcid->GrainSize;
  int s1=state / GrainSize * GrainSize; //column
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

  int nstates=pcid->nstates;
  // all nondiagonal elements 
  // so we'll have to make this the tedious explicit way

  CkArrayIndex4D *elems= new CkArrayIndex4D[nstates/GrainSize];
  int ecount=0;
  for(int s2 =0; s2<=maxpcstateindex; s2+=GrainSize)
    {
      if(s1!=s2)
	{
	  CkArrayIndex4D idx4d(plane,s1,s2,chunk);
	  elems[ecount++]=idx4d;
	}
    }
  int newListStart=state%GrainSize;
  if(newListStart> ecount)
    newListStart= newListStart % ecount;
  bool order=reorder_elem_list_4D( elems, ecount, newListStart);
  CkAssert(order);
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  elems, ecount); 
  delete [] elems;
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId[plane], false, CkCallback(CkCallback::ignore));
  return sectProxy;
}



/**
 * initialize  plane and row wise section reduction for psi->gspace
 */
CProxySection_PairCalculator makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[plane]).ckLocalBranch();       
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  int s2=state/pcid->GrainSize*pcid->GrainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

  int GrainSize=pcid->GrainSize;
  int s2range= (s2==0) ? 1 : GrainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  
									       plane, plane, 1,
									       0, s2, s2range,
									       s2, s2, 1,
									       chunk, chunk, 1);
  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart=state%GrainSize;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list_max( sid._elems, sid._nElems, newListStart);
  CkAssert(order);
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId[plane], false, CkCallback(CkCallback::ignore));
  return sectProxy;
}


/**
 * initialize  plane and column wise section reduction for psi->gspace
 */
CProxySection_PairCalculator makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[plane]).ckLocalBranch();       
  int GrainSize=pcid->GrainSize;
  int s1=state / GrainSize * GrainSize; //column
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

  int nstates=pcid->nstates;
  int s2start=s1+GrainSize;
  s2start= (s2start>maxpcstateindex) ? maxpcstateindex : s2start;
  int s2range= (s2start==maxpcstateindex) ? 1 : GrainSize;
  CkAssert(s2start<nstates);
  CProxySection_PairCalculator sectProxy = 
      CProxySection_PairCalculator::ckNew(pcid->Aid,  
					  plane, plane, 1,
					  s1, s1, 1,
					  s2start, maxpcstateindex, s2range,
					  chunk, chunk, 1);

  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart=state%GrainSize;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list_max( sid._elems, sid._nElems, newListStart);
  CkAssert(order);
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId[plane], false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * initialize the planewise section reduction for Ortho sums across
 * all planes and chunks pass through the orthoX and orthoY so the
 * cookie can be placed in the 2d array
 * (grainSize/orthoGrainSize)^2
 */

void initOneRedSect(int numZ, int* z, int numChunks,  PairCalcID* pcid, CkCallback cb, CkCallback synccb, int s1, int s2, int orthoX, int orthoY, int orthoGrainSize, bool phantom, bool direct, bool commlib)

{
  int ecount=0;
  //  CkPrintf("initOneRedSect for s1 %d s2 %d ortho %d %d sym %d planes %d\n",s1,s2,orthoX, orthoY,pcid->Symmetric, numZ);
  CkArrayIndex4D *elems=new CkArrayIndex4D[numZ*numChunks*2];
  //add chunk loop
  for(int chunk = numChunks-1; chunk >=0; chunk--){
    for(int numX = numZ-1; numX >=0; numX--){
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("initOneRedSect for s1 %d s2 %d ortho %d %d sym %d plane %d\n",s1,s2,orthoX, orthoY,pcid->Symmetric, numX);
#endif
      if(phantom && s1!=s2)
	{
	  CkArrayIndex4D idx4d(z[numX],s1,s2,chunk);
	  elems[ecount++]=idx4d;
	}
      else
	{
	  CkArrayIndex4D idx4d(z[numX],s1,s2,chunk);
	  elems[ecount++]=idx4d;
	  //	  CkPrintf("O [%d,%d] initGred section includes %d %d %d %d sym %d\n",orthoX, orthoY,z[numX],s1,s2,chunk,pcid->Symmetric);
	}
    }
  }
  int numOrthoCol=pcid->GrainSize/orthoGrainSize;
  int maxorthostateindex=(pcid->nstates/orthoGrainSize-1)*orthoGrainSize;
  int orthoIndexX=(orthoX*orthoGrainSize);

  orthoIndexX= (orthoIndexX>maxorthostateindex) ? maxorthostateindex : orthoIndexX;
  int orthoIndexY=(orthoY*orthoGrainSize);
  orthoIndexY= (orthoIndexY>maxorthostateindex) ? maxorthostateindex : orthoIndexY;
  orthoIndexX-=s1;
  orthoIndexY-=s2;
  int orthoIndex=orthoIndexX*numOrthoCol+orthoIndexY;

  int newListStart=orthoIndex;
  if(newListStart> ecount)
    newListStart= newListStart % ecount;
  bool order=reorder_elem_list_4D( elems, ecount, newListStart);
  CkAssert(order);
  // now that we have the section, make the proxy
  CProxySection_PairCalculator sProxy=CProxySection_PairCalculator::ckNew(pcid->Aid,  elems, ecount); 
  CProxySection_PairCalculator *sectProxy=&sProxy;
  delete [] elems;

  // and do delegation
  if(pcid->Symmetric)
    {
      pcid->proxySym = sProxy;
    }
  else
    {
      pcid->proxyAsym = sProxy;
    }
  if(!phantom && !direct) // only delegating nonphantom mcast proxy for reduction
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->orthoRedGrpId).ckLocalBranch();       
      sectProxy->ckSectionDelegate(mcastGrp);
      // send the message to initialize it with the callback and groupid
      setGredProxy(sectProxy, pcid->orthoRedGrpId, cb, false, synccb, orthoX, orthoY);
    }
  else
    {
      if(commlib)
	{
	  CkPrintf("NOTE: Rectangular Send In USE\n");
	  if(pcid->Symmetric)
	    ComlibAssociateProxy(&mcastInstanceCP,*sectProxy);
	  else 
	    ComlibAssociateProxy(&mcastInstanceACP,*sectProxy);	 
 /*
	  if(!pcid->Symmetric)
	    ComlibAssociateProxy(&mcastInstanceACP,*sectProxy);
	  */
	}
      else
	{
	  //CkPrintf("PC: proxy without commlib\n");
	  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->orthomCastGrpId).ckLocalBranch();       
	  sectProxy->ckSectionDelegate(mcastGrp);

	}
    }

  //  return *sectProxy;
}

/**
 * send the multcast message to initialize the ortho section tree and set the cookie
 */
void setGredProxy(CProxySection_PairCalculator *sectProxy, CkGroupID mCastGrpId, CkCallback cb, bool lbsync, CkCallback synccb, int orthoX, int orthoY)
{
  initGRedMsg *gredMsg=new initGRedMsg;
  gredMsg->cb=cb;
  gredMsg->mCastGrpId=mCastGrpId;
  gredMsg->lbsync=lbsync;
  gredMsg->synccb=synccb;
  gredMsg->orthoX=orthoX;
  gredMsg->orthoY=orthoY;
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

void startPairCalcLeft(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV){
#ifdef PC_USE_RDMA
  if(psiV)// RDMA doesn't support PSIV shenanighans
    CkAbort("You must #undef PC_USE_RDMA in ckPairCalculator.h for dynamics. PSIV is broken by RDMA.\n It can be fixed, but won't be done until someone needs production support.");
  else
    startPairCalcLeftRDMA(pcid, n, ptr, myS, myPlane,psiV);
#else
#ifdef _PAIRCALC_NO_MULTI_
  startPairCalcLeftSlow(pcid, n, ptr, myS, myPlane);
#else
  bool flag_dp = pcid->isDoublePacked;
  if(!(pcid->existsLproxy||pcid->existsLNotFromproxy)){
    makeLeftTree(pcid,myS,myPlane);
    CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
    pcid->cproxy= CProxy_PairCalculator(pairCalculatorID);
  }
  //use proxy to send
  if(pcid->existsLproxy)
    {
      int numChunks=pcid->numChunks;
      int chunksize =  n / numChunks;
      int outsize = chunksize;

      for(int chunk=0; chunk < numChunks ; chunk++)
	{
	  if((numChunks > 1) && (chunk == (numChunks - 1)))
	    {// last chunk gets remainder
	      outsize= chunksize + (n % numChunks);
	    }
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	  if(pcid->Symmetric && myPlane==0)
	    dumpMatrixDouble("gspPts",(double *)ptr, 1, n*2,myPlane,myS,0,chunk,pcid->Symmetric);
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	  CkPrintf("L [%d,%d,%d,%d,%d] chunk %d chunksize %d outsize %d for numpoint %d offset will be %d %.12g\n",myPlane,myS, myS, chunk,pcid->Symmetric, chunk,chunksize,outsize,n,chunk*chunksize,ptr[chunk*chunksize].re);
#endif

	  if( pcid->useDirectSend)
	    { // use the ckvec to send
	      CkArrayIndex4D idx;
	      for(int elem=0; elem < pcid->listLFrom.size() ; elem++)
		{ 
		  calculatePairsMsg *msgfromrow=new (outsize, 8* sizeof(int)) calculatePairsMsg;
		  *(int*)CkPriorityPtr(msgfromrow) = pcid->priority;    
		  CkSetQueueing(msgfromrow, CK_QUEUEING_IFIFO);
		  msgfromrow->init(outsize, myS, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
		  idx=pcid->listLFrom[elem];
		  idx.index[3]=chunk;
#ifdef _NAN_CHECK_
		  for(int i=0;i<outsize ;i++)
		    {
		      CkAssert(finite(msgfromrow->points[i].re));
		      CkAssert(finite(msgfromrow->points[i].im));
		    }
#endif

		  pcid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).acceptPairData(msgfromrow);
		}
	    }
	  else
	    {
	      calculatePairsMsg *msgfromrow=new (outsize, 8* sizeof(int)) calculatePairsMsg;
	      *(int*)CkPriorityPtr(msgfromrow) = pcid->priority;    
	      CkSetQueueing(msgfromrow, CK_QUEUEING_IFIFO);
	      msgfromrow->init(outsize, myS, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	      if(pcid->Symmetric && myPlane==0)
		dumpMatrixDouble("pairmsg",(double *)msgfromrow->points, 1, outsize*2,myPlane,myS,0,chunk,pcid->Symmetric);
#endif
#ifdef _NAN_CHECK_
		  for(int i=0;i<outsize ;i++)
		    {
		      CkAssert(finite(msgfromrow->points[i].re));
		      CkAssert(finite(msgfromrow->points[i].im));
		    }
#endif

	      pcid->proxyLFrom[chunk].acceptPairData(msgfromrow);
	    }
	}
      
    }
  else
    {
      CkPrintf("no L proxy for [%d,%d,%d,%d,%d] !!!\n");
    }
  if(pcid->existsLNotFromproxy)
    { //symmetric
      for(int chunk=0; chunk < pcid->numChunks; chunk++)
	{
	  int chunksize=n/pcid->numChunks;
	  int outsize=chunksize;
	  if(pcid->numChunks > 1 && chunk == pcid->numChunks - 1)
	    {// last chunk  gets remainder
	      outsize += n % pcid->numChunks;
	    }
	  if( pcid->useDirectSend)
	    { // use the ckvec to send
	      CkArrayIndex4D idx;
	      for(int elem=0; elem<pcid->listLNotFrom.size();elem++)
		{ idx=pcid->listLNotFrom[elem];
		  idx.index[3]=chunk;
		  calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
		  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
		  *(int*)CkPriorityPtr(msg) = pcid->priority;    
		  msg->init(outsize, myS, false, flag_dp, &ptr[chunk * chunksize], psiV,n);   
#ifdef _NAN_CHECK_
		  for(int i=0;i<outsize ;i++)
		    {
		      CkAssert(finite(msg->points[i].re));
		      CkAssert(finite(msg->points[i].im));
		    }
#endif


		  pcid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).acceptPairData(msg);
		}
	    }
	  else
	    {
	      calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
	      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	      *(int*)CkPriorityPtr(msg) = pcid->priority;    
	      msg->init(outsize, myS, false, flag_dp, &ptr[chunk * chunksize], psiV,n);   
#ifdef _NAN_CHECK_
	      for(int i=0;i<outsize ;i++)
		{
		  CkAssert(finite(msg->points[i].re));
		  CkAssert(finite(msg->points[i].im));
		}
#endif
	      pcid->proxyLNotFrom[chunk].acceptPairData(msg);
	    }
	}
    }
#endif
#endif

}

void startPairCalcLeftRDMA(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)

{
#ifdef PC_USE_RDMA
  int NumGrains=pcid->nstates/pcid->GrainSize;
  for(int chunk=0;chunk<pcid->numChunks;chunk++)
    {
      for(int grain=0;grain<NumGrains;grain++){
#ifdef _PAIRCALC_DEBUG_RDMA_
	CkPrintf("GSP [%d,%d] addr %p RDMA put left handle %d recverNode %d senderbuf %p\n",myS, myPlane, ptr,pcid->RDMAHandlesLeft[chunk][grain].handle.handle,pcid->RDMAHandlesLeft[chunk][grain].handle.recverNode,pcid->RDMAHandlesLeft[chunk][grain].handle.senderBuf);
#endif
	if(pcid->RDMAHandlesLeft[chunk][grain].handle.handle>=0)
	  CmiDirect_put(&(pcid->RDMAHandlesLeft[chunk][grain].handle));
      }
    }
  if(pcid->Symmetric && pcid->existsLNotFromproxy)
    {
#ifdef _PAIRCALC_DEBUG_RDMA_
      CkPrintf("GSP [%d,%d] RDMA put left other\n",myS, myPlane);
#endif
      for(int chunk=0;chunk<pcid->numChunks;chunk++)
	for(int grain=0;grain<NumGrains;grain++)
	  if(pcid->RDMAHandlesRight[chunk][grain].handle.handle>=0)
	    CmiDirect_put(&(pcid->RDMAHandlesRight[chunk][grain].handle));
    }
#endif
}

void startPairCalcRightRDMA(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane)
{
#ifdef PC_USE_RDMA

  //foreach chunk, look it up to get handle and recverproc
  int NumGrains=pcid->nstates/pcid->GrainSize;
  CkAssert(!pcid->Symmetric);
  for(int chunk=0;chunk<pcid->numChunks;chunk++)
    for(int grain=0;grain<NumGrains;grain++)
    {
      
#ifdef _PAIRCALC_DEBUG_RDMA_
	CkPrintf("GSP [%d,%d] addr %p RDMA put right handle %d recverNode %d senderbuf %p\n",myS, myPlane, ptr,pcid->RDMAHandlesRight[chunk][grain].handle.handle,pcid->RDMAHandlesRight[chunk][grain].handle.recverNode,pcid->RDMAHandlesRight[chunk][grain].handle.senderBuf);
#endif
      if(pcid->RDMAHandlesRight[chunk][grain].handle.handle>=0)
	CmiDirect_put(&(pcid->RDMAHandlesRight[chunk][grain].handle));
    }
#endif
}

void dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",x+i,y+j,matrix[i*ydim+j]);
  fclose(loutfile);
}



void loadMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=0;j<ydim;j++)
	  fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i*ydim+j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}

//! NOTE: this uses the evil piny convention
void dumpMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DDouble(const char *infilename, double **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=1;j<=ydim;j++)
	  fscanf(loutfile,"%d %d %lf\n",&junk1,&junk2,&(matrix[i][j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}


//! NOTE: this uses the evil piny convention
void dumpMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=1;j<=ydim;j++)
      fprintf(loutfile,"%d %d %d\n",i,j,matrix[i][j]);
  fclose(loutfile);
}
void loadMatrix2DInt(const char *infilename, int **matrix, int xdim, int ydim,int w,int x,int y, int z, bool symmetric)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, w, x, y, z, symmetric);
  FILE *loutfile = fopen(filename, "r");
  if(loutfile!=NULL)
    {
      int junk1,junk2;
      for(int i=0;i<xdim;i++)
	for(int j=1;j<=ydim;j++)
	  fscanf(loutfile,"%d %d %d\n",&junk1,&junk2,&(matrix[i][j]));
      fclose(loutfile);
    }
  else
    {
      CkAbort(filename);
    }
}



// create multicast proxies
void makeLeftTree(PairCalcID* pcid, int myS, int myPlane){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Left symm=%d\n", pcid->Symmetric);
#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);
  int s1, s2;
  int grainSize = pcid->GrainSize;
  int numChunks =  pcid->numChunks;
  int nstates = pcid->nstates;
  int symmetric = pcid->Symmetric;
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  s1 = (myS/grainSize) * grainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
  if(!(pcid->existsLproxy||pcid->existsLNotFromproxy)){
    //create multicast proxy array section list 
    if(symmetric)
      { // we have this column send in place of a right side
	// 1 proxy for left and 1 for right
	pcid->proxyLNotFrom=new CProxySection_PairCalculator[numChunks];
	pcid->proxyLFrom=new CProxySection_PairCalculator[numChunks];
	for (int chunk = 0; chunk < numChunks; chunk++)  // new proxy for each chunk
	  {
	    CkArrayIndex4D *elems= new CkArrayIndex4D[nstates/grainSize];
	    CkArrayIndex4D *elemsfromrow= new CkArrayIndex4D[nstates/grainSize];
	    int erowcount=0;
	    int ecount=0;
	    CkArrayIndex4D idx(myPlane,0,0,chunk);
	    for(s2 = 0; s2 <= maxpcstateindex; s2 += grainSize){
	      if(s1 <= s2)
		{
		  idx.index[1]=s1;
		  idx.index[2]=s2;
		  elemsfromrow[erowcount++]=idx;
		  if(pcid->useDirectSend && chunk==0)
		    {
		      pcid->listLFrom.push_back(idx);
		    }

		}
	      else // swap s1 : s2 and toggle fromRow
		{
		  idx.index[1]=s2;
		  idx.index[2]=s1;
		  elems[ecount++]=idx;
		  if(pcid->useDirectSend && chunk==0)
		    {
		      pcid->listLNotFrom.push_back(idx);
		    }

		}
	    }
	    if(ecount)
	      {
	      
		pcid->proxyLNotFrom[chunk] = CProxySection_PairCalculator::ckNew(pairCalculatorID, elems, ecount); 
		pcid->existsLNotFromproxy=true;	  
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
		if(pcid->useComlib && _PC_COMMLIB_MULTI_)
		  {
		    ComlibAssociateProxy(&mcastInstanceCP,pcid->proxyLNotFrom[chunk]);
		  }
		else
		  {

		    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[myPlane]).ckLocalBranch();       
		    pcid->proxyLNotFrom[chunk].ckSectionDelegate(mcastGrp);
		  }
#endif
	      }
	    if(erowcount)
	      {
		pcid->proxyLFrom[chunk]  = CProxySection_PairCalculator::ckNew(pairCalculatorID, elemsfromrow, erowcount); 
		pcid->existsLproxy=true;	  
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
		if(pcid->useComlib && _PC_COMMLIB_MULTI_)
		  {
		    ComlibAssociateProxy(&mcastInstanceCP, pcid->proxyLFrom[chunk]);
		  }
		else
		  {
		    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[myPlane]).ckLocalBranch(); 
		    pcid->proxyLFrom[chunk].ckSectionDelegate(mcastGrp);

		  }
#endif
	      }
	    delete [] elemsfromrow;
	    delete [] elems;
	  }
      }
    else
      { //just left, no right
	pcid->proxyLFrom=new CProxySection_PairCalculator[numChunks];
	for (int chunk = 0; chunk < numChunks; chunk++)  // new proxy for each chunk
	  {
	    pcid->proxyLFrom[chunk] = CProxySection_PairCalculator::ckNew(pcid->Aid,  
								      myPlane, myPlane, 1,
								      s1, s1, 1,
								      0, maxpcstateindex, grainSize,
								      chunk, chunk, 1);
	    pcid->existsLproxy=true;      
	    if(pcid->useDirectSend && chunk==0)
	      {
		for(s2 = 0; s2 <= maxpcstateindex; s2 += grainSize){
		  pcid->listLFrom.push_back(CkArrayIndex4D(myPlane,s1,s2,chunk));
		}
	      }
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
	    if(pcid->useComlib && _PC_COMMLIB_MULTI_ )
	      {
		ComlibAssociateProxy(&mcastInstanceCP,pcid->proxyLFrom[chunk]);
	      }
	    else
	      {
		CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[myPlane]).ckLocalBranch(); 
		pcid->proxyLFrom[chunk].ckSectionDelegate(mcastGrp);
	      }
#endif
	  }
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

void startPairCalcRight(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane){
#ifdef PC_USE_RDMA
  startPairCalcRightRDMA(pcid,n,ptr,myS,myPlane);
#else
#ifdef _PAIRCALC_NO_MULTI_
  startPairCalcRightSlow(pcid, n, ptr, myS, myPlane);
#else
  bool flag_dp = pcid->isDoublePacked;
  if(!pcid->existsRproxy)
    {
      makeRightTree(pcid,myS,myPlane);

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
      for(int chunk=0; chunk < pcid->numChunks; chunk++)
	{
	  int chunksize=n/pcid->numChunks;
	  int outsize=chunksize;
	  if(pcid->numChunks > 1 && chunk == pcid->numChunks - 1)
	    {// last chunk gets remainder
	      outsize+=n % pcid->numChunks;
	    }
	  if(pcid->useDirectSend)
	    {
	      CkArrayIndex4D idx;
	      for(int elem=0; elem<pcid->listRNotFrom.size();elem++)
		{ idx=pcid->listRNotFrom[elem];
		  idx.index[3]=chunk;
		  calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
		  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
		  *(int*)CkPriorityPtr(msg) = pcid->priority;    
		  msg->init(outsize,myS,false,flag_dp,&ptr[chunk*chunksize],false, n);
#ifdef _NAN_CHECK_
		  for(int i=0;i<outsize ;i++)
		    {
		      CkAssert(finite(msg->points[i].re));
		      CkAssert(finite(msg->points[i].im));
		    }
#endif

		  pcid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).acceptPairData(msg);
		}
	    }
	  else
	    {
	      calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
	      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	      *(int*)CkPriorityPtr(msg) = pcid->priority;    
	      msg->init(outsize,myS,false,flag_dp,&ptr[chunk*chunksize],false, n);
#ifdef _NAN_CHECK_
		  for(int i=0;i<outsize ;i++)
		    {
		      CkAssert(finite(msg->points[i].re));
		      CkAssert(finite(msg->points[i].im));
		    }
#endif

	      pcid->proxyRNotFrom[chunk].acceptPairData(msg);
	    }
	}
    }
  else
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("Warning! No Right proxy ! \n");
#endif
    }
#endif //_NO_MULTI
#endif
}

void makeRightTree(PairCalcID* pcid, int myS, int myPlane){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Right symm=%d\n", pcid->Symmetric);
#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);

  int s2, c;
  int grainSize = pcid->GrainSize;
  int numChunks =  pcid->numChunks;
  int nstates = pcid->nstates;
  bool symmetric = pcid->Symmetric;
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  CkAssert(symmetric == false);
  s2 = (myS/grainSize) * grainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
  //create multicast proxy list 
  if(!pcid->existsRproxy)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("initializing R multicast proxy in %d %d \n",myPlane,s2);
#endif
      pcid->proxyRNotFrom=new CProxySection_PairCalculator[numChunks];
      for (c = 0; c < numChunks; c++)  // new proxy for each chunk
	{
	  pcid->proxyRNotFrom[c] = 
	    CProxySection_PairCalculator::ckNew(pcid->Aid,  
						myPlane, myPlane, 1,
						0, maxpcstateindex, grainSize,
						s2, s2, 1,
						c, c, 1);
	  if(pcid->useDirectSend && c==0)
	    {
	      for(int s1 = 0; s1 <= maxpcstateindex; s1 += grainSize){
		pcid->listRNotFrom.push_back(CkArrayIndex4D(myPlane,s1,s2,c));
	      }
	    }
	  pcid->existsRproxy=true;      
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
	  if(pcid->useComlib && _PC_COMMLIB_MULTI_)
	    {
	      ComlibAssociateProxy(&mcastInstanceCP, pcid->proxyRNotFrom[c]);
	    }
	  else
	    {
	      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[myPlane]).ckLocalBranch(); 
	      pcid->proxyRNotFrom[c].ckSectionDelegate(mcastGrp);
	    }
#endif
	}
    }
}




void finishPairCalcSection(int n, double *ptr, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority) {
  finishPairCalcSection2(n, ptr, NULL, pcid, orthoX, orthoY, actionType, priority);
}


/* This version uses a section multicast to only send the part of the matrix needed by each section */
void finishPairCalcSection2(int n, double *ptr1, double *ptr2, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority) {
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Finish Mcast 2\n");
#endif
  //NOTE need some configuration adherence check here!
  /*
  if(pcid->Symmetric)
    ComlibAssociateProxy(&mcastInstanceCP,pcid->proxySym);
  else
    ComlibAssociateProxy(&mcastInstanceACP,pcid->proxyAsym);
  */

  if(ptr2==NULL){
#ifdef _NAN_CHECK_
    for(int i=0;i<n ;i++)
      {
	if(pcid->Symmetric)  // just so we can discern in the abort
	  CkAssert(finite(ptr1[i]));
	else
	  CkAssert(finite(ptr1[i]));
      }
#endif

    multiplyResultMsg *omsg;

    if(priority>0)
      {
	omsg=new ( n,0,8*sizeof(int) ) multiplyResultMsg;
	*(int*)CkPriorityPtr(omsg) = priority;    
	CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
      }
    else
      {
	omsg=new ( n,0 ) multiplyResultMsg;
      }
    omsg->init1(n, ptr1, orthoX, orthoY, actionType);
#ifdef _NAN_CHECK_
    for(int i=0;i<n ;i++)
      {
	CkAssert(finite(omsg->matrix1[i]));
      }
#endif
    if(pcid->Symmetric)
      pcid->proxySym.multiplyResult(omsg);
    else
      pcid->proxyAsym.multiplyResult(omsg);
  }
  else {
    multiplyResultMsg *omsg;
    if(priority>0)
      {
	omsg=new ( n,n, 8*sizeof(int) ) multiplyResultMsg;
	*(int*)CkPriorityPtr(omsg) = priority;    
	CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
      }
    else
      {
	omsg=new ( n,n ) multiplyResultMsg;
      }
    omsg->init(n, n, ptr1, ptr2, orthoX, orthoY, actionType);
  if(pcid->Symmetric)
    pcid->proxySym.multiplyResult(omsg);
  else
    pcid->proxyAsym.multiplyResult(omsg);
  }
}


/* Send orthoT now that we have it so it will be ready when Asymm
   starts this iteration.  */
void sendMatrix(int n, double *ptr1,PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority) {
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc SendMatrix (orthoT)\n");
#endif

    multiplyResultMsg *omsg;

    if(priority>0)
      {
	omsg=new ( n,0,8*sizeof(int) ) multiplyResultMsg;
	*(int*)CkPriorityPtr(omsg) = priority;    
	CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
      }
    else
      {
	omsg=new ( n,0 ) multiplyResultMsg;
      }
    omsg->init1(n, ptr1, orthoX, orthoY, actionType);
#ifdef _NAN_CHECK_
    for(int i=0;i<n ;i++)
      {
	CkAssert(finite(omsg->matrix1[i]));
      }
#endif
    pcid->proxyAsym.acceptOrthoT(omsg);
}


void initPairCalcRDMA(PairCalcID *pid, int sender, int totalsize,int myPlane)
{
  //for each pid, send the notify for left and right
  int senderProc=CkMyPe();
  int chunksize =  totalsize / pid->numChunks;
  int outsize = chunksize;
  if(!(pid->existsLproxy||pid->existsLNotFromproxy)){
    makeLeftTree(pid,sender,myPlane);
    CkArrayID pairCalculatorID = (CkArrayID)pid->Aid; 
    pid->cproxy= CProxy_PairCalculator(pairCalculatorID);
  }
  if(!pid->Symmetric && !pid->existsRproxy)
    {
      makeRightTree(pid,sender,myPlane);
      
    }
  int numGrains=pid->nstates/pid->GrainSize;
  int numChunks=pid->numChunks;
  CkAssert(numChunks>0);
  CkAssert(numGrains>0);
  pid->RDMAHandlesLeft=new RDMAHandle*[numChunks];
  pid->RDMAHandlesRight=new RDMAHandle*[numChunks];
  for(int chunk=0;chunk<numChunks;chunk++){
    pid->RDMAHandlesLeft[chunk]=new RDMAHandle[numGrains];
    pid->RDMAHandlesRight[chunk]=new RDMAHandle[numGrains];
  }
#ifdef _PAIRCALC_DEBUG_RDMA_
  CkPrintf("init pair RDMA sender %d size %d chunks %d\n",sender,totalsize,pid->numChunks);
#endif
  for (int chunk = 0; chunk < numChunks; chunk++)  // new proxy for
					   // each chunk
    {
      //use proxy to send
      
      if((numChunks > 1) && (chunk == (numChunks - 1)))
	{// last chunk gets remainder
	  outsize= chunksize + (totalsize % numChunks);
	}
      // send to left regardless
      if(pid->existsLproxy)
	{
	  if( pid->useDirectSend)
	    { // use the ckvec to send
	      CkArrayIndex4D idx;
	      for(int elem=0; elem < pid->listLFrom.size() ; elem++)
		{
		  idx=pid->listLFrom[elem];
		  idx.index[3]=chunk;
		  pid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).receiveRDMASenderNotify(senderProc, sender, true, outsize,totalsize);
		}
	    }
	  else{
	    pid->proxyLFrom[chunk].receiveRDMASenderNotify(senderProc, sender, true, outsize,totalsize);
	  }
	}
      if(pid->Symmetric)
	{
	  if(pid->existsLNotFromproxy)
	    {
	      if( pid->useDirectSend)
		{ // use the ckvec to send
		  CkArrayIndex4D idx;
		  for(int elem=0; elem < pid->listLNotFrom.size() ; elem++)
		    {
		      idx=pid->listLNotFrom[elem];
		      idx.index[3]=chunk;
		      pid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).receiveRDMASenderNotify(senderProc, sender, false, outsize,totalsize);
		    }
		}
	      else
		{

		  pid->proxyLNotFrom[chunk].receiveRDMASenderNotify(senderProc, sender, false, outsize,totalsize);
		}
	    }
	}
      else
	{
	  if(pid->existsRproxy)
	    {
	      if( pid->useDirectSend)
		{ // use the ckvec to send
		  CkArrayIndex4D idx;
		  for(int elem=0; elem < pid->listRNotFrom.size() ; elem++)
		    {
		      idx=pid->listRNotFrom[elem];
		      idx.index[3]=chunk;
		      pid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).receiveRDMASenderNotify(senderProc, sender, false, outsize,totalsize);
		    }
		}
	      else{
		
		pid->proxyRNotFrom[chunk].receiveRDMASenderNotify(senderProc, sender, false, outsize, totalsize);
	      }
	    }
	}
    }
}
#ifdef ROTATE_LIST
bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart)
{
  if(newstart>numelems)
    return(false);
  CkArrayIndexMax  swap;
  int where=newstart;  
  for(int i=0;i<newstart;i++,where++)
    {
      if(where>=numelems)
	where=0;
      swap=elems[i];
      elems[i]=elems[where];
      elems[where]=swap;
    }
  return(true);
}

#else

bool reorder_elem_list(CkArrayIndexMax *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}

bool reorder_elem_list_4D(CkArrayIndex4D *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}

bool reorder_elem_list_max(CkArrayIndexMax *elems, int numelems, int newstart)
{
  //  CkPrintf("reordering list of %d elems\n", numelems);
  std::random_shuffle(elems,elems+numelems);
  return(true);
}
#endif
