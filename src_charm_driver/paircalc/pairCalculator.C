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
 * Expected usage begins with createPairCalculator(*,PairCalcID *,*)
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
#include "InputDataHandler.h"

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

  CkArrayOptions paircalcOpts,handlerOpts;
  CProxy_PairCalculator pairCalculatorProxy;
  CProxy_InputDataHandler<leftCollatorType,rightCollatorType> inputHandlerProxy;
  redtypes cpreduce=section;

#ifdef CONVERSE_VERSION_ELAN
  bool machreduce=(s/grainSize * numZ* numChunks>=CkNumNodes()) ? true: false;
#else
  bool machreduce=false;
#endif

  if(machreduce)
    cpreduce=machine;
  
  // If a chare mapping is not available, create an empty array
  if(!mapid) 
  {
    pairCalculatorProxy = CProxy_PairCalculator::ckNew();
  }
  // else, create the array with element locations as specified by the map 
  else 
  {
    paircalcOpts.setMap(*mapid);
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(inputHandlerProxy, sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, gSpaceSum, gpriority, phantomSym, useBWBarrier,
						       gemmSplitFWk, gemmSplitFWm,
						       gemmSplitBW,expectOrthoT, instance,
						       paircalcOpts);
  }

#ifdef _PAIRCALC_CREATE_DEBUG_
	CkPrintf("createPairCalculator: Creating empty inputHandler chare array for asymm(0)/symm(1) loop PCs:%d \n",sym);
#endif 
  /// Create an empty input handler chare array that will accept all incoming messages from GSpace
  handlerOpts.bindTo(pairCalculatorProxy);
  inputHandlerProxy = CProxy_InputDataHandler<leftCollatorType,rightCollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);

  int proc = 0;
  // Initialize the PairCalcID instance
  pcid->Init(pairCalculatorProxy.ckGetArrayID(), grainSize, numChunks, s, sym, comlib_flag, flag_dp, conserveMemory, lbpaircalc,  priority, useDirectSend);
  pcid->orthomCastGrpId=orthomCastGrpId;
  pcid->orthoRedGrpId=orthoRedGrpId;
  pcid->cproxy=pairCalculatorProxy;
  pcid->mCastGrpId=mCastGrpId;
  
  // Setup the appropriate multicast strategy
#ifdef CMK_BLUEGENEL
  //  CharmStrategy *multistrat = new RectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
#else
  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
#endif
  
  // 
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  // 
  if(sym)
    mcastInstanceCP=ComlibRegister(multistrat);
  else
    mcastInstanceACP=ComlibRegister(multistrat);
  CkAssert(mapid);
  // If the symmetric loop PC instances are being created
  if(sym)
	for(int numX = 0; numX < numZ; numX ++)
	{
	  for (int s1 = 0; s1 <= maxpcstateindex; s1 += grainSize) 
      {
      	// If phantomSym is turned on
		int s2start=(phantomSym) ? 0 : s1;
		for (int s2 = s2start; s2 <= maxpcstateindex; s2 += grainSize) 
		{
			for (int c = 0; c < numChunks; c++) 
			{
				if(mapid) 
				{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym);
#endif
					pairCalculatorProxy(z[numX],s1,s2,c).
						insert(inputHandlerProxy, sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, gSpaceSum, gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance );
				}
				else
				{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym);
#endif
					pairCalculatorProxy(z[numX],s1,s2,c).
						insert(inputHandlerProxy, sym, grainSize, s, numChunks, cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, gSpaceSum, gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance, proc);
					proc++;
					if (proc >= CkNumPes()) proc = 0;
				}
			}
		}
	  }
	}
  // else, if the asymmetric loop PC instances are being created
  else
  {
	for(int numX = 0; numX < numZ; numX ++)
	{
		for (int s1 = 0; s1 <= maxpcstateindex; s1 += grainSize)
		{
			for (int s2 = 0; s2 <= maxpcstateindex; s2 += grainSize)
			{
				for (int c = 0; c < numChunks; c++)
				{
					if(mapid)
					{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym);
#endif
					pairCalculatorProxy(z[numX],s1,s2,c).
						insert(inputHandlerProxy, sym, grainSize, s, numChunks, cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc,  cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, gSpaceSum,  gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance);
					}
					else
					{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym);
#endif
						pairCalculatorProxy(z[numX],s1,s2,c).
							insert(inputHandlerProxy, sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, rdma_ep, conserveMemory, lbpaircalc,   cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, gSpaceSum, gpriority, phantomSym, useBWBarrier, gemmSplitFWk, gemmSplitFWm, gemmSplitBW, expectOrthoT, instance,proc);
						proc++;
						if (proc >= CkNumPes()) proc = 0;
					}
				}
			}
		}
	}
  }
  /// Notify the runtime that we're done inserting all the PC elements
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



// create multicast proxies
void makeLeftTree(PairCalcID* pcid, int myS, int myPlane)
{
	#ifdef _PAIRCALC_DEBUG_
		CkPrintf("GSpace [%d,%d] Making symm(%d) PC array section to receive left data \n", myS, myPlane, pcid->Symmetric);
	#endif
	CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid;
	CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);
	int s1, s2,sColMin;
	int grainSize = pcid->GrainSize;
	int numChunks =  pcid->numChunks;
	int nstates = pcid->nstates;
	int symmetric = pcid->Symmetric;
	int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
	/// Find the row index of the PC chare that handles this state
	s1 = (myS/grainSize) * grainSize;
	s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
	/// If the PC is a symmetric instance, then include only the post-diagonal chares on the row s1  
	/// else, include all the PC chares on row s1
	symmetric ? sColMin = s1 : sColMin = 0;

    /// Allocate one section proxy for each chunk
	pcid->sectionGettingLeft=new CProxySection_PairCalculator[numChunks];
	/// Build an array section for each chunk 
	for (int chunk = 0; chunk < numChunks; chunk++)
	{
		pcid->sectionGettingLeft[chunk] = CProxySection_PairCalculator::ckNew(pcid->Aid,
								      myPlane, myPlane, 1,
								      s1, s1, 1,
								      sColMin, maxpcstateindex, grainSize,
								      chunk, chunk, 1);
		pcid->existsLproxy=true;
		if(pcid->useDirectSend && chunk==0)
			for(s2 = sColMin; s2 <= maxpcstateindex; s2 += grainSize)
				pcid->listGettingLeft.push_back(CkArrayIndex4D(myPlane,s1,s2,chunk));
				
		#ifndef _PAIRCALC_DO_NOT_DELEGATE_
			if(pcid->useComlib && _PC_COMMLIB_MULTI_ )
				ComlibAssociateProxy(&mcastInstanceCP,pcid->sectionGettingLeft[chunk]);
			else
			{
				CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[myPlane]).ckLocalBranch();
				pcid->sectionGettingLeft[chunk].ckSectionDelegate(mcastGrp);
			}
		#endif
	  }
}



void makeRightTree(PairCalcID* pcid, int myS, int myPlane){
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("GSpace [%d,%d] Making symm(%d) PC array section to receive right data \n", myS, myPlane, pcid->Symmetric);
#endif
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid;
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);

  int s2, sRowMax, c;
  int grainSize = pcid->GrainSize;
  int numChunks =  pcid->numChunks;
  int nstates = pcid->nstates;
  bool symmetric = pcid->Symmetric;
  int maxpcstateindex=(pcid->nstates/pcid->GrainSize-1)*pcid->GrainSize;
  CkAssert(symmetric == false);
  s2 = (myS/grainSize) * grainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
  /// If the PC is a symmetric instance, then include only the pre-diagonal chares on the column s2 else, include all the PC chares on column s2
  if (symmetric)
    sRowMax = s2-grainSize;
  else
    sRowMax = maxpcstateindex;
  //create multicast proxy list
    	#ifdef _PAIRCALC_DEBUG_
    		CkPrintf("initializing R multicast proxy in %d %d \n",myPlane,s2);
    	#endif
      pcid->sectionGettingRight=new CProxySection_PairCalculator[numChunks];
      for (c = 0; c < numChunks; c++)  // new proxy for each chunk
	{
		// Accomodate the boundary case:  The first PC chare on the top left [0,0] of the array shouldnt receive a right data. So leave its right data proxy empty
		if (sRowMax >=0)
		{
			pcid->sectionGettingRight[c] = CProxySection_PairCalculator::ckNew(pcid->Aid,
						myPlane, myPlane, 1,
						0, sRowMax, grainSize,
						s2, s2, 1,
						c, c, 1);
	  if(pcid->useDirectSend && c==0)
	      for(int s1 = 0; s1 <= sRowMax; s1 += grainSize)
                  pcid->listGettingRight.push_back(CkArrayIndex4D(myPlane,s1,s2,c));
		}
                  
        pcid->existsRproxy=true;
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
	  if(pcid->useComlib && _PC_COMMLIB_MULTI_)
	      ComlibAssociateProxy(&mcastInstanceCP, pcid->sectionGettingRight[c]);
	  else
	    {
	      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[myPlane]).ckLocalBranch();
	      pcid->sectionGettingRight[c].ckSectionDelegate(mcastGrp);
	    }
#endif
	}
}



/** For symmetric instances, this deposits the data as 
 *  	- Left matrix block to the post-diagonal row of PCs that correspond to state myS
 *  	- Right matrix block to the pre-diagonal column of PCs that correspond to state myS
 * 
 * For asymmetric instances, this deposits the data as left matrix data to the whole row of PCs that correspond to state myS
 */
void startPairCalcLeft(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
    #ifdef PC_USE_RDMA
        // RDMA doesn't support PSIV shenanighans
        if(psiV)
            CkAbort("You must #undef PC_USE_RDMA in ckPairCalculator.h for dynamics. PSIV is broken by RDMA.\n It can be fixed, but won't be done until someone needs production support.");
        else
            startPairCalcLeftRDMA(pcid, n, ptr, myS, myPlane,psiV);
    #else
        /// If multicasting is disabled, use a slow point-to-point send version
        #ifdef _PAIRCALC_NO_MULTI_
            startPairCalcLeftSlow(pcid, n, ptr, myS, myPlane);
        #else
            bool flag_dp = pcid->isDoublePacked;
            /// If a destination array section doesnt exist, build one
            if(!pcid->existsLproxy)
            {
                makeLeftTree(pcid,myS,myPlane);
                CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid;
                pcid->cproxy= CProxy_PairCalculator(pairCalculatorID);
            }
            /// If a left matrix destination section exists, send the data as the left matrix block
            if(pcid->existsLproxy)
            {
                int numChunks=pcid->numChunks;
                int chunksize =  n / numChunks;
                int outsize = chunksize;

                for(int chunk=0; chunk < numChunks ; chunk++)
                {
                    // last chunk gets remainder
                    if((numChunks > 1) && (chunk == (numChunks - 1)))
                        outsize= chunksize + (n % numChunks);
                    #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
                    if(pcid->Symmetric && myPlane==0)
                        dumpMatrixDouble("gspPts",(double *)ptr, 1, n*2,myPlane,myS,0,chunk,pcid->Symmetric);
                    CkPrintf("L [%d,%d,%d,%d,%d] chunk %d chunksize %d outsize %d for numpoint %d offset will be %d %.12g\n",myPlane,myS, myS, chunk,pcid->Symmetric, chunk,chunksize,outsize,n,chunk*chunksize,ptr[chunk*chunksize].re);
                    #endif
                    // If sending directly, use the vector of target PC chares
                    if( pcid->useDirectSend)
                    {
                        CkArrayIndex4D idx;
                        for(int elem=0; elem < pcid->listGettingLeft.size() ; elem++)
                        {
                            paircalcInputMsg *msgfromrow=new (outsize, 8* sizeof(int)) paircalcInputMsg;
                            *(int*)CkPriorityPtr(msgfromrow) = pcid->priority;
                            CkSetQueueing(msgfromrow, CK_QUEUEING_IFIFO);
                            msgfromrow->init(outsize, myS, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                            idx=pcid->listGettingLeft[elem];
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
                    // else, use a typical multicast to the destination section
                    else
                    {
                        paircalcInputMsg *msgfromrow=new (outsize, 8* sizeof(int)) paircalcInputMsg;
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
                        pcid->sectionGettingLeft[chunk].acceptPairData(msgfromrow);
                    }
                }
            }
            /// else, if the destination section doesnt exist even after attempting to create one
            else
                CkPrintf("GSpace [%d,%d] No destination symm(%d) PC array section to send left block data [%d,%d,%d,%d,%d] !!!\n",myS,myPlane,pcid->Symmetric);

            /// temporary
            if(pcid->Symmetric)
            	startPairCalcRight(pcid, n, ptr, myS, myPlane, psiV);
        #endif
    #endif
}



/** Deposits the data as right matrix block to the PC chare section associated with state myS
 */
void startPairCalcRight(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
    #ifdef PC_USE_RDMA
        // RDMA doesn't support PSIV shenanighans
        if(psiV)
            CkAbort("You must #undef PC_USE_RDMA in ckPairCalculator.h for dynamics. PSIV is broken by RDMA.\n It can be fixed, but won't be done until someone needs production support.");
        else
            startPairCalcRightRDMA(pcid, n, ptr, myS, myPlane,psiV);
    #else
        /// If multicasting is disabled, use a slow point-to-point send version
        #ifdef _PAIRCALC_NO_MULTI_
            startPairCalcRightSlow(pcid, n, ptr, myS, myPlane);
        #else
            bool flag_dp = pcid->isDoublePacked;
            /// If a destination array section doesnt exist, build one
            if(!pcid->existsRproxy)
                makeRightTree(pcid,myS,myPlane);
            /// If a right matrix destination section exists, send the data as the left matrix block
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
                    /// last chunk gets remainder
                    if(pcid->numChunks > 1 && chunk == pcid->numChunks - 1)
                        outsize+=n % pcid->numChunks;
                    if(pcid->useDirectSend)
                    {
                        CkArrayIndex4D idx;
                        for(int elem=0; elem<pcid->listGettingRight.size();elem++)
                        { 
                            idx=pcid->listGettingRight[elem];
                            idx.index[3]=chunk;
                            paircalcInputMsg *msg= new ( outsize,8*sizeof(int) ) paircalcInputMsg;
                            CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                            *(int*)CkPriorityPtr(msg) = pcid->priority;
                            msg->init(outsize,myS,false,flag_dp,&ptr[chunk*chunksize],psiV, n);
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
                        paircalcInputMsg *msg= new ( outsize,8*sizeof(int) ) paircalcInputMsg;
                        CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                        *(int*)CkPriorityPtr(msg) = pcid->priority;
                        msg->init(outsize,myS,false,flag_dp,&ptr[chunk*chunksize],psiV, n);
                        #ifdef _NAN_CHECK_
                        for(int i=0;i<outsize ;i++)
                        {
                            CkAssert(finite(msg->points[i].re));
                            CkAssert(finite(msg->points[i].im));
                        }
                        #endif
                        pcid->sectionGettingRight[chunk].acceptPairData(msg);
                    }
                }
            }
            /// else, if the destination section doesnt exist even after attempting to create one
            else
                CkPrintf("GSpace [%d,%d] No destination symm(%d) PC array section to send right block data [%d,%d,%d,%d,%d] !!!\n",myS,myPlane,pcid->Symmetric);
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
  if(pcid->Symmetric && pcid->existsRproxy)
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
  if(!pid->existsLproxy){
    makeLeftTree(pid,sender,myPlane);
    CkArrayID pairCalculatorID = (CkArrayID)pid->Aid;
    pid->cproxy= CProxy_PairCalculator(pairCalculatorID);
  }
  if(!pid->existsRproxy)
      makeRightTree(pid,sender,myPlane);

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
	      for(int elem=0; elem < pid->listGettingLeft.size() ; elem++)
		{
                  idx=pid->listGettingLeft[elem];
		  idx.index[3]=chunk;
		  pid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).receiveRDMASenderNotify(senderProc, sender, true, outsize,totalsize);
		}
	    }
	  else{
	    pid->sectionGettingLeft[chunk].receiveRDMASenderNotify(senderProc, sender, true, outsize,totalsize);
	  }
	}
      if(pid->Symmetric)
	{
	  if(pid->existsRproxy)
	    {
	      if( pid->useDirectSend)
		{ // use the ckvec to send
		  CkArrayIndex4D idx;
		  for(int elem=0; elem < pid->listGettingRight.size() ; elem++)
		    {
		      idx=pid->listGettingRight[elem];
		      idx.index[3]=chunk;
		      pid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).receiveRDMASenderNotify(senderProc, sender, false, outsize,totalsize);
		    }
		}
	      else
		{

		  pid->sectionGettingRight[chunk].receiveRDMASenderNotify(senderProc, sender, false, outsize,totalsize);
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
		  for(int elem=0; elem < pid->listGettingRight.size() ; elem++)
		    {
		      idx=pid->listGettingRight[elem];
		      idx.index[3]=chunk;
		      pid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).receiveRDMASenderNotify(senderProc, sender, false, outsize,totalsize);
		    }
		}
	      else{

		pid->sectionGettingRight[chunk].receiveRDMASenderNotify(senderProc, sender, false, outsize, totalsize);
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

