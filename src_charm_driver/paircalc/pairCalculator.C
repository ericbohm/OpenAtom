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
extern ComlibInstanceHandle mcastInstanceCP;
extern ComlibInstanceHandle gAsymInstance;
extern ComlibInstanceHandle gSymInstance;


void createPairCalculator(bool sym, int s, int grainSize, int numZ, int* z, 
			  CkCallback cb,  PairCalcID* pcid, int cb_ep, 
			  int cb_ep_tol, 
			  CkArrayID cb_aid, int comlib_flag, CkGroupID *mapid,
			  int flag_dp, bool conserveMemory, bool lbpaircalc, 
			  int priority, CkVec <CkGroupID> mCastGrpId, int numChunks, int orthoGrainSize, int useEtoM, bool collectTiles, bool streamBWout, bool delayBWSend, int streamFW) {

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
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, options);
  }

  int proc = 0;

  pcid->Init(pairCalculatorProxy.ckGetArrayID(), grainSize, numChunks, s, sym, comlib_flag, flag_dp, conserveMemory, lbpaircalc,  priority, useEtoM);
  pcid->cproxy=pairCalculatorProxy;
  pcid->mCastGrpId=mCastGrpId;
  CharmStrategy *multistrat = new DirectMulticastStrategy(pairCalculatorProxy.ckGetArrayID());
  if(sym)// cheap hack to only do this once 
    mcastInstanceCP=ComlibRegister(multistrat);

  if(sym)
    for(int numX = 0; numX < numZ; numX ++){
      for (int s1 = 0; s1 < s; s1 += grainSize) {
	for (int s2 = s1; s2 < s; s2 += grainSize) {
	  for (int c = 0; c < numChunks; c++) {
	    if(mapid) {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
	      pairCalculatorProxy(z[numX],s1,s2,c).
		insert(sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW );
	    }
	    else
	      {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, numChunks, cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc, cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW,  proc);
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
	for (int s1 = 0; s1 < s; s1 += grainSize) {
	  for (int s2 = 0; s2 < s; s2 += grainSize) {
	    for (int c = 0; c < numChunks; c++) {
	      if(mapid) {
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, numChunks, cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc,  cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW );
	      }
	      else{
#ifdef _PAIRCALC_CREATE_DEBUG_
	      CkPrintf("inserting [%d %d %d %d %d]\n",z[numX],s1,s2,c,sym); 
#endif
		pairCalculatorProxy(z[numX],s1,s2,c).
		  insert(sym, grainSize, s, numChunks,  cb, cb_aid, cb_ep, cb_ep_tol, conserveMemory, lbpaircalc,   cpreduce, orthoGrainSize, collectTiles, streamBWout, delayBWSend, streamFW, proc);
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
  int s2=state/pcid->GrainSize*pcid->GrainSize;
  int nstates=pcid->nstates;
  int GrainSize=pcid->GrainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  
									       plane, plane, 1,
									       0, nstates-GrainSize, GrainSize,
									       s2, s2, 1,
									       chunk, chunk,1);
  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart=state%GrainSize + chunk;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list( sid._elems, sid._nElems, newListStart);
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
  int nstates=pcid->nstates;
  // all nondiagonal elements 
  // so we'll have to make this the tedious explicit way
  CkArrayIndexMax *elems= new CkArrayIndexMax[nstates/GrainSize-1];
  int ecount=0;
  for(int s2 =0; s2<nstates; s2+=GrainSize)
    {
      if(s1!=s2)
	{
	  CkArrayIndex4D idx4d(plane,s1,s2,chunk);
	  elems[ecount++]=idx4d;
	}
    }
  int newListStart=state%GrainSize + chunk;
  if(newListStart> ecount)
    newListStart= newListStart % ecount;
  bool order=reorder_elem_list( elems, ecount, newListStart);
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
  int s2=state/pcid->GrainSize*pcid->GrainSize; //row
  int GrainSize=pcid->GrainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  
									       plane, plane, 1,
									       0, s2, GrainSize,
									       s2, s2, 1,
									       chunk, chunk, 1);
  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart=state%GrainSize + chunk;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list( sid._elems, sid._nElems, newListStart);
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
  int s1=state/pcid->GrainSize*pcid->GrainSize; //column
  int nstates=pcid->nstates;
  int GrainSize=pcid->GrainSize;
  CkAssert(s1+GrainSize<nstates);
  CProxySection_PairCalculator sectProxy = 
      CProxySection_PairCalculator::ckNew(pcid->Aid,  
					  plane, plane, 1,
					  s1, s1, 1,
					  s1+GrainSize, nstates-GrainSize, GrainSize,
					  chunk, chunk, 1);

  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart=state%GrainSize + chunk;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list( sid._elems, sid._nElems, newListStart);
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
CProxySection_PairCalculator initOneRedSect(int numZ, int* z, int numChunks,  PairCalcID* pcid, CkCallback cb, int s1, int s2, int orthoX, int orthoY, int orthoGrainSize)
{
  int ecount=0;
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("initGred for s1 %d s2 %d ortho %d %d sym %d\n",s1,s2,orthoX, orthoY,pcid->Symmetric);
#endif
  CkArrayIndexMax *elems= new CkArrayIndexMax[numZ*numChunks];
  //add chunk loop
  for(int chunk = numChunks-1; chunk >=0; chunk--){
    for(int numX = numZ-1; numX >=0; numX--){
      CkArrayIndex4D idx4d(z[numX],s1,s2,chunk);
      elems[ecount++]=idx4d;
    }
  }
  int numOrtho=pcid->GrainSize/orthoGrainSize;
  int orthoIndexX=(orthoX*orthoGrainSize-s1)/orthoGrainSize;
  int orthoIndexY=(orthoY*orthoGrainSize-s2)/orthoGrainSize;
  int orthoIndex=orthoIndexX*numOrtho+orthoIndexY;

  int newListStart=orthoIndex;
  if(newListStart> ecount)
    newListStart= newListStart % ecount;
  bool order=reorder_elem_list( elems, ecount, newListStart);
  CkAssert(order);

  // now that we have the section, make the proxy and do delegation
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcid->Aid,  elems, ecount); 
  delete [] elems;

  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId[0]).ckLocalBranch();       
  sectProxy.ckSectionDelegate(mcastGrp);

  // send the message to initialize it with the callback and groupid
  setGredProxy(&sectProxy, pcid->mCastGrpId[0], cb, false, CkCallback(CkCallback::ignore), orthoX, orthoY);
  return sectProxy;
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
#ifdef _PAIRCALC_NO_MULTI_
  startPairCalcLeftSlow(pcid, n, ptr, myS, myPlane);
#else
  bool flag_dp = pcid->isDoublePacked;
  if(!(pcid->existsLproxy||pcid->existsLNotFromproxy)){
    makeLeftTree(pcid,myS,myPlane);
    CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
    pcid->cproxy= CProxy_PairCalculator(pairCalculatorID);
    if(pcid->useEtoM)
      if(pcid->Symmetric)
	ComlibAssociateProxy(&gSymInstance,pcid->cproxy);
      else
	ComlibAssociateProxy(&gAsymInstance,pcid->cproxy);
  }
  //use proxy to send
  if(pcid->useEtoM)
    if(pcid->Symmetric)
      gSymInstance.beginIteration();
    else
      gAsymInstance.beginIteration();

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
	  if(symmetric && myPlane==0)
	    dumpMatrixDouble("gspPts",(double *)ptr, 1, n*2,myPlane,myS,myS,chunk,symmetric);
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	  CkPrintf("L [%d,%d,%d,%d,%d] chunk %d chunksize %d outsize %d for numpoint %d offset will be %d %.12g\n",myPlane,myS, myS, chunk,symmetric, chunk,chunksize,outsize,n,chunk*chunksize,ptr[chunk*chunksize].re);
#endif

	  if(pcid->useEtoM)
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
	      if(symmetric && myPlane==0)
		dumpMatrixDouble("pairmsg",(double *)msgfromrow->points, 1, outsize*2,myPlane,myS,myS,chunk,symmetric);
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
	  if(pcid->useEtoM)
	    { // use the ckvec to send
	      CkArrayIndex4D idx;
	      for(int elem=0; elem<pcid->listLNotFrom.size();elem++)
		{ idx=pcid->listLNotFrom[elem];
		  idx.index[3]=chunk;
		  calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
		  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
		  *(int*)CkPriorityPtr(msg) = pcid->priority;    
		  msg->init(outsize, myS, false, flag_dp, &ptr[chunk * chunksize], psiV,n);   

		  pcid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).acceptPairData(msg);
		}
	    }
	  else
	    {
	      calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
	      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	      *(int*)CkPriorityPtr(msg) = pcid->priority;    
	      msg->init(outsize, myS, false, flag_dp, &ptr[chunk * chunksize], psiV,n);   
	      pcid->proxyLNotFrom[chunk].acceptPairData(msg);
	    }
	}
    }


  if(pcid->useComlib && _PC_COMMLIB_MULTI_) {
  }  
#endif

  if(pcid->useEtoM)
    if(pcid->Symmetric)
      gSymInstance.endIteration();
  //    else
  //      gAsymInstance.endIteration(); in startRight


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
      fprintf(loutfile,"%d %d %.12g\n",i,j,matrix[i*ydim+j]);
  fclose(loutfile);
}

// create multicast proxies
void makeLeftTree(PairCalcID* pcid, int myS, int myPlane){
  CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
  CProxy_PairCalculator pairCalculatorProxy(pairCalculatorID);
  int s1, s2;
  int grainSize = pcid->GrainSize;
  int numChunks =  pcid->numChunks;
  int nstates = pcid->nstates;
  int symmetric = pcid->Symmetric;

  s1 = (myS/grainSize) * grainSize;
  if(!(pcid->existsLproxy||pcid->existsLNotFromproxy)){
    //create multicast proxy array section list 
    if(symmetric)
      { // we have this column send in place of a right side
	// 1 proxy for left and 1 for right
	pcid->proxyLNotFrom=new CProxySection_PairCalculator[numChunks];
	pcid->proxyLFrom=new CProxySection_PairCalculator[numChunks];
	for (int chunk = 0; chunk < numChunks; chunk++)  // new proxy for each chunk
	  {
	    CkArrayIndexMax *elems= new CkArrayIndexMax[nstates/grainSize];
	    CkArrayIndexMax *elemsfromrow= new CkArrayIndexMax[nstates/grainSize];
	    int erowcount=0;
	    int ecount=0;
	    CkArrayIndex4D idx(myPlane,0,0,chunk);
	    for(s2 = 0; s2 < nstates; s2 += grainSize){
	      if(s1 <= s2)
		{
		  idx.index[1]=s1;
		  idx.index[2]=s2;
		  elemsfromrow[erowcount++]=idx;
		  if(pcid->useEtoM && chunk==0)
		    {
		      pcid->listLFrom.push_back(idx);
		    }

		}
	      else // swap s1 : s2 and toggle fromRow
		{
		  idx.index[1]=s2;
		  idx.index[2]=s1;
		  elems[ecount++]=idx;
		  if(pcid->useEtoM && chunk==0)
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
								      0, nstates-grainSize, grainSize,
								      chunk, chunk, 1);
	    pcid->existsLproxy=true;      
	    if(pcid->useEtoM && chunk==0)
	      {
		for(s2 = 0; s2 < nstates; s2 += grainSize){
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
#ifdef _PAIRCALC_NO_MULTI_
  startPairCalcRightSlow(pcid, n, ptr, myS, myPlane);
#else
  bool flag_dp = pcid->isDoublePacked;
  if(!pcid->existsRproxy)
    {
      makeRightTree(pcid,myS,myPlane);
      /*    CkArrayID pairCalculatorID = (CkArrayID)pcid->Aid; 
    pcid->cproxy= CProxy_PairCalculator(pairCalculatorID);
    if(pcid->useEtoM)
      if(pcid->Symmetric)
	ComlibAssociateProxy(&gSymInstance,pcid->cproxy);
      else
	ComlibAssociateProxy(&gAsymInstance,pcid->cproxy);
      */

    }
  /*  if(pcid->useEtoM)
    if(pcid->Symmetric)
      gSymInstance.beginIteration();
    else
      gAsymInstance.beginIteration();
  */
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
	  if(pcid->useEtoM)
	    {
	      CkArrayIndex4D idx;
	      for(int elem=0; elem<pcid->listRNotFrom.size();elem++)
		{ idx=pcid->listRNotFrom[elem];
		  idx.index[3]=chunk;
		  calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
		  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
		  *(int*)CkPriorityPtr(msg) = pcid->priority;    
		  msg->init(outsize,myS,false,flag_dp,&ptr[chunk*chunksize],false, n);

		  pcid->cproxy(idx.index[0],idx.index[1],idx.index[2],idx.index[3]).acceptPairData(msg);
		}
	    }
	  else
	    {
	      calculatePairsMsg *msg= new ( outsize,8*sizeof(int) ) calculatePairsMsg;
	      CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	      *(int*)CkPriorityPtr(msg) = pcid->priority;    
	      msg->init(outsize,myS,false,flag_dp,&ptr[chunk*chunksize],false, n);

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
  if(pcid->useEtoM)
      gAsymInstance.endIteration();

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

  CkAssert(symmetric == false);
  s2 = (myS/grainSize) * grainSize;
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
						0, nstates-grainSize, grainSize,
						s2, s2, 1,
						c, c, 1);
	  if(pcid->useEtoM && c==0)
	    {

	      for(int s1 = 0; s1 < nstates; s1 += grainSize){
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




void finishPairCalcSection(int n, double *ptr, CProxySection_PairCalculator sectionProxy, int orthoX, int orthoY, int actionType, int priority) {
  finishPairCalcSection2(n, ptr, NULL, sectionProxy, orthoX, orthoY, actionType, priority);
}


/* This version uses a section multicast to only send the part of the matrix needed by each section */
void finishPairCalcSection2(int n, double *ptr1, double *ptr2, CProxySection_PairCalculator sectionProxy, int orthoX, int orthoY, int actionType, int priority) {
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("     Calc Finish Mcast 2\n");
#endif

  if(ptr2==NULL){
    multiplyResultMsg *omsg;
    if(priority>0)
      {
	omsg=new ( n,0,8*sizeof(int) ) multiplyResultMsg;
	*(int*)CkPriorityPtr(omsg) = priority;    
	CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
      }
    else
      {
	omsg=new ( n,0,0 ) multiplyResultMsg;
      }
    omsg->init1(n, ptr1, orthoX, orthoY, actionType);
    sectionProxy.multiplyResult(omsg);
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
	omsg=new ( n,n, 0 ) multiplyResultMsg;
      }
    omsg->init(n, n, ptr1, ptr2, orthoX, orthoY, actionType);
    sectionProxy.multiplyResult(omsg);
  }
}

//
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

