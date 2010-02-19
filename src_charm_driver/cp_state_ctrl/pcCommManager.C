#include "pcCommManager.h"
#include "paircalc/pairCalculator.h"
#include "paircalc/InputDataHandler.h"
#include "ckcomplex.h"

namespace cp {
    namespace gspace {

void PCCommManager::createPCarray(const pc::pcConfig pcCfg, PairCalcID* pcid, CkGroupID *mapid)
{
    traceRegisterUserEvent("calcpairDGEMM", 210);
    traceRegisterUserEvent("calcpairContrib", 220);
    traceRegisterUserEvent("multiplyResultDGEMM1", 230);
    traceRegisterUserEvent("multiplyResultDGEMM2", 240);
    traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

    CkArrayOptions paircalcOpts,handlerOpts;
    CProxy_PairCalculator pairCalculatorProxy;
    CProxy_InputDataHandler<CollatorType,CollatorType> inputHandlerProxy;

    // If a chare mapping is not available, just create an empty array
    if(!mapid)
        pairCalculatorProxy = CProxy_PairCalculator::ckNew();
    // else, create an empty array but specify element locations using the map
    else
    {
        paircalcOpts.setMap(*mapid);
        pairCalculatorProxy = CProxy_PairCalculator::ckNew(inputHandlerProxy, pcCfg, paircalcOpts);
    }

    #ifdef DEBUG_CP_PAIRCALC_CREATION
        CkPrintf("createPairCalculator: Creating empty PairCalculator and InputDataHandler chare arrays for %d loop: asymm(0)/symm(1)\n",pcCfg.isSymmetric);
    #endif
    /// Create an empty input handler chare array that will accept all incoming messages from GSpace
    handlerOpts.bindTo(pairCalculatorProxy);
    inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);

    int proc = 0;
    // Initialize the PairCalcID instance
    pcid->Init(pairCalculatorProxy.ckGetArrayID(), inputHandlerProxy.ckGetArrayID());
    pcid->mCastGrpId = CProxy_CkMulticastMgr::ckNew(pcCfg.inputSpanningTreeFactor);

    #ifdef USE_COMLIB
        // Setup the appropriate multicast strategy
        Strategy *multistrat = new DirectMulticastStrategy();
        if(pcCfg.isSymmetric)
            mcastInstanceCP=ComlibRegister(multistrat);
        else
            mcastInstanceACP=ComlibRegister(multistrat);
    #endif

    int maxpcstateindex = (pcCfg.numStates / pcCfg.grainSize-1) * pcCfg.grainSize;
    CkAssert(mapid);
    // If the symmetric loop PC instances are being created
    if(pcCfg.isSymmetric)
        for(int numX = 0; numX < pcCfg.numPlanes; numX ++)
        {
            for (int s1 = 0; s1 <= maxpcstateindex; s1 += pcCfg.grainSize)
            {
                // If phantomSym is turned on
                int s2start=(pcCfg.arePhantomsOn) ? 0 : s1;
                for (int s2 = s2start; s2 <= maxpcstateindex; s2 += pcCfg.grainSize)
                {
                    for (int c = 0; c < pcCfg.numChunks; c++)
                    {
                        if(mapid)
                        {
                            #ifdef DEBUG_CP_PAIRCALC_CREATION
                                CkPrintf("Inserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,pcCfg.isSymmetric);
                            #endif
                            pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg);
                        }
                        else
                        {
                            #ifdef DEBUG_CP_PAIRCALC_CREATION
                                CkPrintf("Inserting PC element [%d %d %d %d %d] at PE %d\n",numX,s1,s2,c,pcCfg.isSymmetric,proc);
                            #endif
                            pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg, proc);
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
        for(int numX = 0; numX < pcCfg.numPlanes; numX ++)
        {
            for (int s1 = 0; s1 <= maxpcstateindex; s1 += pcCfg.grainSize)
            {
                for (int s2 = 0; s2 <= maxpcstateindex; s2 += pcCfg.grainSize)
                {
                    for (int c = 0; c < pcCfg.numChunks; c++)
                    {
                        if(mapid)
                        {
                            #ifdef DEBUG_CP_PAIRCALC_CREATION
                                CkPrintf("Inserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,pcCfg.isSymmetric);
                            #endif
                            pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg);
                        }
                        else
                        {
                            #ifdef DEBUG_CP_PAIRCALC_CREATION
                                CkPrintf("Inserting PC element [%d %d %d %d %d] on PE %d\n",numX,s1,s2,c,pcCfg.isSymmetric,proc);
                            #endif
                            pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg, proc);
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
        CkPrintf("    Finished init {grain=%d, sym=%d, blk=%d, Z=%d, S=%d}\n", pcCfg.grainSize, pcCfg.isSymmetric, pcCfg.numChunks, pcCfg.numPlanes, pcCfg.numStates);
    #endif
}




void PCCommManager::makeLeftTree(PairCalcID* pcid, int myS, int myPlane)
{
	#ifdef DEBUG_CP_PAIRCALC_CREATION
		CkPrintf("GSpace[%d,%d] Making symm(%d) PC array section to receive left data \n", myS, myPlane, pcCfg.isSymmetric);
	#endif

	int grainSize = pcCfg.grainSize;
	int numChunks = pcCfg.numChunks;
	/// Compute the max index along the state dimensions of the PC array 
	int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1) * pcCfg.grainSize;
	/// Find the row index of the PC chare that handles this state
	int s1 = (myS/pcCfg.grainSize) * pcCfg.grainSize;
	s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
	/// If the PC is a symmetric instance, then include only the post-diagonal chares on the row s1, else, include all the PC chares on row s1
	int sColMin = (pcCfg.isSymmetric) ? s1 : 0;
	
	#ifdef DEBUG_CP_PAIRCALC_COMM
		CkPrintf("GSpace[%d,%d] will send left matrix data to symm(%d) PC chares on: Row %d, Cols %d to %d\n", 
															myS, myPlane, pcCfg.isSymmetric,s1,sColMin,maxpcstateindex);
	#endif

	/// If GSpace to PC comm is point to point direct msging
	if(!pcCfg.isInputMulticast)
	{
		/// simply create a list of PC chare array indices, with chunk=0 (as the comm list is the same for all chunks)
		for(int s2 = sColMin; s2 <= maxpcstateindex; s2 += grainSize)
			pcid->listGettingLeft.push_back(CkArrayIndex4D(myPlane,s1,s2,0));
	}
	/// else, if communication is through section multicasts
	else
	{
		/// Allocate one section proxy for each chunk
		pcid->sectionGettingLeft=new CProxySection_InputDataHandler<CollatorType,CollatorType>[numChunks];
		/// Build an array section for each chunk
		for (int chunk = 0; chunk < numChunks; chunk++)
		{
			pcid->sectionGettingLeft[chunk] = CProxySection_InputDataHandler<CollatorType,CollatorType>::ckNew(pcid->ipHandlerID,
																myPlane, myPlane, 1,
																s1, s1, 1,
																sColMin, maxpcstateindex, grainSize,
																chunk, chunk, 1);
			/// Delegate the multicast work to an appropriate library
			#ifndef _PAIRCALC_DO_NOT_DELEGATE_
#ifdef USE_COMLIB
				if(pcid->useComlib && _PC_COMMLIB_MULTI_ )
					ComlibAssociateProxy(mcastInstanceCP,pcid->sectionGettingLeft[chunk]);
				else
#endif
				{
					CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();
					pcid->sectionGettingLeft[chunk].ckSectionDelegate(mcastGrp);
				}
			#endif
		}
	}
	/// PC chares receiving data have been identified and memorized (either as an array section or a vector of IDs)
	pcid->existsLproxy=true;
	/** @todo: Removing this bit of code crashes the application. It seems necessary to get fresh proxies at 
	 * this point, instead of simply creating them in the pcid object's Init() method. Perhaps, getting a proxy
	 * to an empty chare array is not a robust operation. Or there could be some other issue. Look into this. 
	 */
	pcid->handlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> (pcid->ipHandlerID);
}




void PCCommManager::makeRightTree(PairCalcID* pcid, int myS, int myPlane)
{
	#ifdef DEBUG_CP_PAIRCALC_CREATION
		CkPrintf("GSpace[%d,%d] Making symm(%d) PC array section to receive right data \n", myS, myPlane, pcCfg.isSymmetric);
	#endif
	
	int grainSize = pcCfg.grainSize;
	int numChunks =  pcCfg.numChunks;
	/// Compute the max index along the state dimensions of the PC array
	int maxpcstateindex=(pcCfg.numStates/grainSize-1) * grainSize;
	int s2 = (myS/grainSize) * grainSize;
	s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
	/// If the PC is a symmetric instance, then include only the pre-diagonal chares on the column s2 else, include all the PC chares on column s2
	int sRowMax = (pcCfg.isSymmetric) ? s2-grainSize : maxpcstateindex;
    #ifdef DEBUG_CP_PAIRCALC_COMM
		CkPrintf("GSpace[%d,%d] will send left matrix data to symm(%d) PC chares on: Col %d, Rows %d to %d\n", 
															myS, myPlane, pcCfg.isSymmetric,s2,0,sRowMax);
	#endif
    
	// Accomodate the boundary case: PC chares on the top left [*,0,0,*] of the array shouldnt receive any right data. So dont build any proxies to them
	if (sRowMax >=0)
	{
		/// If GSpace to PC comm is point to point direct msging
		if(!pcCfg.isInputMulticast)
		{
			/// simply create a list of PC chare array indices, with chunk=0 (as the comm list is the same for all chunks)
			for(int s1 = 0; s1 <= sRowMax; s1 += grainSize)
				pcid->listGettingRight.push_back(CkArrayIndex4D(myPlane,s1,s2,0));
		}
		/// else, if communication is through section multicasts
		else
		{
			/// Allocate one section proxy for each chunk
			pcid->sectionGettingRight=new CProxySection_InputDataHandler<CollatorType,CollatorType>[numChunks];
			/// Build an array section for each chunk
			for (int c = 0; c < numChunks; c++)
			{
				pcid->sectionGettingRight[c] = CProxySection_InputDataHandler<CollatorType,CollatorType>::ckNew(pcid->ipHandlerID,
																myPlane, myPlane, 1,
																0, sRowMax, grainSize,
																s2, s2, 1,
																c, c, 1);
				/// Delegate the multicast work to an appropriate library
				#ifndef _PAIRCALC_DO_NOT_DELEGATE_
#ifdef USE_COMLIB
					if(pcid->useComlib && _PC_COMMLIB_MULTI_)
						ComlibAssociateProxy(mcastInstanceCP, pcid->sectionGettingRight[c]);
					else
#endif
					{
						CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();
						pcid->sectionGettingRight[c].ckSectionDelegate(mcastGrp);
					}
				#endif
			}
		}
		/// PC chares receiving data have been identified and memorized (either as an array section or a vector of IDs)
		pcid->existsRproxy=true;
	}
}




/** Deposits data as left matrix block with InputHandler chare array
 * For symmetric instances,  sends to the post-diagonal row of PCs that correspond to state myS (including the chare on the chare array diagonal)
 * For asymmetric instances, sends to the whole row of PCs that correspond to state myS
 */
void PCCommManager::sendLeftDataMcast(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
    #ifdef PC_USE_RDMA
		/// If RDMA is enabled, we should be here ONLY during PsiV updates
		CkAssert(psiV);
		#ifdef DEBUG_CP_PAIRCALC_RDMA
			CkPrintf("GSpace[%d,%d] Using traditional channels (not RDMA) for psiV left data.\n",myS,myPlane);
		#endif
    #endif

    bool flag_dp = pcCfg.isDoublePackOn;
    /// If a destination array section doesnt exist, build one
    if(!pcid->existsLproxy)
    {
        makeLeftTree(pcid,myS,myPlane);
    }
    /// If a left matrix destination section exists, send the data as the left matrix block
    if(pcid->existsLproxy)
    {
        int numChunks=pcCfg.numChunks;
        int chunksize =  n / numChunks;
        int outsize = chunksize;

        for(int chunk=0; chunk < numChunks ; chunk++)
        {
            // last chunk gets remainder
            if((numChunks > 1) && (chunk == (numChunks - 1)))
                outsize= chunksize + (n % numChunks);
            #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
            if(pcCfg.isSymmetric && myPlane==0)
                dumpMatrixDouble("gspPts",(double *)ptr, 1, n*2,myPlane,myS,0,chunk,pcCfg.isSymmetric);
            CkPrintf("L [%d,%d,%d,%d,%d] chunk %d chunksize %d outsize %d for numpoint %d offset will be %d %.12g\n",myPlane,myS, myS, chunk,pcCfg.isSymmetric, chunk,chunksize,outsize,n,chunk*chunksize,ptr[chunk*chunksize].re);
            #endif
            // If sending directly, use the vector of target PC chares
            if( !pcCfg.isInputMulticast)
            {
                CkArrayIndex4D idx;
                for(int elem=0; elem < pcid->listGettingLeft.size() ; elem++)
                {
                    paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, myS, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                    *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                    idx=pcid->listGettingLeft[elem];
                    reinterpret_cast<short*> (idx.data() )[3]=chunk;
                    #ifdef _NAN_CHECK_
                    for(int i=0;i<outsize ;i++)
                    {
                        CkAssert(finite(msg->points[i].re));
                        CkAssert(finite(msg->points[i].im));
                    }
                    #endif
                    pcid->handlerProxy(idx).acceptLeftData(msg);
                }
            }
            // else, use a typical multicast to the destination section
            else
            {
                paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, myS, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
                if(pcCfg.isSymmetric && myPlane==0)
                    dumpMatrixDouble("pairmsg",(double *)msg->points, 1, outsize*2,myPlane,myS,0,chunk,pcCfg.isSymmetric);
                #endif
                #ifdef _NAN_CHECK_
                for(int i=0;i<outsize ;i++)
                {
                    CkAssert(finite(msg->points[i].re));
                    CkAssert(finite(msg->points[i].im));
                }
                #endif
                pcid->sectionGettingLeft[chunk].acceptLeftData(msg);
            }
        }
    }
		/// else, if the destination section doesnt exist even after attempting to create one
		else
			CkPrintf("GSpace[%d,%d] No destination symm(%d) PC array section to send left block data [%d,%d,%d,%d,%d] !!!\n",myS,myPlane,pcCfg.isSymmetric);
}




/** Deposits data as right matrix block with InputHandler chare array
 * For symmetric instances,  sends to the strictly pre-diagonal column of PCs that correspond to state myS
 * For asymmetric instances, sends to the whole row of PCs that correspond to state myS
 */
void PCCommManager::sendRightDataMcast(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
    #ifdef PC_USE_RDMA
		/// If RDMA is enabled, we should be here ONLY during PsiV updates
		CkAssert(psiV);
		#ifdef DEBUG_CP_PAIRCALC_RDMA
			CkPrintf("GSpace[%d,%d] Using traditional channels (not RDMA) for psiV right data.\n",myS,myPlane);
		#endif
    #endif

    bool flag_dp = pcCfg.isDoublePackOn;
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
        for(int chunk=0; chunk < pcCfg.numChunks; chunk++)
        {
            int chunksize=n/pcCfg.numChunks;
            int outsize=chunksize;
            /// last chunk gets remainder
            if(pcCfg.numChunks > 1 && chunk == pcCfg.numChunks - 1)
                outsize+=n % pcCfg.numChunks;
            if(!pcCfg.isInputMulticast)
            {
                CkArrayIndex4D idx;
                for(int elem=0; elem<pcid->listGettingRight.size();elem++)
                {
                    idx=pcid->listGettingRight[elem];
                    reinterpret_cast<short*> (idx.data() )[3]=chunk;
                    paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, myS, false, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                    *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                    #ifdef _NAN_CHECK_
                    for(int i=0;i<outsize ;i++)
                    {
                        CkAssert(finite(msg->points[i].re));
                        CkAssert(finite(msg->points[i].im));
                    }
                    #endif
                    pcid->handlerProxy(idx).acceptRightData(msg);
                }
            }
            else
            {
                paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, myS, false, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                #ifdef _NAN_CHECK_
                for(int i=0;i<outsize ;i++)
                {
                    CkAssert(finite(msg->points[i].re));
                    CkAssert(finite(msg->points[i].im));
                }
                #endif
                pcid->sectionGettingRight[chunk].acceptRightData(msg);
            }
        }
    }
    /// else, if the destination section doesnt exist even after attempting to create one
    else
        CkPrintf("GSpace[%d,%d] No destination symm(%d) PC array section to send right block data [%d,%d,%d,%d,%d] !!!\n",myS,myPlane,pcCfg.isSymmetric);
}




void PCCommManager::sendLeftDataRDMA(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to send data to paircalcs via RDMA when RDMA is not enabled\n");
	#else
		if(!psiV)
		{
			/// Trigger an RDMA send for every rdma handle associated with all the PCs getting my data as left matrix
			for (int i=0; i< pcid->leftDestinationHandles.size();i++)
				if (pcid->leftDestinationHandles[i].handle >=0)
				{
					#ifdef DEBUG_CP_PAIRCALC_RDMA
						CkPrintf("GSpace[%d,%d] Sending left data to PC via RDMA.\n",myS,myPlane);
					#endif
					CmiDirect_put( &(pcid->leftDestinationHandles[i]) );
				}
		}
		/// else, if it is a PsiV update step, send the data via traditional messaging
		else
		    sendLeftData(pcid, n, ptr, myS, myPlane, psiV);
	#endif // PC_USE_RDMA
}




void PCCommManager::sendRightDataRDMA(PairCalcID* pcid, int n, complex* ptr, int myS, int myPlane, bool psiV)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to send data to paircalcs via RDMA when RDMA is not enabled\n");
	#else
		if (!psiV)
		{
			/// Trigger an RDMA send for every rdma handle associated with all the PCs getting my data as right matrix
			for (int i=0; i< pcid->rightDestinationHandles.size();i++)
				if (pcid->rightDestinationHandles[i].handle >=0)
				{
					#ifdef DEBUG_CP_PAIRCALC_RDMA
						CkPrintf("GSpace[%d,%d] Sending right data to PC via RDMA.\n",myS,myPlane);
					#endif
					CmiDirect_put( &(pcid->rightDestinationHandles[i]) );
				}
		}
		/// else, if it is a PsiV update step, send the data via traditional messaging
		else
		    sendRightData(pcid, n, ptr, myS, myPlane, psiV);
	#endif // PC_USE_RDMA
}




void PCCommManager::sendLeftRDMARequest(PairCalcID *pid, RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to setup RDMA when RDMA is not enabled\n");
	#else
	#ifdef DEBUG_CP_PAIRCALC_RDMA
		CkPrintf("GSpace[%d,%d] Sending out RDMA setup requests to PCs getting left matrix data from me.\n",idTkn.gspIndex.x,idTkn.gspIndex.y);
	#endif
	/// If the destination PC chares are not known, determine them
	if(!pid->existsLproxy)
		makeLeftTree(pid,idTkn.gspIndex.x,idTkn.gspIndex.y);
	/// If there exist any destination PC chares
	if(pid->existsLproxy)
	{
		/// Verify
		CkAssert(pid->numChunks > 0);
		/// Compute the size of the chunk of data to be sent out in terms of the number of doubles (as PC treats them) and not complex
		int chunksize  = 2 * (totalsize / pid->numChunks);

		/// Send an RDMA setup request to each destination PC
		for (int chunk=0; chunk < pid->numChunks; chunk++)
		{
			/// The last chunk gets the remainder of the points
			if( (pid->numChunks > 1) && (chunk == pid->numChunks-1) )
				chunksize += 2 * (totalsize % pid->numChunks);
			/// If the communication is through a direct p2p send
			if(pid->useDirectSend)
			{
				CkArrayIndex4D idx;
				for(int elem=0; elem < pid->listGettingLeft.size() ; elem++)
				{
					idx=pid->listGettingLeft[elem];
					reinterpret_cast<short*> (idx.data() )[3]=chunk;
					RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
					pid->handlerProxy(idx).setupRDMALeft(msg);
				}
			}
			/// else, if we're multicasting
			else
			{
				RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
				pid->sectionGettingLeft[chunk].setupRDMALeft(msg);
			}
		}
	}
	#endif // PC_USE_RDMA
}




void PCCommManager::sendRightRDMARequest(PairCalcID *pid, RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to setup RDMA when RDMA is not enabled\n");
	#else
	#ifdef DEBUG_CP_PAIRCALC_RDMA
		CkPrintf("GSpace[%d,%d] Sending out RDMA setup requests to PCs getting right matrix data from me.\n",idTkn.gspIndex.x,idTkn.gspIndex.y);
	#endif
	/// If the destination PC chares are not known, determine them
	if(!pid->existsRproxy)
		makeRightTree(pid,idTkn.gspIndex.x,idTkn.gspIndex.y);
	/// If there exist any destination PC chares
	if(pid->existsRproxy)
	{
		/// Verify
		CkAssert(pid->numChunks > 0);
		/// Compute the size of the chunk of data to be sent out in terms of the number of doubles (as PC treats them) and not complex
		int chunksize  = 2 * (totalsize / pid->numChunks);

		/// Send an RDMA setup request to each destination PC
		for (int chunk=0; chunk < pid->numChunks; chunk++)
		{
			/// The last chunk gets the remainder of the points
			if( (pid->numChunks > 1) && (chunk == pid->numChunks-1) )
				chunksize += 2 * (totalsize % pid->numChunks);
			/// If the communication is through a direct p2p send
			if(pid->useDirectSend)
			{
				CkArrayIndex4D idx;
				for(int elem=0; elem < pid->listGettingRight.size() ; elem++)
				{
					idx=pid->listGettingRight[elem];
					reinterpret_cast<short*> (idx.data() )[3]=chunk;
					RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
					pid->handlerProxy(idx).setupRDMARight(msg);
				}
			}
			/// else, if we're multicasting
			else
			{
				RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
				pid->sectionGettingRight[chunk].setupRDMARight(msg);
			}
		}
	}
	#endif // PC_USE_RDMA
}




/**
 * send the multcast message to initialize the section tree and set the cookie
 */
void PCCommManager::setResultProxy(CProxySection_PairCalculator *sectProxy, int state, int GrainSize, CkGroupID mCastGrpId, bool lbsync, CkCallback synccb)
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
CProxySection_PairCalculator PCCommManager::makeOneResultSection_asym(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();
  int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
  int s2=state/pcCfg.grainSize*pcCfg.grainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
  int nstates=pcCfg.numStates;
  int GrainSize=pcCfg.grainSize;
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
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * initialize  plane and column wise section reduction for lambda->gspace
 */
CProxySection_PairCalculator PCCommManager::makeOneResultSection_asym_column(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();
  int GrainSize=pcCfg.grainSize;
  int s1=state / GrainSize * GrainSize; //column
  int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

  int nstates=pcCfg.numStates;
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
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}



/**
 * initialize  plane and row wise section reduction for psi->gspace
 */
CProxySection_PairCalculator PCCommManager::makeOneResultSection_sym1(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();
  int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
  int s2=state/pcCfg.grainSize*pcCfg.grainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

  int GrainSize=pcCfg.grainSize;
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
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}


/**
 * initialize  plane and column wise section reduction for psi->gspace
 */
CProxySection_PairCalculator PCCommManager::makeOneResultSection_sym2(PairCalcID* pcid, int state, int plane, int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcid->mCastGrpId).ckLocalBranch();
  int GrainSize=pcCfg.grainSize;
  int s1=state / GrainSize * GrainSize; //column
  int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

  int nstates=pcCfg.numStates;
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
  setResultProxy(&sectProxy, state, GrainSize, pcid->mCastGrpId, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

    } // end namespace gspace
} // end namespace cp

