#include "pcCommManager.h"

#include "paircalc/paircalcMessages.h"
#include "ckPairCalculator.decl.h"
#include "paircalc/InputDataHandler.h"
#include "paircalc/pairCalculator.h" //< Just for the reorder_elem declarations
#include "paircalc/pcMaps.h"
#include "utility/MapFile.h"

#include "ckmulticast.h"
#include "ckcomplex.h"

// Delegated paircalc proxies perform like fermented dung on BG/L
#ifdef CMK_BLUEGENEL
#define _PAIRCALC_DO_NOT_DELEGATE_ 1
#endif
// Do not use comlib for multicasts within paircalc
#define _PC_COMMLIB_MULTI_ 0

extern CkVec <MapType4> AsymScalcImaptable;
extern CkVec <MapType4> SymScalcImaptable;
extern CkVec <MapType2> GSImaptable;

namespace cp {
    namespace gspace {

PCCommManager::PCCommManager(const CkIndex2D gspaceIdx, const pc::pcConfig &_cfg):
    gspaceIndex(gspaceIdx), pcCfg(_cfg),
    sectionGettingLeft(0), sectionGettingRight(0),
    existsLproxy(false), existsRproxy(false)
{}




void PCCommManager::createPCarray()
{
    traceRegisterUserEvent("calcpairDGEMM", 210);
    traceRegisterUserEvent("calcpairContrib", 220);
    traceRegisterUserEvent("multiplyResultDGEMM1", 230);
    traceRegisterUserEvent("multiplyResultDGEMM2", 240);
    traceRegisterUserEvent("multiplyResultDGEMM1R", 250);

    CkArrayOptions paircalcOpts,handlerOpts;
    CProxy_PairCalculator pairCalculatorProxy;
    CProxy_InputDataHandler<CollatorType,CollatorType> inputHandlerProxy;

    // Create an empty array but specify element locations using the map
    paircalcOpts.setMap(mapperGID);
    pairCalculatorProxy = CProxy_PairCalculator::ckNew(inputHandlerProxy, pcCfg, paircalcOpts);

    #ifdef DEBUG_CP_PAIRCALC_CREATION
        CkPrintf("createPairCalculator: Creating empty PairCalculator and InputDataHandler chare arrays for %d loop: asymm(0)/symm(1)\n",pcCfg.isSymmetric);
    #endif
    /// Create an empty input handler chare array that will accept all incoming messages from GSpace
    handlerOpts.bindTo(pairCalculatorProxy);
    inputHandlerProxy = CProxy_InputDataHandler<CollatorType,CollatorType> ::ckNew(pairCalculatorProxy,handlerOpts);

    int proc = 0;
    // Initialize my set of array / group IDs
    pcAID = pairCalculatorProxy.ckGetArrayID();
    ipHandlerAID = inputHandlerProxy.ckGetArrayID();
    mCastMgrGID = CProxy_CkMulticastMgr::ckNew(pcCfg.inputSpanningTreeFactor);

    #ifdef USE_COMLIB
        // Setup the appropriate multicast strategy
        Strategy *multistrat = new DirectMulticastStrategy();
        if(pcCfg.isSymmetric)
            mcastInstanceCP=ComlibRegister(multistrat);
        else
            mcastInstanceACP=ComlibRegister(multistrat);
    #endif

    int maxpcstateindex = (pcCfg.numStates / pcCfg.grainSize - 1) * pcCfg.grainSize;
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
                        #ifdef DEBUG_CP_PAIRCALC_CREATION
                            CkPrintf("Inserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,pcCfg.isSymmetric);
                        #endif
                        pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg);
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
                        #ifdef DEBUG_CP_PAIRCALC_CREATION
                            CkPrintf("Inserting PC element [%d %d %d %d %d]\n",numX,s1,s2,c,pcCfg.isSymmetric);
                        #endif
                        pairCalculatorProxy(numX,s1,s2,c).insert(inputHandlerProxy, pcCfg);
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




void PCCommManager::makeLeftTree()
{
	#ifdef DEBUG_CP_PAIRCALC_CREATION
		CkPrintf("GSpace[%d,%d] Making symm(%d) PC array section to receive left data \n", gspaceIndex.x, gspaceIndex.y, pcCfg.isSymmetric);
	#endif

	/// Compute the max index along the state dimensions of the PC array 
	int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1) * pcCfg.grainSize;
	/// Find the row index of the PC chare that handles this state
	int s1 = (gspaceIndex.x/pcCfg.grainSize) * pcCfg.grainSize;
	s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;
	/// If the PC is a symmetric instance, then include only the post-diagonal chares on the row s1, else, include all the PC chares on row s1
	int sColMin = (pcCfg.isSymmetric) ? s1 : 0;
	
	#ifdef DEBUG_CP_PAIRCALC_COMM
		CkPrintf("GSpace[%d,%d] will send left matrix data to symm(%d) PC chares on: Row %d, Cols %d to %d\n", 
															gspaceIndex.x, gspaceIndex.y, pcCfg.isSymmetric,s1,sColMin,maxpcstateindex);
	#endif

	/// If GSpace to PC comm is point to point direct msging
	if(!pcCfg.isInputMulticast)
	{
		/// simply create a list of PC chare array indices, with chunk=0 (as the comm list is the same for all chunks)
		for(int s2 = sColMin; s2 <= maxpcstateindex; s2 += pcCfg.grainSize)
			listGettingLeft.push_back(CkArrayIndex4D(gspaceIndex.y,s1,s2,0));
	}
	/// else, if communication is through section multicasts
	else
	{
		/// Allocate one section proxy for each chunk
		sectionGettingLeft = new CProxySection_InputDataHandler<CollatorType,CollatorType>[pcCfg.numChunks];
		/// Build an array section for each chunk
		for (int chunk = 0; chunk < pcCfg.numChunks; chunk++)
		{
			sectionGettingLeft[chunk] = CProxySection_InputDataHandler<CollatorType,CollatorType>::ckNew(ipHandlerAID,
																gspaceIndex.y, gspaceIndex.y, 1,
																s1, s1, 1,
																sColMin, maxpcstateindex, pcCfg.grainSize,
																chunk, chunk, 1);
			/// Delegate the multicast work to an appropriate library
			#ifndef _PAIRCALC_DO_NOT_DELEGATE_
#ifdef USE_COMLIB
				if(_PC_COMMLIB_MULTI_ )
					ComlibAssociateProxy(mcastInstanceCP,sectionGettingLeft[chunk]);
				else
#endif
				{
					CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastMgrGID).ckLocalBranch();
					sectionGettingLeft[chunk].ckSectionDelegate(mcastGrp);
				}
			#endif
		}
	}
	/// PC chares receiving data have been identified and memorized (either as an array section or a vector of IDs)
	existsLproxy=true;
}




void PCCommManager::makeRightTree()
{
	#ifdef DEBUG_CP_PAIRCALC_CREATION
		CkPrintf("GSpace[%d,%d] Making symm(%d) PC array section to receive right data \n", gspaceIndex.x, gspaceIndex.y, pcCfg.isSymmetric);
	#endif
	
	/// Compute the max index along the state dimensions of the PC array
	int maxpcstateindex=(pcCfg.numStates/ pcCfg.grainSize - 1) * pcCfg.grainSize;
	int s2 = (gspaceIndex.x / pcCfg.grainSize) * pcCfg.grainSize;
	s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;
	/// If the PC is a symmetric instance, then include only the pre-diagonal chares on the column s2 else, include all the PC chares on column s2
	int sRowMax = (pcCfg.isSymmetric) ? s2 - pcCfg.grainSize : maxpcstateindex;
    #ifdef DEBUG_CP_PAIRCALC_COMM
		CkPrintf("GSpace[%d,%d] will send left matrix data to symm(%d) PC chares on: Col %d, Rows %d to %d\n", 
															gspaceIndex.x, gspaceIndex.y, pcCfg.isSymmetric,s2,0,sRowMax);
	#endif
    
	// Accomodate the boundary case: PC chares on the top left [*,0,0,*] of the array shouldnt receive any right data. So dont build any proxies to them
	if (sRowMax >=0)
	{
		/// If GSpace to PC comm is point to point direct msging
		if(!pcCfg.isInputMulticast)
		{
			/// simply create a list of PC chare array indices, with chunk=0 (as the comm list is the same for all chunks)
			for(int s1 = 0; s1 <= sRowMax; s1 += pcCfg.grainSize)
				listGettingRight.push_back(CkArrayIndex4D(gspaceIndex.y,s1,s2,0));
		}
		/// else, if communication is through section multicasts
		else
		{
			/// Allocate one section proxy for each chunk
			sectionGettingRight=new CProxySection_InputDataHandler<CollatorType,CollatorType>[pcCfg.numChunks];
			/// Build an array section for each chunk
			for (int c = 0; c < pcCfg.numChunks; c++)
			{
				sectionGettingRight[c] = CProxySection_InputDataHandler<CollatorType,CollatorType>::ckNew(ipHandlerAID,
																gspaceIndex.y, gspaceIndex.y, 1,
																0, sRowMax, pcCfg.grainSize,
																s2, s2, 1,
																c, c, 1);
				/// Delegate the multicast work to an appropriate library
				#ifndef _PAIRCALC_DO_NOT_DELEGATE_
#ifdef USE_COMLIB
					if(_PC_COMMLIB_MULTI_)
						ComlibAssociateProxy(mcastInstanceCP, sectionGettingRight[c]);
					else
#endif
					{
						CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastMgrGID).ckLocalBranch();
						sectionGettingRight[c].ckSectionDelegate(mcastGrp);
					}
				#endif
			}
		}
		/// PC chares receiving data have been identified and memorized (either as an array section or a vector of IDs)
		existsRproxy=true;
	}
}




/** Deposits data as left matrix block with InputHandler chare array
 * For symmetric instances,  sends to the post-diagonal row of PCs that correspond to state gspaceIndex.x (including the chare on the chare array diagonal)
 * For asymmetric instances, sends to the whole row of PCs that correspond to state gspaceIndex.x
 */
void PCCommManager::sendLeftDataMcast(int n, complex* ptr, bool psiV)
{
    #ifdef PC_USE_RDMA
		/// If RDMA is enabled, we should be here ONLY during PsiV updates
		CkAssert(psiV);
		#ifdef DEBUG_CP_PAIRCALC_RDMA
			CkPrintf("GSpace[%d,%d] Using traditional channels (not RDMA) for psiV left data.\n",gspaceIndex.x,gspaceIndex.y);
		#endif
    #endif

    bool flag_dp = pcCfg.isDoublePackOn;
    /// If a destination array section doesnt exist, build one
    if(!existsLproxy)
    {
        makeLeftTree();
    }
    /// If a left matrix destination section exists, send the data as the left matrix block
    if(existsLproxy)
    {
        int chunksize =  n / pcCfg.numChunks;
        int outsize = chunksize;

        for(int chunk=0; chunk < pcCfg.numChunks ; chunk++)
        {
            // last chunk gets remainder
            if((pcCfg.numChunks > 1) && (chunk == (pcCfg.numChunks - 1)))
                outsize= chunksize + (n % pcCfg.numChunks);
            #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
            if(pcCfg.isSymmetric && gspaceIndex.y==0)
                dumpMatrixDouble("gspPts",(double *)ptr, 1, n*2,gspaceIndex.y,gspaceIndex.x,0,chunk,pcCfg.isSymmetric);
            CkPrintf("L [%d,%d,%d,%d,%d] chunk %d chunksize %d outsize %d for numpoint %d offset will be %d %.12g\n",gspaceIndex.y,gspaceIndex.x, gspaceIndex.x, chunk,pcCfg.isSymmetric, chunk,chunksize,outsize,n,chunk*chunksize,ptr[chunk*chunksize].re);
            #endif
            // If sending directly, use the vector of target PC chares
            if( !pcCfg.isInputMulticast)
            {
                CkArrayIndex4D idx;
                for(int elem=0; elem < listGettingLeft.size() ; elem++)
                {
                    paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                    *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                    idx=listGettingLeft[elem];
                    reinterpret_cast<short*> (idx.data() )[3]=chunk;
                    #ifdef _NAN_CHECK_
                    for(int i=0;i<outsize ;i++)
                    {
                        CkAssert(finite(msg->points[i].re));
                        CkAssert(finite(msg->points[i].im));
                    }
                    #endif
                    CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(ipHandlerAID);
                    handlerProxy(idx).acceptLeftData(msg);
                }
            }
            // else, use a typical multicast to the destination section
            else
            {
                paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, true, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
                if(pcCfg.isSymmetric && gspaceIndex.y==0)
                    dumpMatrixDouble("pairmsg",(double *)msg->points, 1, outsize*2,gspaceIndex.y,gspaceIndex.x,0,chunk,pcCfg.isSymmetric);
                #endif
                #ifdef _NAN_CHECK_
                for(int i=0;i<outsize ;i++)
                {
                    CkAssert(finite(msg->points[i].re));
                    CkAssert(finite(msg->points[i].im));
                }
                #endif
                sectionGettingLeft[chunk].acceptLeftData(msg);
            }
        }
    }
    /// else, if the destination section doesnt exist even after attempting to create one
    else
        CkPrintf("GSpace[%d,%d] No destination symm(%d) PC array section to send left block data [%d,%d,%d,%d,%d] !!!\n",gspaceIndex.x,gspaceIndex.y,pcCfg.isSymmetric);
}




/** Deposits data as right matrix block with InputHandler chare array
 * For symmetric instances,  sends to the strictly pre-diagonal column of PCs that correspond to state gspaceIndex.x
 * For asymmetric instances, sends to the whole row of PCs that correspond to state gspaceIndex.x
 */
void PCCommManager::sendRightDataMcast(int n, complex* ptr, bool psiV)
{
    #ifdef PC_USE_RDMA
		/// If RDMA is enabled, we should be here ONLY during PsiV updates
		CkAssert(psiV);
		#ifdef DEBUG_CP_PAIRCALC_RDMA
			CkPrintf("GSpace[%d,%d] Using traditional channels (not RDMA) for psiV right data.\n",gspaceIndex.x,gspaceIndex.y);
		#endif
    #endif

    bool flag_dp = pcCfg.isDoublePackOn;
    /// If a destination array section doesnt exist, build one
    if(!existsRproxy)
        makeRightTree();
    /// If a right matrix destination section exists, send the data as the left matrix block
    if(existsRproxy)
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
                for(int elem=0; elem<listGettingRight.size();elem++)
                {
                    idx=listGettingRight[elem];
                    reinterpret_cast<short*> (idx.data() )[3]=chunk;
                    paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, false, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                    *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                    #ifdef _NAN_CHECK_
                    for(int i=0;i<outsize ;i++)
                    {
                        CkAssert(finite(msg->points[i].re));
                        CkAssert(finite(msg->points[i].im));
                    }
                    #endif
                    CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(ipHandlerAID);
                    handlerProxy(idx).acceptRightData(msg);
                }
            }
            else
            {
                paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, false, flag_dp, &(ptr[chunk * chunksize]), psiV, n);
                CkSetQueueing(msg, CK_QUEUEING_IFIFO);
                *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
                #ifdef _NAN_CHECK_
                for(int i=0;i<outsize ;i++)
                {
                    CkAssert(finite(msg->points[i].re));
                    CkAssert(finite(msg->points[i].im));
                }
                #endif
                sectionGettingRight[chunk].acceptRightData(msg);
            }
        }
    }
    /// else, if the destination section doesnt exist even after attempting to create one
    else
        CkPrintf("GSpace[%d,%d] No destination symm(%d) PC array section to send right block data [%d,%d,%d,%d,%d] !!!\n",gspaceIndex.x,gspaceIndex.y,pcCfg.isSymmetric);
}




void PCCommManager::sendLeftDataRDMA(int n, complex* ptr, bool psiV)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to send data to paircalcs via RDMA when RDMA is not enabled\n");
	#else
		if(!psiV)
		{
			/// Trigger an RDMA send for every rdma handle associated with all the PCs getting my data as left matrix
			for (int i=0; i< leftDestinationHandles.size();i++)
				if (leftDestinationHandles[i].handle >=0)
				{
					#ifdef DEBUG_CP_PAIRCALC_RDMA
						CkPrintf("GSpace[%d,%d] Sending left data to PC via RDMA.\n",gspaceIndex.x,gspaceIndex.y);
					#endif
					CmiDirect_put( &(leftDestinationHandles[i]) );
				}
		}
		/// else, if it is a PsiV update step, send the data via traditional messaging
		else
		    sendLeftDataMcast(n, ptr, psiV);
	#endif // PC_USE_RDMA
}




void PCCommManager::sendRightDataRDMA(int n, complex* ptr, bool psiV)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to send data to paircalcs via RDMA when RDMA is not enabled\n");
	#else
		if (!psiV)
		{
			/// Trigger an RDMA send for every rdma handle associated with all the PCs getting my data as right matrix
			for (int i=0; i< rightDestinationHandles.size();i++)
				if (rightDestinationHandles[i].handle >=0)
				{
					#ifdef DEBUG_CP_PAIRCALC_RDMA
						CkPrintf("GSpace[%d,%d] Sending right data to PC via RDMA.\n",gspaceIndex.x,gspaceIndex.y);
					#endif
					CmiDirect_put( &(rightDestinationHandles[i]) );
				}
		}
		/// else, if it is a PsiV update step, send the data via traditional messaging
		else
		    sendRightDataMcast(n, ptr, psiV);
	#endif // PC_USE_RDMA
}




void PCCommManager::sendLeftRDMARequest(RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to setup RDMA when RDMA is not enabled\n");
	#else
	#ifdef DEBUG_CP_PAIRCALC_RDMA
		CkPrintf("GSpace[%d,%d] Sending out RDMA setup requests to PCs getting left matrix data from me.\n",idTkn.gspIndex.x,idTkn.gspIndex.y);
	#endif
	/// If the destination PC chares are not known, determine them
	if(!existsLproxy)
		makeLeftTree();
	/// If there exist any destination PC chares
	if(existsLproxy)
	{
		/// Verify
		CkAssert(pcCfg.numChunks > 0);
		/// Compute the size of the chunk of data to be sent out in terms of the number of doubles (as PC treats them) and not complex
		int chunksize  = 2 * (totalsize / pcCfg.numChunks);

		/// Send an RDMA setup request to each destination PC
		for (int chunk=0; chunk < pcCfg.numChunks; chunk++)
		{
			/// The last chunk gets the remainder of the points
			if( (pcCfg.numChunks > 1) && (chunk == pcCfg.numChunks-1) )
				chunksize += 2 * (totalsize % pcCfg.numChunks);
			/// If the communication is through a direct p2p send
			if(!pcCfg.isInputMulticast)
			{
				CkArrayIndex4D idx;
				for(int elem=0; elem < listGettingLeft.size() ; elem++)
				{
					idx=listGettingLeft[elem];
					reinterpret_cast<short*> (idx.data() )[3]=chunk;
					RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
                    CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(ipHandlerAID);
					handlerProxy(idx).setupRDMALeft(msg);
				}
			}
			/// else, if we're multicasting
			else
			{
				RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
				sectionGettingLeft[chunk].setupRDMALeft(msg);
			}
		}
	}
	#endif // PC_USE_RDMA
}




void PCCommManager::sendRightRDMARequest(RDMApair_GSP_PC idTkn, int totalsize, CkCallback cb)
{
	#ifndef PC_USE_RDMA
		CkAbort("GSpace[,] Trying to setup RDMA when RDMA is not enabled\n");
	#else
	#ifdef DEBUG_CP_PAIRCALC_RDMA
		CkPrintf("GSpace[%d,%d] Sending out RDMA setup requests to PCs getting right matrix data from me.\n",idTkn.gspIndex.x,idTkn.gspIndex.y);
	#endif
	/// If the destination PC chares are not known, determine them
	if(!existsRproxy)
		makeRightTree();
	/// If there exist any destination PC chares
	if(existsRproxy)
	{
		/// Verify
		CkAssert(pcCfg.numChunks > 0);
		/// Compute the size of the chunk of data to be sent out in terms of the number of doubles (as PC treats them) and not complex
		int chunksize  = 2 * (totalsize / pcCfg.numChunks);

		/// Send an RDMA setup request to each destination PC
		for (int chunk=0; chunk < pcCfg.numChunks; chunk++)
		{
			/// The last chunk gets the remainder of the points
			if( (pcCfg.numChunks > 1) && (chunk == pcCfg.numChunks-1) )
				chunksize += 2 * (totalsize % pcCfg.numChunks);
			/// If the communication is through a direct p2p send
			if(!pcCfg.isInputMulticast)
			{
				CkArrayIndex4D idx;
				for(int elem=0; elem < listGettingRight.size() ; elem++)
				{
					idx=listGettingRight[elem];
					reinterpret_cast<short*> (idx.data() )[3]=chunk;
					RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
                    CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(ipHandlerAID);
					handlerProxy(idx).setupRDMARight(msg);
				}
			}
			/// else, if we're multicasting
			else
			{
				RDMASetupRequestMsg<RDMApair_GSP_PC> *msg = new RDMASetupRequestMsg<RDMApair_GSP_PC> (idTkn,idTkn.gspIndex.x,CkMyPe(),chunksize,cb);
				sectionGettingRight[chunk].setupRDMARight(msg);
			}
		}
	}
	#endif // PC_USE_RDMA
}




/**
 * send the multcast message to initialize the section tree and set the cookie
 */
void PCCommManager::setResultProxy(CProxySection_PairCalculator *sectProxy, bool lbsync, CkCallback synccb)
{
    int offset = gspaceIndex.x % pcCfg.grainSize;
    int dest = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize; //row or column
    initResultMsg *redMsg=new initResultMsg;
    redMsg->mCastGrpId = mCastMgrGID;
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
CProxySection_PairCalculator PCCommManager::makeOneResultSection_asym(int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastMgrGID).ckLocalBranch();
  int maxpcstateindex = (pcCfg.numStates/pcCfg.grainSize-1) * pcCfg.grainSize;
  int s2 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcAID,
									       gspaceIndex.y, gspaceIndex.y, 1,
									       0, maxpcstateindex, pcCfg.grainSize,
									       s2, s2, 1,
									       chunk, chunk,1);
  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart = gspaceIndex.x % pcCfg.grainSize;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list_max( sid._elems, sid._nElems, newListStart);
  CkAssert(order);
  sectProxy.ckSectionDelegate(mcastGrp);
  //initialize proxy
  setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}

/**
 * initialize  plane and column wise section reduction for lambda->gspace
 */
CProxySection_PairCalculator PCCommManager::makeOneResultSection_asym_column(int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastMgrGID).ckLocalBranch();
  int s1 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize; //column
  int maxpcstateindex = (pcCfg.numStates/pcCfg.grainSize-1) * pcCfg.grainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

  // all nondiagonal elements
  // so we'll have to make this the tedious explicit way

  CkArrayIndex4D *elems= new CkArrayIndex4D[pcCfg.numStates/ pcCfg.grainSize];
  int ecount=0;
  for(int s2 =0; s2<=maxpcstateindex; s2+=pcCfg.grainSize)
    {
      if(s1!=s2)
	{
	  CkArrayIndex4D idx4d(gspaceIndex.y,s1,s2,chunk);
	  elems[ecount++]=idx4d;
	}
    }
  int newListStart = gspaceIndex.x % pcCfg.grainSize;
  if(newListStart> ecount)
    newListStart= newListStart % ecount;
  bool order=reorder_elem_list_4D( elems, ecount, newListStart);
  CkAssert(order);
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcAID,  elems, ecount);
  delete [] elems;
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}



/**
 * initialize  plane and row wise section reduction for psi->gspace
 */
CProxySection_PairCalculator PCCommManager::makeOneResultSection_sym1(int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastMgrGID).ckLocalBranch();
  int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
  int s2 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize;
  s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

  int s2range= (s2==0) ? 1 : pcCfg.grainSize;
  CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcAID,
									       gspaceIndex.y, gspaceIndex.y, 1,
									       0, s2, s2range,
									       s2, s2, 1,
									       chunk, chunk, 1);
  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart = gspaceIndex.x % pcCfg.grainSize;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list_max( sid._elems, sid._nElems, newListStart);
  CkAssert(order);
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}


/**
 * initialize  plane and column wise section reduction for psi->gspace
 */
CProxySection_PairCalculator PCCommManager::makeOneResultSection_sym2(int chunk)
{
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastMgrGID).ckLocalBranch();
  int s1 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize; //column
  int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
  s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

  int s2start = s1 + pcCfg.grainSize;
  s2start= (s2start>maxpcstateindex) ? maxpcstateindex : s2start;
  int s2range= (s2start==maxpcstateindex) ? 1 : pcCfg.grainSize;
  CkAssert(s2start<pcCfg.numStates);
  CProxySection_PairCalculator sectProxy =
      CProxySection_PairCalculator::ckNew(pcAID,
					  gspaceIndex.y, gspaceIndex.y, 1,
					  s1, s1, 1,
					  s2start, maxpcstateindex, s2range,
					  chunk, chunk, 1);

  CkSectionID sid=sectProxy.ckGetSectionID();
  int newListStart = gspaceIndex.x % pcCfg.grainSize;
  if(newListStart> sid._nElems)
    newListStart= newListStart % sid._nElems;
  bool order=reorder_elem_list_max( sid._elems, sid._nElems, newListStart);
  CkAssert(order);
  sectProxy.ckSectionDelegate(mcastGrp);
  setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
  return sectProxy;
}




/**
 * Create the map for placing the paircalculator chare array elements. Also perform other housekeeping chores like dumping the maps to files etc.
 */
void PCCommManager::createMap(const int boxSize, PeListFactory getPeList, UberCollection thisInstance)
{
    int instanceNum = thisInstance.getPO();
    bool maptype = pcCfg.isSymmetric;
    int achunks = config.numChunksAsym;
    if(pcCfg.isSymmetric && pcCfg.arePhantomsOn)
    { // evil trickery to use asym map code for phantom sym
        maptype=false;
        achunks=config.numChunksSym;
    }

    // Generate a map name
    std::string mapName = pcCfg.isSymmetric ? "SymScalcMap" : "AsymScalcMap";
    // Use the appropriate map table
    MapType4 &mapTable  = pcCfg.isSymmetric ? SymScalcImaptable[instanceNum] : AsymScalcImaptable[instanceNum];
    /// Get an appropriately constructed PeList from the supplied factory functor
    PeList *availGlobG = getPeList();
    availGlobG->reset();

    // Compute num PEs along the states dimension of GSpace
    int pl = pcCfg.numStates / config.Gstates_per_pe;
    // Compute num PEs along the planes dimension of GSpace
    int pm = config.numPesPerInstance / pl;
    // Compute the num of GSpace planes per PE
    int planes_per_pe = pcCfg.numPlanes / pm;

    int size[4];
    size[0] = pcCfg.numPlanes;
    size[1] = pcCfg.numStates/pcCfg.grainSize;
    size[2] = pcCfg.numStates/pcCfg.grainSize;
    size[3] = achunks;
    //-------------------------------------------------------------
    // Populate maptable for the PC array

    // Start timing the map creation
    double mapCreationTime = CmiWallTimer();

    mapTable.buildMap(pcCfg.numPlanes, pcCfg.numStates/pcCfg.grainSize, pcCfg.numStates/pcCfg.grainSize, achunks, pcCfg.grainSize);

    int success = 0;
    if(config.loadMapFiles)
    {
        MapFile *mf = new MapFile(mapName.c_str(), 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, pcCfg.grainSize);
        success = mf->loadMap(mapName.c_str(), &mapTable);
        delete mf;
    }

    // If loading the map from a file failed, create a maptable
    if(success == 0)
    {
        SCalcMapTable symTable = SCalcMapTable(&mapTable, availGlobG, pcCfg.numStates, pcCfg.numPlanes, pcCfg.grainSize, maptype, config.scalc_per_plane,
                        planes_per_pe, achunks, config.numChunksSym, &GSImaptable[instanceNum], config.useCuboidMap, config.useCentroidMap, boxSize);
    }

    /// Create a map group that will read and use this map table
    CProxy_SCalcMap pcMapGrp = CProxy_SCalcMap::ckNew(pcCfg.isSymmetric, thisInstance);

    mapCreationTime = CmiWallTimer() - mapCreationTime;
    CkPrintf("PairCalculator[%dx%dx%dx%d,%d] map created in %g\n", size[0], size[1], size[2], pcCfg.numChunks, pcCfg.isSymmetric, mapCreationTime);

    // If the user wants map dumps
    if(config.dumpMapFiles)
    {
        MapFile *mf = new MapFile(mapName.c_str(), 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, pcCfg.grainSize);
        mf->dumpMap(&mapTable, instanceNum);
        delete mf;
    }

    // If the user wants map coordinate dumps
    if(config.dumpMapCoordFiles)
    {
        MapFile *mf = new MapFile((mapName+"_coord").c_str(), 4, size, config.numPes, "TXYZ", 2, 1, 1, 1, pcCfg.grainSize);
        mf->dumpMapCoords(&mapTable, instanceNum);
        delete mf;
    }

    // Record the group that will provide the procNum mapping function
    mapperGID  = pcMapGrp.ckGetGroupID();
    delete availGlobG;
}

    } // end namespace gspace
} // end namespace cp

#include "RDMAMessages.def.h"
#include "pcMaps.def.h"

