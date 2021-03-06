#include "pcCommManager.h"

#include "paircalc/pcMessages.h"
#include "ckPairCalculator.decl.h"
#include "paircalc/InputDataHandler.h"
#include "utility/matrix2file.h"

#include "ckmulticast.h"
#include "ckcomplex.h"

#include <algorithm>

// Do not use comlib for multicasts within paircalc
#define _PC_COMMLIB_MULTI_ 0
/** @addtogroup GSpaceState
  @{
 */
namespace cp {
  namespace gspace {

    PCCommManager::PCCommManager(const CkIndex2D gspaceIdx, const pc::pcConfig &_cfg, const pc::InstanceIDs _pcHandle):
      gspaceIndex(gspaceIdx), pcCfg(_cfg), pcHandle(_pcHandle),
      sectionGettingLeft(0), sectionGettingRight(0),
      existsLproxy(false), existsRproxy(false)
    {}




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
          sectionGettingLeft[chunk] = CProxySection_InputDataHandler<CollatorType,CollatorType>::ckNew(pcHandle.handlerAID,
              gspaceIndex.y, gspaceIndex.y, 1,
              s1, s1, 1,
              sColMin, maxpcstateindex, pcCfg.grainSize,
              chunk, chunk, 1);
          /// Delegate the multicast work to an appropriate library
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
#ifndef _AUTO_DELEGATE_MCASTMGR_ON_
          CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcHandle.mCastMgrGID).ckLocalBranch();
          sectionGettingLeft[chunk].ckSectionDelegate(mcastGrp);
#endif
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
            sectionGettingRight[c] = CProxySection_InputDataHandler<CollatorType,CollatorType>::ckNew(pcHandle.handlerAID,
                gspaceIndex.y, gspaceIndex.y, 1,
                0, sRowMax, pcCfg.grainSize,
                s2, s2, 1,
                c, c, 1);
            /// Delegate the multicast work to an appropriate library
#ifndef _PAIRCALC_DO_NOT_DELEGATE_
#ifndef _AUTO_DELEGATE_MCASTMGR_ON_
            CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcHandle.mCastMgrGID).ckLocalBranch();
            sectionGettingRight[c].ckSectionDelegate(mcastGrp);
#endif
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
    void PCCommManager::sendLeftDataMcast(int numPoints, complex* ptr, bool psiV)
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
        int chunksize =  numPoints / pcCfg.numChunks;
        int outsize = chunksize;

        for(int chunk=0; chunk < pcCfg.numChunks ; chunk++)
        {
          // last chunk gets remainder
          if((pcCfg.numChunks > 1) && (chunk == (pcCfg.numChunks - 1)))
            outsize= chunksize + (numPoints % pcCfg.numChunks);
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
          if(pcCfg.isSymmetric && gspaceIndex.y==0)
            dumpMatrix("gspPts",(double *)ptr, 1, numPoints*2,gspaceIndex.y,gspaceIndex.x,0,chunk,pcCfg.isSymmetric);
          CkPrintf("L [%d,%d,%d,%d,%d] chunk %d chunksize %d outsize %d for numpoint %d offset will be %d %.12g\n",gspaceIndex.y,gspaceIndex.x, gspaceIndex.x, chunk,pcCfg.isSymmetric, chunk,chunksize, outsize, numPoints, chunk*chunksize, ptr[chunk*chunksize].re);
#endif
          // If sending directly, use the vector of target PC chares
          if( !pcCfg.isInputMulticast)
          {
            CkArrayIndex4D idx;
            for(int elem=0; elem < listGettingLeft.size() ; elem++)
            {
              paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, true, flag_dp, &(ptr[chunk * chunksize]), psiV, numPoints);
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
              CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(pcHandle.handlerAID);
              handlerProxy(idx).acceptLeftData(msg);
            }
          }
          // else, use a typical multicast to the destination section
          else
          {
            paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, true, flag_dp, &(ptr[chunk * chunksize]), psiV, numPoints);
            *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
            CkSetQueueing(msg, CK_QUEUEING_IFIFO);
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
            if(pcCfg.isSymmetric && gspaceIndex.y==0)
              dumpMatrix("pairmsg",(double *)msg->points, 1, outsize*2,gspaceIndex.y,gspaceIndex.x,0,chunk,pcCfg.isSymmetric);
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
    void PCCommManager::sendRightDataMcast(int numPoints, complex* ptr, bool psiV)
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
        for(int i=0;i<numPoints;i++)
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
          int chunksize = numPoints / pcCfg.numChunks;
          int outsize   = chunksize;
          /// last chunk gets remainder
          if(pcCfg.numChunks > 1 && chunk == pcCfg.numChunks - 1)
            outsize += numPoints % pcCfg.numChunks;
          if(!pcCfg.isInputMulticast)
          {
            CkArrayIndex4D idx;
            for(int elem=0; elem<listGettingRight.size();elem++)
            {
              idx=listGettingRight[elem];
              reinterpret_cast<short*> (idx.data() )[3]=chunk;
              paircalcInputMsg *msg=new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, false, flag_dp, &(ptr[chunk * chunksize]), psiV, numPoints);
              CkSetQueueing(msg, CK_QUEUEING_IFIFO);
              *(int*)CkPriorityPtr(msg) = pcCfg.inputMsgPriority;
#ifdef _NAN_CHECK_
              for(int i=0;i<outsize ;i++)
              {
                CkAssert(finite(msg->points[i].re));
                CkAssert(finite(msg->points[i].im));
              }
#endif
              CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(pcHandle.handlerAID);
              handlerProxy(idx).acceptRightData(msg);
            }
          }
          else
          {
            paircalcInputMsg *msg = new (outsize, 8* sizeof(int)) paircalcInputMsg(outsize, gspaceIndex.x, false, flag_dp, &(ptr[chunk * chunksize]), psiV, numPoints);
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




    void PCCommManager::sendLeftDataRDMA(int numPoints, complex* ptr, bool psiV)
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
        sendLeftDataMcast(numPoints, ptr, psiV);
#endif // PC_USE_RDMA
    }




    void PCCommManager::sendRightDataRDMA(int numPoints, complex* ptr, bool psiV)
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
        sendRightDataMcast(numPoints, ptr, psiV);
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
              CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(pcHandle.handlerAID);
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
              CProxy_InputDataHandler<CollatorType,CollatorType> handlerProxy(pcHandle.handlerAID);
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
      redMsg->mCastGrpId = pcHandle.mCastMgrGID;
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
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcHandle.mCastMgrGID).ckLocalBranch();
      int maxpcstateindex = (pcCfg.numStates/pcCfg.grainSize-1) * pcCfg.grainSize;
      int s2 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize;
      s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

      CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcHandle.pcAID,
          gspaceIndex.y, gspaceIndex.y, 1,
          0, maxpcstateindex, pcCfg.grainSize,
          s2, s2, 1,
          chunk, chunk,1);
      CkSectionID sid=sectProxy.ckGetSectionID();
      std::random_shuffle(sid._elems, sid._elems + sid._nElems);
#ifndef _AUTO_DELEGATE_MCASTMGR_ON_
      sectProxy.ckSectionDelegate(mcastGrp);
#endif
      //initialize proxy
      setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
      return sectProxy;
    }

    /**
     * initialize  plane and column wise section reduction for lambda->gspace
     */
    CProxySection_PairCalculator PCCommManager::makeOneResultSection_asym_column(int chunk)
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcHandle.mCastMgrGID).ckLocalBranch();
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
      std::random_shuffle(elems, elems + ecount);
      CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcHandle.pcAID,  elems, ecount);
      delete [] elems;
#ifndef _AUTO_DELEGATE_MCASTMGR_ON_
      sectProxy.ckSectionDelegate(mcastGrp);
#endif
      setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
      return sectProxy;
    }



    /**
     * initialize  plane and row wise section reduction for psi->gspace
     */
    CProxySection_PairCalculator PCCommManager::makeOneResultSection_sym1(int chunk)
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcHandle.mCastMgrGID).ckLocalBranch();
      int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
      int s2 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize;
      s2 = (s2>maxpcstateindex) ? maxpcstateindex :s2;

      int s2range= (s2==0) ? 1 : pcCfg.grainSize;
      CProxySection_PairCalculator sectProxy = CProxySection_PairCalculator::ckNew(pcHandle.pcAID,
          gspaceIndex.y, gspaceIndex.y, 1,
          0, s2, s2range,
          s2, s2, 1,
          chunk, chunk, 1);
      CkSectionID sid=sectProxy.ckGetSectionID();
      std::random_shuffle(sid._elems, sid._elems + sid._nElems);
#ifndef _AUTO_DELEGATE_MCASTMGR_ON_
      sectProxy.ckSectionDelegate(mcastGrp);
#endif
      setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
      return sectProxy;
    }


    /**
     * initialize  plane and column wise section reduction for psi->gspace
     */
    CProxySection_PairCalculator PCCommManager::makeOneResultSection_sym2(int chunk)
    {
      CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(pcHandle.mCastMgrGID).ckLocalBranch();
      int s1 = gspaceIndex.x / pcCfg.grainSize * pcCfg.grainSize; //column
      int maxpcstateindex=(pcCfg.numStates/pcCfg.grainSize-1)*pcCfg.grainSize;
      s1 = (s1>maxpcstateindex) ? maxpcstateindex :s1;

      int s2start = s1 + pcCfg.grainSize;
      s2start= (s2start>maxpcstateindex) ? maxpcstateindex : s2start;
      int s2range= (s2start==maxpcstateindex) ? 1 : pcCfg.grainSize;
      CkAssert(s2start<pcCfg.numStates);
      CProxySection_PairCalculator sectProxy =
        CProxySection_PairCalculator::ckNew(pcHandle.pcAID,
            gspaceIndex.y, gspaceIndex.y, 1,
            s1, s1, 1,
            s2start, maxpcstateindex, s2range,
            chunk, chunk, 1);

      CkSectionID sid=sectProxy.ckGetSectionID();
      std::random_shuffle(sid._elems, sid._elems + sid._nElems);
#ifndef _AUTO_DELEGATE_MCASTMGR_ON_
      sectProxy.ckSectionDelegate(mcastGrp);
#endif
      setResultProxy(&sectProxy, false, CkCallback(CkCallback::ignore));
      return sectProxy;
    }

  } // end namespace gspace
} // end namespace cp

#include "RDMAMessages.def.h"
/*@}*/
