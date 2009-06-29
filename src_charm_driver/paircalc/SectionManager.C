#include "SectionManager.h"

namespace cp {
    namespace paircalc {

void SectionManager::setupArraySection(int numZ, int* z, int numChunks,  PairCalcID* pcid, CkCallback cb, CkCallback synccb, int s1, int s2, int orthoX, int orthoY, int orthoGrainSize, bool phantom, bool direct, bool commlib)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("SectionManager::setupArraySection called\n");
    #endif
    int ecount=0;
    CkArrayIndex4D *elems=new CkArrayIndex4D[numZ*numChunks*2];
    //add chunk loop
    for(int chunk = numChunks-1; chunk >=0; chunk--)
    {
        for(int numX = numZ-1; numX >=0; numX--)
        {
            #ifdef VERBOSE_SECTIONMANAGER
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
    pcSection = sProxy;
    
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




void SectionManager::sendResults(int n, double *ptr1, double *ptr2, PairCalcID *pcid, int orthoX, int orthoY, int actionType, int priority)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("SectionManager::sendResults()\n");
    #endif

    /// Allocate a msg of the right size
    multiplyResultMsg *omsg;
    int size2 = (ptr2)? n : 0;
    if(priority>0)
    {
        omsg=new (n, size2, 8*sizeof(int) ) multiplyResultMsg;
        *(int*)CkPriorityPtr(omsg) = priority;
        CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
    }
    else
        omsg=new (n, size2) multiplyResultMsg;

    /// Fill it with results
    if(ptr2==NULL)
        omsg->init1(n, ptr1, orthoX, orthoY, actionType);
    else 
        omsg->init(n, n, ptr1, ptr2, orthoX, orthoY, actionType);
    #ifdef _NAN_CHECK_
        for(int i=0;i<n ;i++)
        {
            CkAssert(finite(ptr1[i]));
            CkAssert(finite(omsg->matrix1[i]));
        }
    #endif

    /// Trigger the backward path for my paircalc section
    pcSection.multiplyResult(omsg);
}

    } // end namespace paircalc
} // end namespace cp
