#include "SectionManager.h"

namespace cp {
    namespace paircalc {

void SectionManager::pup(PUP::er &p) 
{
    p | pcSection;
}




void SectionManager::setupArraySection(const InstanceInfo &info, int numZ, int* z, CkCallback cb, CkCallback synccb, int s1, int s2, bool direct, bool commlib)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("SectionManager::setupArraySection called\n");
    #endif
    int ecount=0;
    CkArrayIndex4D *elems=new CkArrayIndex4D[numZ*info.numChunks*2];
    //add chunk loop
    for(int chunk = info.numChunks-1; chunk >=0; chunk--)
    {
        for(int numX = numZ-1; numX >=0; numX--)
        {
            #ifdef VERBOSE_SECTIONMANAGER
                CkPrintf("initOneRedSect for s1 %d s2 %d ortho %d %d sym %d plane %d\n",s1,s2,orthoIndex.x, orthoIndex.y,info.isSymmetric, numX);
            #endif
            if(info.arePhantomsOn && s1!=s2)
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
    int numOrthoCol= info.pcGrainSize / info.orthoGrainSize;
    int maxorthostateindex=(info.numStates / info.orthoGrainSize - 1) * info.orthoGrainSize;
    int orthoIndexX=(info.orthoIndex.x * info.orthoGrainSize);
    
    orthoIndexX= (orthoIndexX>maxorthostateindex) ? maxorthostateindex : orthoIndexX;
    int orthoIndexY=(info.orthoIndex.y * info.orthoGrainSize);
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
    CProxySection_PairCalculator sProxy=CProxySection_PairCalculator::ckNew(info.pcArrayID,  elems, ecount);
    CProxySection_PairCalculator *sectProxy=&sProxy;
    delete [] elems;
    
    // and do delegation
    pcSection = sProxy;
    
    if(!info.arePhantomsOn && !direct) // only delegating nonphantom mcast proxy for reduction
    {
        CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(info.orthoRedGrpID).ckLocalBranch();
        sectProxy->ckSectionDelegate(mcastGrp);
        // Send the multcast message to initialize the ortho section tree and set the cookie
        initGRedMsg *gredMsg=new initGRedMsg;
        gredMsg->cb=cb;
        gredMsg->mCastGrpId= info.orthoRedGrpID;
        gredMsg->lbsync=false;
        gredMsg->synccb=synccb;
        gredMsg->orthoX=info.orthoIndex.x;
        gredMsg->orthoY=info.orthoIndex.y;
        sectProxy->initGRed(gredMsg);
    }
    else
    {
        if(commlib)
        {
            CkPrintf("NOTE: Rectangular Send In USE\n");
            if(info.isSymmetric)
                ComlibAssociateProxy(&mcastInstanceCP,*sectProxy);
            else
                ComlibAssociateProxy(&mcastInstanceACP,*sectProxy);
            /*
            if(!info.isSymmetric)
                ComlibAssociateProxy(&mcastInstanceACP,*sectProxy);
            */
        }
        else
        {
            //CkPrintf("PC: proxy without commlib\n");
            CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(info.orthomCastGrpID).ckLocalBranch();
            sectProxy->ckSectionDelegate(mcastGrp);
        }
    }
    //  return *sectProxy;
}




void SectionManager::sendResults(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority)
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


void SectionManager::sendMatrix(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("SectionManager::sendMatrix()\n");
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
    pcSection.acceptOrthoT(omsg);
}

    } // end namespace paircalc
} // end namespace cp
