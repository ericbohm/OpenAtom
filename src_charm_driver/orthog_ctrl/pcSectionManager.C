#include "pcSectionManager.h"
#include "paircalc/pairCalculator.h" ///< @note: Just for the definition of PairCalcID. Eliminate

extern ComlibInstanceHandle mcastInstanceCP;
extern ComlibInstanceHandle mcastInstanceACP;

namespace cp {
    namespace ortho {

void PCSectionManager::pup(PUP::er &p) 
{
    p | pcSection;
    p | numStates;
    p | pcGrainSize;
    p | orthoGrainSize;
    p | numChunks;
    p | pcArrayID;
    p | isSymmetric;
    p | arePhantomsOn; 
    p | orthoIndex;
    p | orthomCastGrpID;
    p | orthoRedGrpID;
}




void PCSectionManager::init(const CkIndex2D orthoIdx, const PairCalcID &pcid,const int orthoGrSize)
{
    numStates       = pcid.nstates;
    numChunks       = pcid.numChunks;
    pcGrainSize     = pcid.GrainSize;
    orthoGrainSize  = orthoGrSize;

    pcArrayID       = pcid.Aid;
    isSymmetric     = pcid.Symmetric;
    //arePhantomsOn   = phantom;

    orthoIndex      = orthoIdx;
    orthomCastGrpID = pcid.orthomCastGrpId;
    orthoRedGrpID   = pcid.orthoRedGrpId;
}




/**
 * Initialize the planewise section reduction for Ortho sums across all planes and chunks 
 * pass through the the owning Ortho chare so the cookie can be placed in the 2d array
 * (grainSize/orthoGrainSize)^2
 */
void PCSectionManager::setupArraySection(int numZ, int* z, CkCallback cb, CkCallback synccb, int s1, int s2, bool arePhantomsOn, bool direct, bool commlib)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("PCSectionManager::setupArraySection called\n");
    #endif
    int ecount=0;
    CkArrayIndex4D *elems=new CkArrayIndex4D[numZ*numChunks*2];
    //add chunk loop
    for(int chunk = numChunks-1; chunk >=0; chunk--)
    {
        for(int numX = numZ-1; numX >=0; numX--)
        {
            #ifdef VERBOSE_SECTIONMANAGER
                CkPrintf("initOneRedSect for s1 %d s2 %d ortho %d %d sym %d plane %d\n",s1,s2,orthoIndex.x, orthoIndex.y,isSymmetric, numX);
            #endif
            if(arePhantomsOn && s1!=s2)
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
    int numOrthoCol= pcGrainSize / orthoGrainSize;
    int maxorthostateindex=(numStates / orthoGrainSize - 1) * orthoGrainSize;
    int orthoIndexX=(orthoIndex.x * orthoGrainSize);
    
    orthoIndexX= (orthoIndexX>maxorthostateindex) ? maxorthostateindex : orthoIndexX;
    int orthoIndexY=(orthoIndex.y * orthoGrainSize);
    orthoIndexY= (orthoIndexY>maxorthostateindex) ? maxorthostateindex : orthoIndexY;
    orthoIndexX-=s1;
    orthoIndexY-=s2;
    int orthoArrIndex=orthoIndexX*numOrthoCol+orthoIndexY;
    
    int newListStart=orthoArrIndex;
    if(newListStart> ecount)
        newListStart= newListStart % ecount;
    bool order=reorder_elem_list_4D( elems, ecount, newListStart);
    CkAssert(order);
    // now that we have the section, make the proxy
    CProxySection_PairCalculator sProxy=CProxySection_PairCalculator::ckNew(pcArrayID,  elems, ecount);
    CProxySection_PairCalculator *sectProxy=&sProxy;
    delete [] elems;
    
    // and do delegation
    pcSection = sProxy;
    
    if(!arePhantomsOn && !direct) // only delegating nonphantom mcast proxy for reduction
    {
        CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(orthoRedGrpID).ckLocalBranch();
        sectProxy->ckSectionDelegate(mcastGrp);
        // Send the multcast message to initialize the ortho section tree and set the cookie
        initGRedMsg *gredMsg=new initGRedMsg;
        gredMsg->cb=cb;
        gredMsg->mCastGrpId= orthoRedGrpID;
        gredMsg->lbsync=false;
        gredMsg->synccb=synccb;
        gredMsg->orthoX=orthoIndex.x;
        gredMsg->orthoY=orthoIndex.y;
        sectProxy->initGRed(gredMsg);
    }
    else
    {
        if(commlib)
        {
            CkPrintf("NOTE: Rectangular Send In USE\n");
            if(isSymmetric)
                ComlibAssociateProxy(&mcastInstanceCP,*sectProxy);
            else
                ComlibAssociateProxy(&mcastInstanceACP,*sectProxy);
            /*
            if(!isSymmetric)
                ComlibAssociateProxy(&mcastInstanceACP,*sectProxy);
            */
        }
        else
        {
            //CkPrintf("PC: proxy without commlib\n");
            CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(orthomCastGrpID).ckLocalBranch();
            sectProxy->ckSectionDelegate(mcastGrp);
        }
    }
    //  return *sectProxy;
}




void PCSectionManager::sendResults(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("PCSectionManager::sendResults()\n");
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


void PCSectionManager::sendMatrix(int n, double *ptr1, double *ptr2, int orthoX, int orthoY, int actionType, int priority)
{
    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("PCSectionManager::sendMatrix()\n");
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

    } // end namespace ortho
} // end namespace cp
