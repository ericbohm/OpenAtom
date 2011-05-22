#include "pcSectionManager.h"
#include "paircalc/pcConfig.h"

#include <algorithm>

#ifdef USE_COMLIB
extern ComlibInstanceHandle mcastInstanceCP;
extern ComlibInstanceHandle mcastInstanceACP;
#endif

namespace cp {
    namespace ortho {

void PCSectionManager::pup(PUP::er &p) 
{
    p | numPlanes;
    p | numStates;
    p | numChunks;
    p | pcGrainSize;
    p | orthoGrainSize;

    p | pcArrayID;
    p | isSymmetric;
    p | pcSection;

    p | orthoIndex;
    p | orthomCastGrpID;
    p | orthoRedGrpID;
    p | msgPriority;
}



/**
 * The section manager now finds most of its init data from the global config class. If we need support for differently
 * configured PC instances, we should make the section managers init themselves from an instance config object and not 
 * a global config object. But first, we need to implement the concept of a config class for an instance :)
 */
void PCSectionManager::init(const CkIndex2D orthoIdx, const pc::pcConfig &pcCfg, CkArrayID pcAID, CkGroupID oMCastGID, CkGroupID oRedGID)
{
    pcArrayID       = pcAID;
    isSymmetric     = pcCfg.isSymmetric;

    numPlanes       = pcCfg.numPlanes;
    numStates       = pcCfg.numStates;
    numChunks       = pcCfg.numChunks;
    pcGrainSize     = pcCfg.grainSize;
    orthoGrainSize  = pcCfg.orthoGrainSize;

    orthoIndex      = orthoIdx;
    orthomCastGrpID = oMCastGID;
    orthoRedGrpID   = oRedGID;
    msgPriority     = pcCfg.inputMsgPriority;
}




/**
 * Ortho chares are a 2D (nstates x nstates) array. ortho[sx,sy] talks to all the paircalc chares which handle 
 * the ordered pair of states (sx,sy). This will be pc[p,s1,s2,c] where p ranges across all planes and c across 
 * all chunks. 
 */
void PCSectionManager::createPCsection(const int s1, const int s2)
{
    int ecount=0;
    CkArrayIndex4D *elems=new CkArrayIndex4D[numPlanes*numChunks*2];
    for(int chunk = numChunks-1; chunk >=0; chunk--)
        for(int numX = numPlanes-1; numX >=0; numX--)
        {
            CkArrayIndex4D idx4d(numX,s1,s2,chunk);
            elems[ecount++]=idx4d;
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
    
    //    std::random_shuffle(elems, elems + ecount);

    /// Create and save this paircalc section
    pcSection = CProxySection_PairCalculator::ckNew(pcArrayID,  elems, ecount);
    delete [] elems;
}




/**
 * Initialize the planewise section reduction for Ortho sums across all planes and chunks 
 * pass through the the owning Ortho chare so the cookie can be placed in the 2d array
 * (grainSize/orthoGrainSize)^2
 * 
 * Ortho chares talk to paircalc chares based on their state indices. That is, results for an ordered pair of states
 * (s1,s2) from an ortho will end up at the paircalc(s) responsible for that state pair. As paircalcs are decomposed
 * along two other dimensions (planes & points), this results in a section of paircalcs which need to get data for 
 * any given state pair (s1,s2) via a multicast from the ortho responsible for that state pair (s1,s2).
 *
 * The same logic holds for input data from the paircalcs to ortho. As orthos are 2D and oblivious of any plane-wise
 * or point-wise decomposition, all data pertaining to a state pair (s1,s2) from all the paircalcs mut be collated and
 * delivered to the ortho that handles (s1,s2). This results in a reduction across a paircalc section that spans all 
 * planes and chunks of the paircalc array.
 *
 * An extra twist in these straightforward section rules happens because of the presence of phantom chares in symmetric
 * paircalc instance. Phantom chares do not participate in the forward path and hence do not send any input data. Hence
 * they should not be included in the reductions to ortho. However, they work in the backward path off the results from 
 * ortho and are hence included in the multicasts back to the paircalcs. 
 *
 * We compound this a bit further in our greed to avoid extra work when possible. It should be noted that because of the 
 * limitations of the underlying matrix multiply libraries, we canNOT perform a triangular multiply on what is essentially
 * symmetric input. Hence orthos orchestrate a square matrix multiply and end up with almost identical results (differing
 * only by a transpose) in mirror chares across the array diagonal. Without going into the underlying math, we can say
 * that the orthos that talk to phantom or on-diagonal pc chares end up with the exact data required, whereas the ortho 
 * chares that need to talk to the non-phantom sections, have to perform an additional transpose before they can pack the
 * data off to their pc sections. This is where we get lazy and try to avoid this extra transpose if we can. Basically,
 * when phantoms are turned off we rig the section creation so that the orthos which should have spoken to the phantoms
 * instead talk to their non-phantom mirror sections. As these orthos have the data in the correct form already, no one 
 * has to perform any extra transposes.
 *
 * Hence orthos whose indices correspond to those of phantom paircalc chares, will talk to:
 *  - their original phantom section if the user turns on phantoms
 *  - a mirror non-phantom section if the user turns off phantoms 
 *
 */
void PCSectionManager::setupArraySection(CkCallback cb, bool arePhantomsOn, bool useComlibForOrthoToPC)
{
    int s1, s2;
    /// Find the states indices of the paircalcs this ortho *should* be talking to. 
    CkIndex2D pc = computePCStateIndices(orthoIndex.x,orthoIndex.y);
    /// When phantoms are off, ortho chares that should talk to a phantom section will instead talk to  a mirror section
    if (!arePhantomsOn && isSymmetric && pc.y<pc.x)
    {   s1 = pc.y;  s2 = pc.x;  }
    else
    {   s1 = pc.x;  s2 = pc.y;  }

    #ifdef VERBOSE_SECTIONMANAGER
        CkPrintf("Ortho[%d,%d] PCSectionManager setting up a paircalc section that includes PC[%d-%d,%d,%d,%d-%d,%d]\n",orthoIndex.x, orthoIndex.y,0,numPlanes-1,s1,s2,0,numChunks-1,isSymmetric);
    #endif

    /// Create the paircalc section that this ortho chare will actually talk to
    createPCsection(s1,s2);
   
    /// Paircalcs end their forward path by sending data to the orthos. Irrespective of their type (symm/asymm) or 
    /// whether phantoms are turned on or not, they always talk to the orthos corresponding to their state indices.
    /// All the mirrors and section switching tricks that happen above only apply to the multicast back from ortho to
    /// paircalc. Also, phantom paircalcs do not participate in the forward path and do not have any data to send. 
    /// Hence, only orthos whose indices correspond to the non-phantoms will register with their pc sections to get data
    if ( !(isSymmetric && pc.y<pc.x) )
    {
        /// Delegate the pc section --> ortho reduction to CkMulticast
        CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(orthoRedGrpID).ckLocalBranch();
        pcSection.ckSectionDelegate(mcastGrp);
        /// Register this ortho chare with all the paircalcs in this section 
        initGRedMsg *gredMsg=new initGRedMsg;
        gredMsg->cb=cb;
        gredMsg->mCastGrpId= orthoRedGrpID;
        gredMsg->lbsync=false;
        gredMsg->orthoX=orthoIndex.x;
        gredMsg->orthoY=orthoIndex.y;
        pcSection.initGRed(gredMsg);
    }

    /// Delegate the ortho --> pc section multicast to the appropriate library
    if(useComlibForOrthoToPC)
    {
        #ifdef USE_COMLIB
            CkPrintf("NOTE: Rectangular Send In USE\n");
            if(isSymmetric)
                ComlibAssociateProxy(mcastInstanceCP,pcSection);
            else
                ComlibAssociateProxy(mcastInstanceACP,pcSection);
        #endif
    }
    else
    {
        CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(orthomCastGrpID).ckLocalBranch();
        pcSection.ckSectionDelegate(mcastGrp);
    }
}




void PCSectionManager::sendResults(int n, internalType *ptr1, internalType *ptr2, int orthoX, int orthoY, int actionType, int priority)
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
            CkAssert( isfinite(ptr1[i]) );
            CkAssert( isfinite(omsg->matrix1[i]) );
        }
    #endif

    /// Trigger the backward path for my paircalc section
    pcSection.multiplyResult(omsg);
}


void PCSectionManager::sendMatrix(int n, internalType *ptr1, internalType *ptr2, int orthoX, int orthoY, int actionType, int priority)
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
            CkAssert( isfinite(ptr1[i]) );
            CkAssert( isfinite(omsg->matrix1[i]) );
        }
    #endif

    /// Trigger the backward path for my paircalc section
    pcSection.acceptOrthoT(omsg);
}

    } // end namespace ortho
} // end namespace cp
