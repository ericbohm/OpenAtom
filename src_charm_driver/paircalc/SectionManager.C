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

    } // end namespace paircalc
} // end namespace cp
