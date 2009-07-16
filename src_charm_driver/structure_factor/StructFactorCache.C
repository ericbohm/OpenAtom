//==============================================================================
/** \file StructFactorCache.C
 *
 */
//==============================================================================
 
#include "utility/util.h"
#include "main/groups.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "StructFactorCache.h"
#include "gParticlePlane.decl.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
extern CkVec <CProxy_CP_State_ParticlePlane> UparticlePlaneProxy;
extern Config config;

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int StructFactCache::existStructFact(int planeIndex) {
    for(int i=0; i<planeCountList.length(); i++){
	PlaneCount thisCount = planeCountList[i];
	if(thisCount.plane == planeIndex) {
	    return i;
	}
    }
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("[%d] lacks plane %d\n",CkMyPe(),planeIndex);
#endif

    return -1;
}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int StructFactCache::existStructFactGrp(int planeIndex, int atmGrp) {
    for(int i=0; i<planeCountList.length(); i++){
	PlaneCount thisCount = planeCountList[i];
	if(thisCount.plane == planeIndex) {
	  if(structFactorAtmGrps[i][atmGrp]==1)
	    {
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("[%d] has plane %d atom %d\n",CkMyPe(),planeIndex,atmGrp);
#endif
	      return i;
	    }
	}
    }
#ifdef _CP_DEBUG_SF_CACHE_
    CkPrintf("[%d] lacks plane %d atom %d\n",CkMyPe(),planeIndex,atmGrp);
#endif
    return -1;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int StructFactCache::incCountStructFact(int planeIndex) {
    
    for(int i=0; i<planeCountList.length(); i++){
	PlaneCount* thisCount = &(planeCountList[i]);
	if(thisCount->plane == planeIndex) {
	    thisCount->count ++;
	    return thisCount->count;
	}
    }
 
    planeCountList.insert(planeCountList.length(), PlaneCount(planeIndex, 1, 1));
    return 1;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int StructFactCache::decCountStructFact(int planeIndex) {
    for(int i=0; i<planeCountList.length(); i++){
	PlaneCount* thisCount = &(planeCountList[i]);
	if(thisCount->plane == planeIndex) {
	    thisCount->count--;
	    if(thisCount->count < 0) {
		thisCount->count = 0;
		CkPrintf("Over Decrementing counters!!!\n");
	    }
	    return thisCount->count;
	}
    }
    CkPrintf("Trying to Decrement non-existing counters!!!\n");
    return -1;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void StructFactCache::printCountStructFact() {
    CkPrintf("=====================================\n");
    for(int i=0; i<planeCountList.length(); i++){
	PlaneCount thisCount = planeCountList[i];
	CkPrintf(" plane -- %d,  count -- %d\n", thisCount.plane, thisCount.count);
    }    
    CkPrintf("=====================================\n");
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void StructFactCache::setZero(int  planeIndex) {
    int i = existStructFact(planeIndex);
    if(i>=0)
      {
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("[%d] zeros plane %d\n",CkMyPe(),planeIndex);
#endif
	planeCountList[i].updated = 0;
	for(int j=0;j<numSfGrps;j++)
	  structFactorAtmGrps[i][j]=0;
      }
}
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int StructFactCache::getStructFact(int planeIndex, int atmGrpIndex, complex** sf, 
                      complex** sf_x, complex** sf_y, complex** sf_z) {
    int i = existStructFact(planeIndex);
    if(i < 0) return -1;
    int arroffset=structFactorSize[i]*atmGrpIndex*natm_nl_grp_max;
    *sf = structFactorList[i]+arroffset;
    *sf_x = structFactorfxList[i]+arroffset;
    *sf_y = structFactorfyList[i]+arroffset;
    *sf_z = structFactorfzList[i]+arroffset;
    return i;
}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void StructFactCache::getStructFactIdx(int sfindex, int atmGrpIndex, complex** sf, 
                      complex** sf_x, complex** sf_y, complex** sf_z) {
    int arroffset=structFactorSize[sfindex]*atmGrpIndex*natm_nl_grp_max;
    *sf = structFactorList[sfindex]+arroffset;
    *sf_x = structFactorfxList[sfindex]+arroffset;
    *sf_y = structFactorfyList[sfindex]+arroffset;
    *sf_z = structFactorfzList[sfindex]+arroffset;
}
//==============================================================================


void StructFactCache::acceptStructFact(StructFactorMsg *msg)
{
    // check if we have the plane

    int atmIndex=msg->atmGrpIndex;  
    int planeIndex=msg->planeIndex;
    int arroffset=msg->gsSize*atmIndex*natm_nl_grp_max;
    int sfindex=existStructFact(planeIndex);
#ifdef _CP_DEBUG_SF_CACHE_
    CkPrintf("[%d] received SF for plane %d atom %d\n",CkMyPe(),planeIndex, atmIndex);
#endif

    if(sfindex<0)
    { // create the entry
      //allocate the arrays
	incCountStructFact(planeIndex);
	int structSize=msg->gsSize*numSfGrps*natm_nl_grp_max;
	complex *structFactor    = (complex *)fftw_malloc(structSize*sizeof(complex));
	complex *structFactor_fx = (complex *)fftw_malloc(structSize*sizeof(complex));
	complex *structFactor_fy = (complex *)fftw_malloc(structSize*sizeof(complex));
	complex *structFactor_fz = (complex *)fftw_malloc(structSize*sizeof(complex));
	int *atmGrpArrive        = new int[numSfGrps];
	bzero(atmGrpArrive, numSfGrps*sizeof(int));
	structFactorAtmGrps.insert(structFactorAtmGrps.length(),atmGrpArrive);
	structFactorList.insert(structFactorList.length(),structFactor);
	structFactorfxList.insert(structFactorfxList.length(), structFactor_fx);
	structFactorfyList.insert(structFactorfyList.length(), structFactor_fy);
	structFactorfzList.insert(structFactorfzList.length(), structFactor_fz);
	structFactorSize.insert(structFactorSize.length(), msg->gsSize);
	//set list index
	sfindex=existStructFact(planeIndex);
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("[%d] created SF entry %d for plane %d atom %d\n",CkMyPe(),sfindex, planeIndex, atmIndex);
#endif
    }
    CkAssert(sfindex>=0);
    planeCountList[sfindex].updated      =  1; // possible semantic misuse here
    structFactorAtmGrps[sfindex][atmIndex] = 1;
    CmiMemcpy(structFactorList[sfindex]+arroffset,msg->structFactor,msg->datalen*sizeof(complex));
    CmiMemcpy(structFactorfxList[sfindex]+arroffset,msg->structFactor_fx,msg->datalen*sizeof(complex));
    CmiMemcpy(structFactorfyList[sfindex]+arroffset,msg->structFactor_fy,msg->datalen*sizeof(complex));
    CmiMemcpy(structFactorfzList[sfindex]+arroffset,msg->structFactor_fz,msg->datalen*sizeof(complex));

    // launch the computeZ's registered for this atom index
    int ppregindex;
    if((ppregindex=existsPP(planeIndex,atmIndex))>=0)
    {
	PlaneAtom regPPs=ppList[ppregindex];
	CkArrayIndex2D idx2d;
	for(int i=0;i<regPPs.particles.length();i++)
	{

	  // we could set up a section multicast for this, but they're
	  // all processor local so it would only save us some copies
	  // of a tiny message.  precious little savings
	  
	  // these are of course not made as cklocal branch calls
	  // because we want to lower the priority of this work so it
	  // happens when the system is otherwise idle

	    idx2d=regPPs.particles[i];
	    // call up computeZ
	    PPDummyMsg *pmsg = new (8*sizeof(int)) PPDummyMsg;
	    pmsg->atmGrp=atmIndex;
	    pmsg->sfindex=sfindex;
	    CkSetQueueing(pmsg, CK_QUEUEING_IFIFO);
	    *(int*)CkPriorityPtr(pmsg) = config.sfpriority+atmIndex+numSfGrps; //lower than sf and sfcache
	    UparticlePlaneProxy[thisInstance.proxyOffset](idx2d.index[0], idx2d.index[1]).computeZ(pmsg);
	    //	    CkPrintf("triggering computeZ for %d on %d %d\n", pmsg->atmGrp, idx2d.index[0], idx2d.index[1]);
	}
    }
    else
    {
	char out[100];
	snprintf(out, 100, "[%d] accepted SF for %d %d without any registered particles",CkMyPe(), planeIndex, atmIndex);
	CkAbort(out);
    }
    delete msg; //do not delete nokeep
}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** local particle planes register themselves with the cache so they can be launched by the arrival of 
 * an update for the (plane, atom)
 */
int StructFactCache::registerPP(int state, int plane, int atom) {
  // check for existing entry
  PlaneAtom pa(plane,atom);
  CkArrayIndex2D idx2d(state,plane);
  for(int i=0; i<ppList.length(); i++){
    if(ppList[i]==pa)
      {
	ppList[i].particles.push_back(idx2d);
#ifdef _CP_DEBUG_SF_CACHE_
	CkPrintf("[%d] registers PP[%d,%d] atom %d at %d\n",CkMyPe(),state,plane,atom, ppList[i].particles.length()-1);
#endif
	return i;
      }
  }
  // no entry add new one
  ppList.push_back(pa);
  ppList[ppList.length()-1].particles.push_back(idx2d);
#ifdef _CP_DEBUG_SF_CACHE_
  CkPrintf("SFC [%d] registers PP[%d,%d] atom %d at %d\n",CkMyPe(),state,plane,atom,ppList[ppList.length()-1].particles.length()-1);
#endif
  return ppList.length()-1;
}
//==============================================================================

int StructFactCache::existsPP(int plane, int atom) {
  // find existing entry
  PlaneAtom pa(plane,atom);
  for(int i=0; i<ppList.length(); i++)
    if(ppList[i]==pa)
      return i;
  return -1;
}

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void StructFactCache::removeAll() {

    for(int i=structFactorList.length()-1; i>=0; i--){
	fftw_free(structFactorList[i]);
	fftw_free(structFactorfxList[i]);
	fftw_free(structFactorfyList[i]);
	fftw_free(structFactorfzList[i]);
	delete []structFactorAtmGrps[i];
    }
    
    structFactorList.removeAll();
    structFactorfxList.removeAll();
    structFactorfyList.removeAll();
    structFactorfzList.removeAll();
    structFactorSize.removeAll();
    gSpaceSlabs.removeAll();
}
//==============================================================================

#include "structureFactorCache.def.h"

