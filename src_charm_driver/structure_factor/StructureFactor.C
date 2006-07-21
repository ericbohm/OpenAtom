//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file StructureFactor.C
 *
 */
//============================================================================== 

#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "../../include/debug_flags.h"
#include "StructureFactor.h"
#include "StructFactorCache.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"

//==============================================================================

extern CProxy_StructureFactor sfCompProxy;
extern CProxy_StructFactCache sfCacheProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_EnergyGroup egroupProxy; //energy group proxy
extern CProxy_CPcharmParaInfoGrp      scProxy;
StructureFactor::StructureFactor(CkMigrateMessage *m){ }

//#define _CP_DEBUG_SF_CALC_

//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void StructureFactor::computeSF(SFDummyMsg *msg){
//==============================================================================
// 
   int iteration_src = msg->iteration_src;
   delete msg; // prioritized trigger
   CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   if(sim->ees_nloc_on==1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("No structure factors under EES nonlocal\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
   }//endif
   //    CkPrintf("[%d %d %d] compute %d\n",thisIndex.x,thisIndex.y,thisIndex.z,numdest);
   // The guy who called us is up to date. Are we? The caller is one ahead
   // of the energy dude and the atoms because it flipped its iteration counter
   // at the top of the loop.
   if(numdest){
     if(atomsGrpProxy.ckLocalBranch()->iteration != 
        egroupProxy.ckLocalBranch()->iteration_gsp || 
        atomsGrpProxy.ckLocalBranch()->iteration != (iteration_src-1)){
       CkPrintf("Flow of Control Warning  in computeSF : atoms slow\n");
       SFDummyMsg *newMsg = new(8*sizeof(int)) SFDummyMsg;
       CkSetQueueing(newMsg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(newMsg) = config.sfpriority;
       newMsg->iteration_src = iteration_src;
       sfCompProxy(thisIndex.x,thisIndex.y,thisIndex.z).computeSF(newMsg);
     }//endif
   }//endif

//==============================================================================
  if(numdest){ //we have work to do
  //----------------------------------------------------------------------------
  // Manage memory and compute SF
      CkAssert(gsSize>0);
#ifdef _CP_DEBUG_SF_CALC_
      CkPrintf("[%d %d %d] compute\n",thisIndex.x,thisIndex.y,thisIndex.z);
#endif
      // allow for clean slate 
      if(structFactor==NULL) {
	structFactor    = (complex *)fftw_malloc(natm_nl_grp_max*gsSize*sizeof(complex));
      }//endif
      if(structFactor_fx==NULL) {
	structFactor_fx    = (complex *)fftw_malloc(natm_nl_grp_max*gsSize*sizeof(complex));
	structFactor_fy    = (complex *)fftw_malloc(natm_nl_grp_max*gsSize*sizeof(complex));
	structFactor_fz    = (complex *)fftw_malloc(natm_nl_grp_max*gsSize*sizeof(complex));
      }//endif
      AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy
      CPNONLOCAL::CP_calc_Struct_Fact(gsSize,k_x, k_y,k_z, 
				      structFactor,structFactor_fx,structFactor_fy,
				      structFactor_fz,ag->atoms, config.doublePack, 
				      numSfGrps,thisIndex.x);
  //----------------------------------------------------------------------------
  // Communicate the results
      int totalsize=gsSize*natm_nl_grp_max;
      for(int i=0;i<numdest;i++){
	  // create message 
#ifdef _CP_DEBUG_SF_CALC_
	  CkPrintf("[%d %d %d] sending %d to %d\n",thisIndex.x,thisIndex.y,thisIndex.z,i,destinations[i]);
#endif
	  StructFactorMsg *msg = new (totalsize, totalsize, totalsize, totalsize, 8*sizeof(int)) StructFactorMsg;
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.sfpriority+thisIndex.x;
	  msg->datalen = totalsize;
	  msg->gsSize=gsSize;
	  msg->atmGrpIndex = thisIndex.x;
	  msg->planeIndex=thisIndex.y;
	  CmiMemcpy(msg->structFactor,structFactor,totalsize*sizeof(complex));
	  CmiMemcpy(msg->structFactor_fx,structFactor_fx,totalsize*sizeof(complex));
	  CmiMemcpy(msg->structFactor_fy,structFactor_fy,totalsize*sizeof(complex));
	  CmiMemcpy(msg->structFactor_fz,structFactor_fz,totalsize*sizeof(complex));
	  // send message
	  sfCacheProxy[destinations[i]].acceptStructFact(msg);      
      }//endfor
//==============================================================================
// no work to do
  }else{
#ifdef _CP_DEBUG_SF_CALC_
      CkPrintf("[%d %d %d] redundant\n",thisIndex.x,thisIndex.y,thisIndex.z);
#endif
  }//endif

//==============================================================================
   }//end routine
//==============================================================================

