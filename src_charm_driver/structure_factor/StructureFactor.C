/** \file StructureFactor.C
 *
 */
 
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "../../include/debug_flags.h"
#include "StructureFactor.h"
#include "StructFactorCache.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"

extern CProxy_StructFactCache sfCacheProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
StructureFactor::StructureFactor(CkMigrateMessage *m){ }

void StructureFactor::computeSF(SFDummyMsg *msg)
{
  delete msg; // prioritized trigger

  if(numdest)
    { //we have work to do
      CkAssert(gsSize>0);
#ifdef _CP_DEBUG_SF_CALC_
      CkPrintf("[%d %d %d] compute\n",thisIndex.x,thisIndex.y,thisIndex.z);
#endif
      // allow for clean slate 
      if(structFactor==NULL) {
	structFactor    = new complex[natm_nl_grp_max*gsSize];
      }
      if(structFactor_fx==NULL) {
	structFactor_fx    = new complex[natm_nl_grp_max*gsSize];
	structFactor_fy    = new complex[natm_nl_grp_max*gsSize];
	structFactor_fz    = new complex[natm_nl_grp_max*gsSize];
      }

      // compute
      AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy


      CPNONLOCAL::CP_calc_Struct_Fact(gsSize, 
				      k_x, k_y,k_z, 
				      structFactor,structFactor_fx,structFactor_fy,
				      structFactor_fz,ag->atoms, config.doublePack, 
				      numSfGrps, thisIndex.x);
      int totalsize=gsSize*natm_nl_grp_max;
      int totalbytesize=gsSize*natm_nl_grp_max*sizeof(complex);
      for(int i=0;i<numdest;i++)
	{
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
	}
    }
  else
    {
#ifdef _CP_DEBUG_SF_CALC_
      CkPrintf("[%d %d %d] redundant\n",thisIndex.x,thisIndex.y,thisIndex.z);
#endif
    }
}
