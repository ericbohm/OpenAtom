
#include "charm++.h"
#include "cpaimd.h"
#include "src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "PhysScratchCache.h"
extern CP  readonly_cp;


PhysScratchCache::PhysScratchCache()
{
  // allocate the scratches using the dimensions from CP
  psscratch=new PSSCRATCH(&readonly_cp.cppseudo.nonlocal,  &readonly_cp.cpatom_maps);

}

#include "PhysScratchCache.def.h"
