/** \file atomsCache.h
 */
#ifndef PHYSSCRATCHCACHE_H
#define PHYSSCRATCHCACHE_H
#include "src_piny_physics_v1.0/include/class_defs/CP/class_psnonlocal.h"
#include "PhysScratchCache.decl.h"



/** PhysScratchCache class.
 *
 * PhysScratchCache is an entirely passive structure for thread safe
 * data scratch spaces used in for various computations.  In some
 * future, post C++11, world this could be replaced by thread local
 * storage.
 *
 * This cannot be a nodegroup unless the scratch structures are
 * protected from data races.  Some savings in memory footprint could
 * be realized by a pooling scheme to satisfy simultaneous usage, but
 * the effort implementing that has not yet been sufficiently
 * motivated.
 */



class PhysScratchCache: public Group {
 public:
  PSSCRATCH *psscratch;
  PhysScratchCache();
  PhysScratchCache(CkMigrateMessage *m) {}
  ~PhysScratchCache(){}

};


#endif // PHYSSCRATCHCACHE_H
