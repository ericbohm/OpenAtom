#include "load_balance/IntMap.h"
#include "charm++.h"

#ifndef ORTHO_MAP_H
#define ORTHO_MAP_H
/** @addtogroup mapping
  @{
 */

/// Centroid based ortho map (actual map creation in MapTable.C)
class OrthoMap : public CkArrayMap
{
  private:
    MapType2 *maptable;

  public:
    OrthoMap(MapType2 map)
    {
      maptable = new MapType2(map);
    }

    ~OrthoMap() { }

    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|*maptable;
    }

    inline int procNum(int, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
      int proc;
      proc=maptable->get(index[0],index[1]);
      CkAssert(proc>=0);
      if(numPes != CkNumPes())
        return(proc%CkNumPes());
      else
        return(proc);
    }
};



extern bool fakeTorus; ///< readonly defined in cpaimd.C

/// Map group for placing OrthoHelper chares
class OrthoHelperMap : public CkArrayMap
{
  private:
    MapType2 *maptable;

  public:
    OrthoHelperMap(MapType2 map)
    {
      maptable = new MapType2(map);
    }

    ~OrthoHelperMap() { }

    void pup(PUP::er &p)
    {
      CkArrayMap::pup(p);
      p|*maptable;
    }

    inline int procNum(int, const CkArrayIndex &iIndex)
    {
      int *index=(int *) iIndex.data();
      int proc;
      proc=maptable->get(index[0],index[1]);
      if(fakeTorus)
        return(proc%CkNumPes());
      else
        return(proc);
    }
};
/*@}*/
#endif // ORTHO_MAP_H

