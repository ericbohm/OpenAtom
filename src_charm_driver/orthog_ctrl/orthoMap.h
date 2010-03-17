#include "main/cpaimd.h"

#ifndef ORTHO_MAP_H
#define ORTHO_MAP_H

#define USE_INT_MAP

/// Centroid based ortho map (actual map creation in MapTable.C)
class OrthoMap : public CkArrayMap
{
    private:
        MapType2 *maptable;

    public:
        OrthoMap(MapType2 map)
        {
            #ifdef USE_INT_MAP
                maptable = new MapType2(map);
            #else
                maptable= &Orthomaptable;
            #endif
        }

        ~OrthoMap() { }

        void pup(PUP::er &p)
        {
            CkArrayMap::pup(p);
            #ifdef USE_INT_MAP
                p|*maptable;
            #else
                maptable= &Orthomaptable;
            #endif
        }

        inline int procNum(int, const CkArrayIndex &iIndex)
        {
            int *index=(int *) iIndex.data();
            int proc;
            #ifdef USE_INT_MAP
                proc=maptable->get(index[0],index[1]);
            #else
                proc=maptable->get(intdual(index[0],index[1]));
            #endif
            CkAssert(proc>=0);
            if(numPes != CkNumPes())
                return(proc%CkNumPes());
            else
                return(proc);
        }
};




/// Map group for placing OrthoHelper chares
class OrthoHelperMap : public CkArrayMap
{
    private:
        MapType2 *maptable;

    public:
        OrthoHelperMap(MapType2 map)
        {
            #ifdef USE_INT_MAP
                maptable = new MapType2(map);
            #else
                maptable= &OrthoHelpermaptable;
            #endif
        }

        ~OrthoHelperMap() { }

        void pup(PUP::er &p)
        {
            CkArrayMap::pup(p);
            #ifdef USE_INT_MAP
                p|*maptable;
            #else
                maptable= &OrthoHelpermaptable;
            #endif
        }

        inline int procNum(int, const CkArrayIndex &iIndex)
        {
            int *index=(int *) iIndex.data();
            int proc;
            #ifdef USE_INT_MAP
                proc=maptable->get(index[0],index[1]);
            #else
                proc=maptable->get(intdual(index[0],index[1]));
            #endif
            if(fakeTorus)
                return(proc%CkNumPes());
            else
                return(proc);
        }
};

#endif // ORTHO_MAP_H

