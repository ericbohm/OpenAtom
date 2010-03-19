#include "pcMaps.decl.h"
#include "charm++.h"

#ifndef PC_MAPS_H
#define PC_MAPS_H

/// Class used for instantiation of pair-calculator group objects.
class SCalcMap : public CkArrayMap
{
    private:
        MapType4 maptable;

    public:
        SCalcMap(const MapType4 _mtable): maptable(_mtable) { }

        void pup(PUP::er &p)
        {
            CkArrayMap::pup(p);
            p|maptable;
        }

        inline int procNum(int, const CkArrayIndex &iIndex)
        {
            int proc;
            short *sindex=(short *) iIndex.data();
            proc = maptable.get(sindex[0], sindex[1], sindex[2], sindex[3]);
            CkAssert(proc>=0);
            if(numPes!=CkNumPes())
                return(proc%CkNumPes());
            else
                return(proc);
        }

        /// Let local code access the paircalc maptable I store. (for use by Ortho mapping code)
        MapType4* getMapTable() { return &maptable; }
};

#endif // PC_MAPS_H

