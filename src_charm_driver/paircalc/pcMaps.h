#include "load_balance/IntMap.h"
#include "pcMaps.decl.h"
#include "charm++.h"

extern CkVec <MapType4> AsymScalcImaptable;
extern CkVec <MapType4> SymScalcImaptable;

#ifndef PC_MAPS_H
#define PC_MAPS_H

/// Class used for instantiation of pair-calculator group objects.
class SCalcMap : public CkArrayMap
{
    private:
        UberCollection thisInstance;
        MapType4 *maptable;
        bool symmetric;

    public:
        SCalcMap(bool _symmetric, UberCollection _instance): symmetric(_symmetric), thisInstance(_instance)
        {
            if(symmetric)
                maptable= &SymScalcImaptable[thisInstance.getPO()];
            else
                maptable= &AsymScalcImaptable[thisInstance.getPO()];
        }

        void pup(PUP::er &p)
        {
            CkArrayMap::pup(p);
            p|symmetric;
            p|thisInstance;
            if(symmetric)
                maptable= &SymScalcImaptable[thisInstance.getPO()];
            else
                maptable= &AsymScalcImaptable[thisInstance.getPO()];
        }

        inline int procNum(int, const CkArrayIndex &iIndex)
        {
            int proc;
            short *sindex=(short *) iIndex.data();
            proc=maptable->get(sindex[0], sindex[1], sindex[2], sindex[3]);
            CkAssert(proc>=0);
            if(numPes!=CkNumPes())
                return(proc%CkNumPes());
            else
                return(proc);
        }
};

#endif // PC_MAPS_H

