#include "load_balance/PeList.h"
#include "uber/Uber.h"

#ifndef ORTHO_ARRAY_BUILDER_H
#define ORTHO_ARRAY_BUILDER_H

namespace cp {
    namespace ortho {

class ArrayBuilder
{
    public:
        static void build(int nstates, PeListFactory getPeList, UberCollection thisInstance);
};

    } // end namespace ortho
} // end namespace cp

#endif // ORTHO_ARRAY_BUILDER_H

