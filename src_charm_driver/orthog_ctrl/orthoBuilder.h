#include "load_balance/PeList.h"
#include "uber/Uber.h"

#ifndef ORTHO_BUILDER_H
#define ORTHO_BUILDER_H

namespace cp {
    namespace ortho {

class Builder
{
    public:
        /// Construct an ortho world given the configs
        static CkArrayID build(int nstates, PeListFactory getPeList, UberCollection thisInstance);
};

    } // end namespace ortho
} // end namespace cp

#endif // ORTHO_BUILDER_H

