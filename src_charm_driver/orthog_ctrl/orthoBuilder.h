#include "orthoConfig.h"
#include "paircalc/pcInstanceIDs.h"
#include "load_balance/PeList.h"
#include "uber/Uber.h"

#ifndef ORTHO_BUILDER_H
#define ORTHO_BUILDER_H

namespace cp {
    namespace ortho {

class Builder
{
    public:
        /// A builder will always create the same ortho as per the supplied configs
        Builder(const orthoConfig &_cfg): cfg(_cfg) {}
        /// Construct an ortho world given the configs
        CkArrayID build(int nstates, cp::paircalc::InstanceIDs &asymmHandle, PeListFactory getPeList, UberCollection thisInstance);

    private:
        /// The configurations for the ortho that should be instantiated
        orthoConfig cfg;
};

    } // end namespace ortho
} // end namespace cp

#endif // ORTHO_BUILDER_H

