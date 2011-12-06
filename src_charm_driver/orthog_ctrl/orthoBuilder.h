#include "debug_flags.h"
#include "orthoConfig.h"
#include "paircalc/pcInstanceIDs.h"
#include "paircalc/pcMapConfig.h"
#include "load_balance/PeList.h"

#ifndef ORTHO_BUILDER_H
#define ORTHO_BUILDER_H

namespace cp {
    namespace ortho {

/// A class that orchestrates the mapping and creation of one ortho array and accompanying chares like OrthoHelper, CLA_Matrix etc.
class Builder
{
    public:
        /// A builder will always create the same ortho as per the supplied configs
        Builder(const orthoConfig &_cfg): cfg(_cfg) {}
        /// Construct an ortho world given the configs
        CkArrayID build(cp::paircalc::InstanceIDs &asymmHandle, const startup::PCMapConfig mapCfg);

    private:
        /// The configurations for the ortho that should be instantiated
        orthoConfig cfg;
};

    } // end namespace ortho
} // end namespace cp

#endif // ORTHO_BUILDER_H

