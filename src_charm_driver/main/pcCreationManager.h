#include "debug_flags.h"
#include "startupMessages.h"
#include "paircalc/pcConfig.h"
#include "orthog_ctrl/orthoConfig.h"
#include "paircalc/pcMapConfig.h"
#include "load_balance/PeList.h"
#include "load_balance/IntMap.h"
#ifndef PC_CREATION_MANAGER_H
#define PC_CREATION_MANAGER_H

namespace cp {
  namespace startup {

    /**
     * Manages the creation of a complete paircalc bubble that includes
     * two paircalc instances (symmetric and asymmetric), an ortho instance
     * and all accompanying helper entities (map groups, InputDataHandler,
     * OrthoHelper, CLA_Matrix etc.)
     */
    class PCCreationManager
    {
      public:
        PCCreationManager(const paircalc::pcConfig &_symmCfg, const paircalc::pcConfig &_asymmCfg, const ortho::orthoConfig &_orthoCfg);
        void build(CkCallback cb, const PCMapConfig &mapCfg, bool lsda);

      private:
        /// The configs for the symmetric and asymmetric paircalc instances
        paircalc::pcConfig symmCfg, asymmCfg;
        /// The configurations for the ortho instance shared by the two pc arrays
        ortho::orthoConfig orthoCfg;

    };

  } // end namespace startup
} // end namespace cp

#endif // PC_CREATION_MANAGER_H

