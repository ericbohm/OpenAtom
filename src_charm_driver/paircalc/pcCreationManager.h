#include "pcConfig.h"
#include "paircalcMessages.h"
#include "load_balance/PeList.h"

#ifndef PC_CREATION_MANAGER_H
#define PC_CREATION_MANAGER_H

class MapType2;

namespace cp {
    namespace paircalc {

/**
 * Manages the creation of a complete paircalc bubble that includes
 * two paircalc instances (symmetric and asymmetric), an ortho instance
 * and all accompanying helper entities (map groups, InputDataHandler,
 * OrthoHelper, CLA_Matrix etc.)
 */
class CreationManager
{
    public:
        CreationManager(const pcConfig &_symmCfg, const pcConfig &_asymmCfg);
        void build(CkCallback cb, const int boxSize, PeListFactory getPeList, MapType2 *gSpaceMap);

    private:
        pcConfig symmCfg, asymmCfg;

};

    } // end namespace paircalc
} // end namespace cp

#endif // PC_CREATION_MANAGER_H

