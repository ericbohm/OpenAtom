#include "load_balance/PeList.h"
#include "load_balance/MapTable.h"
#include "load_balance/IntMap.h"

#ifndef PC_MAP_CONFIG_H
#define PC_MAP_CONFIG_H

namespace cp {
    namespace startup {

/// A container for assorted mapping inputs to pass around easily
struct PCMapConfig
{
    public:
        int boxSize;
        PeListFactory getPeList;
        MapType2 *gSpaceMap;
        bool isTorusMap;
        bool isTorusFake;
        inttriple mapOffset;
};

    } // end namespace startup
} // end namespace cp

#endif // PC_MAP_CONFIG_H

