#include "load_balance/PeList.h"
#include "load_balance/MapTable.h"
#include "load_balance/IntMap.h"

#ifndef PC_MAP_CONFIG_H
#define PC_MAP_CONFIG_H
/** @addtogroup PairCalculator
    @{
*/

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
  PCMapConfig(){ CkAbort("why are we here?");}
  PCMapConfig(int _boxSize, PeListFactory _getPeList, MapType2 *_gSpaceMap, bool _isTorusMap, bool _isTorusFake, inttriple _mapOffset): boxSize(_boxSize), getPeList(_getPeList), gSpaceMap(_gSpaceMap), isTorusMap(_isTorusMap), isTorusFake(_isTorusFake), mapOffset(_mapOffset)
  {}
};

    } // end namespace startup
} // end namespace cp
/*@}*/
#endif // PC_MAP_CONFIG_H

