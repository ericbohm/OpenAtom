#include "debug_flags.h"
#include "pcConfig.h"
#include "pcInstanceIDs.h"
#include "load_balance/PeList.h"
#include "uber/Uber.h"

#include "charm++.h"

#ifndef PC_BUILDER_H
#define PC_BUILDER_H

namespace cp {
    namespace paircalc {

/// A class that orchestrates the mapping and creation of one paircalc array and accompanying classes
class Builder
{
    public:
        /// Always charge a builder with creating one paircalc array
        Builder(const pcConfig &_cfg): cfg(_cfg) {}
        /// Trigger the creation of a pc array with the given configs, within the given pes/boxes etc
        InstanceIDs build(const int boxSize, PeListFactory getPeList, UberCollection thisInstance);

    private:
        /// Create the mapping required to instantiate a PC array
        void createMap(const int boxSize, PeListFactory getPeList, UberCollection thisInstance);
        /// Create a paircalc array using info in the supplied pcConfig object
        void createPairCalcs();


        /// The configs for the paircalc array that I am charged with building
        const pcConfig cfg;
        /// The result of an array build
        InstanceIDs pcHandle;
};

    } // end namespace paircalc
} // end namespace cp

#endif // PC_BUILDER_H

