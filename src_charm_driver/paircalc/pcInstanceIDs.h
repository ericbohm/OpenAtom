#include "debug_flags.h"
#include "charm++.h"

#ifndef PC_INSTANCE_IDS_H
#define PC_INSTANCE_IDS_H

namespace cp {
  namespace paircalc {

    /// A tiny structure to hold the relevant IDs/ proxies required to interact with a paircalc instance
    class InstanceIDs
    {
      public:
        /// The group providing procNum() for placing the objects of this paircalc array instance
        CkGroupID mapperGID;
        /// The array IDs of the paircalc and its servant input handler arrays
        CkArrayID pcAID, handlerAID;
        /// The CkMulticast group that will handle gspace <--> pc comm
        CkGroupID mCastMgrGID;

        void pup(PUP::er &p)
        {
          p|mapperGID;
          p|pcAID;
          p|handlerAID;
          p|mCastMgrGID;
        }
    };

  } // end namespace paircalc
} // end namespace cp

#endif // PC_INSTANCE_IDS_H

