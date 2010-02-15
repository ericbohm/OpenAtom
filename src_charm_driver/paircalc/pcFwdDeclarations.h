#include "MessageDataCollator.h"

#ifndef PC_FWD_DECLARATIONS_H
#define PC_FWD_DECLARATIONS_H

// Class for paircalc config data
namespace cp { namespace paircalc { class pcConfig; } }
/// A shorter name for the namespace
namespace pc = cp::paircalc;

// The msg carrying input data to the paircalcs
class paircalcInputMsg;
// A template message collator
template <class msgType, typename dataType> class MessageDataCollator;
/// The type of the input msg collator
typedef pc::MessageDataCollator<paircalcInputMsg,double> CollatorType;

class PairCalcID; ///< @note: Should be temporary

// PC chare array proxies
class CProxySection_PairCalculator;

// Input handler chare array proxies
template < class leftHandlerType, class rightHandlerType >  class CProxy_InputDataHandler;
template < class leftHandlerType, class rightHandlerType >  class CProxySection_InputDataHandler;

/// Forward declaration of the RDMA handshake token
struct RDMApair_GSP_PC;

#endif // PC_FWD_DECLARATIONS_H

