#include "debug_flags.h"
#include "MessageDataCollator.h"
#include "ckcomplex.h"

#ifndef PC_FWD_DECLARATIONS_H
#define PC_FWD_DECLARATIONS_H

// Class for paircalc config data
namespace cp { namespace paircalc { class pcConfig; } }
/// A shorter name for the namespace
namespace pc = cp::paircalc;


/// The type of input data as perceived by the paircalc world
typedef complex inputType;
/// Compile time decision on whether paircalc crunches complex numbers or real
//#define CP_PAIRCALC_USES_COMPLEX_MATH
#ifdef CP_PAIRCALC_USES_COMPLEX_MATH
/// The representation of the input data internal to the paircalc world
typedef complex internalType;
#else
/// The representation of the input data internal to the paircalc world
typedef double internalType;
#endif
// The ratio between the input and internal data type sizes
#define pcDataSizeFactor ( sizeof(inputType) / sizeof(internalType) )


// The msg carrying input data to the paircalcs
class paircalcInputMsg;
// A template message collator
template <class msgType, typename dataType> class MessageDataCollator;
/// The type of the input msg collator
typedef pc::MessageDataCollator<paircalcInputMsg,inputType> CollatorType;


// PC chare array proxies
class CProxySection_PairCalculator;

// Input handler chare array proxies
template < class leftHandlerType, class rightHandlerType >  class CProxy_InputDataHandler;
template < class leftHandlerType, class rightHandlerType >  class CProxySection_InputDataHandler;

/// Forward declaration of the RDMA handshake token
struct RDMApair_GSP_PC;

#endif // PC_FWD_DECLARATIONS_H

