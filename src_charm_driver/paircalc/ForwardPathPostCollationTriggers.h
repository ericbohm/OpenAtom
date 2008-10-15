/** @file This file defines a set of forward path computation trigger functors that are used to start
 * things up after message collation is complete. To be used in tandem with instances of MessageDataCollators
 */ 

#include "ckPairCalculator.h"

#ifndef FORWARD_PATH_POSTCOLLATION_TRIGGERS_H
#define FORWARD_PATH_POSTCOLLATION_TRIGGERS_H

namespace cp {
namespace paircalc {


/** Functor to signal the waiting PC chare that the "left block is ready". 
 * Defined in an attempt at a blind callback to keep MessageDataCollator somewhat generic.
 */
class LeftBlockReadyTrigger
{
	public:
		LeftBlockReadyTrigger(PairCalculator *paircalc): waitingChare(paircalc) {}
		inline void operator() (const double *ptr, const int numMsgs, const int numUnitsInMsg) const 
		{ 
			waitingChare->acceptLeftData(ptr,numMsgs,numUnitsInMsg); 
		}
	private:
		PairCalculator *waitingChare;
};



/** Functor to signal the waiting PC chare that the "right block is ready". 
 * Defined in an attempt at a blind callback to keep MessageDataCollator somewhat generic.
 */
class RightBlockReadyTrigger
{
	public:
		RightBlockReadyTrigger(PairCalculator *paircalc): waitingChare(paircalc) {}
		inline void operator() (const double *ptr, const int numMsgs, const int numUnitsInMsg) const 
		{ 
			waitingChare->acceptRightData(ptr,numMsgs,numUnitsInMsg); 
		}
	private:
		PairCalculator *waitingChare;
};


} // end namespace paircalc
} // end namespace cp

#endif // FORWARD_PATH_POSTCOLLATION_TRIGGERS_H
