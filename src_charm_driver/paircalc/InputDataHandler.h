class CProxy_PairCalculator; 	///< Forward declaration required before including the decl.h file
#include "inputDataHandler.decl.h"
#include "../../include/debug_flags.h"	

#ifndef INPUT_DATA_HANDLER_H
#define INPUT_DATA_HANDLER_H

//namespace cp {
//namespace paircalc {

/// Forward declaration of the message type
class paircalcInputMsg;

/** A very thin wrapper class. Instantiated as a chare array that provides GSpace chares with an API to send the input data in.
 * The intent is to permit a somewhat generic message handling behavior that doesnt require GSpace to jump through hoops.
 * 
 * Currently it is desired to collate messages and compute on them later. Requirements might change later to require streamed 
 * computing etc. Perhaps even things like RDMA on the forward path can be folded into this single interface by configuring 
 * this class with different message handlers.
 *   
 * This is a template class only so that it does not have to pay the penalties associated with using base class pointers.
 * The left and right message handlers could simply have been pointers to a base class MessageHandler. However,
 * this would involve virtual function call overhead and would almost certainly prevent inlining of the handlers' operator(). 
 * To avoid this and *yet* retain the ability to support new behavior at a later date, this has been made a template class. 
 * This should provide a relatively stable API for the GSpace chares to use, and all this template muck is only compile 
 * time overhead anyway.
 * 
 * @note: The only behavior expected of handlerType classes is that they act as functors that accept pointers to 
 * messages. i.e. they should define void operator() (paircalcInputMsg *msg). This will be used to pass in the messages 
 * for processing.   
 */
template <class leftHandlerType, class rightHandlerType>
class InputDataHandler: public ArrayElement4D // CBaseT< ArrayElementT<CkIndex4D>,CProxy_InputDataHandler<leftHandlerType,rightHandlerType> > 
{
	public: 
		/// @entry An instance needs to be configured with the actual message handlers. 
		InputDataHandler(CProxy_PairCalculator pcProxy)
		{
			#ifdef DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
				CkPrintf("InputHandler[%d,%d,%d,%d] I was just born. Now I am going to try and get a direct handle to my parent PC\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
			#endif
			/// Get a handle on your customer PairCalc chare array element
			PairCalculator *myCustomerPC = pcProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).ckLocal();
			/// If your customer is on this same processor (which it MUST be) Ask it for pointers to the message handlers
			if (myCustomerPC)
			{
				leftHandler = myCustomerPC->leftHandler();
				rightHandler = myCustomerPC->rightHandler();
				#ifdef DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
					CkPrintf("InputHandler[%d,%d,%d,%d] My left and right handlers: %p %p \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,leftHandler,rightHandler);
				#endif
				
			}
			else
				CkAbort("InputHandler: My customer PairCalc object is not (yet?) on this processor. What am I doing here? I am aborting... \n");
				//CkAbort("\nInputHandler[%d,%d,%d,%d] My customer PairCalc object is not (yet?) on this processor. What am I doing here? I am aborting...",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
		}
		/// @entry migration constructor
		InputDataHandler(CkMigrateMessage *msg) {}
		
		/// @entry Deposit some data belonging to the left matrix block
		inline void acceptLeftData(paircalcInputMsg *msg) 		{ (*leftHandler)(msg); }
		/// @entry Deposit some data belonging to the right matrix block
		inline void acceptRightData(paircalcInputMsg *msg) 		{ (*rightHandler)(msg); } 
		
	private:
		/// Messages are passed on to handler objects that actually handle the data appropriately
		leftHandlerType  *leftHandler;
		rightHandlerType *rightHandler;		
};

//} // end namespace paircalc
//} // end namespace cp

#endif // INPUT_DATA_HANDLER_H


#define CK_TEMPLATES_ONLY
	#include "inputDataHandler.def.h"
#undef CK_TEMPLATES_ONLY
