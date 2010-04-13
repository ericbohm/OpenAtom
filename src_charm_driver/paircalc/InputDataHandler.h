#include "debug_flags.h"
#include "inputDataHandler.decl.h"
#include "RDMAMessages.h"
#include "ckPairCalculator.h"


#ifndef INPUT_DATA_HANDLER_H
#define INPUT_DATA_HANDLER_H

//namespace cp {
//namespace paircalc {

/// Forward declaration of the message types
class paircalcInputMsg;

/** A very thin wrapper class. Instantiated as a chare array that provides GSpace chares with an API to send the input data in.
 * The intent is to permit a somewhat generic message handling behavior that doesnt require GSpace to jump through hoops.
 * 
 * Currently it is desired to collate messages and compute on them later. Requirements might change later to require streamed 
 * computing etc which can be met by configuring this class with different message handlers.
 *   
 * This is a template class only so that it does not have to pay the penalties associated with using base class pointers.
 * The left and right message handlers could simply have been pointers to a base class MessageHandler. However,
 * this would involve virtual function call overhead and would almost certainly prevent inlining of the handlers' operator(). 
 * To avoid this and *yet* retain the ability to support new behavior at a later date, this has been made a template class. 
 * This should provide a relatively stable API for the GSpace chares to use, and all this template muck is only compile 
 * time overhead anyway.
 * 
 * @note: Behavior expected of handlerType classes:
 *  - void handlerType::operator() (paircalcInputMsg *msg)   : Used to pass in the messages for processing
 *  - void handlerType::setupRDMA(RDMASetupRequestMsg<RDMApair_GSP_PC> *msg,RDMApair_GSP_PC *tkn) 	: Used to pass on a RDMA setup request from a sender
 */
template <class leftHandlerType, class rightHandlerType>
class InputDataHandler: public ArrayElement4D // CBaseT< ArrayElementT<CkIndex4D>,CProxy_InputDataHandler<leftHandlerType,rightHandlerType> > 
{
	public:
		// ----------------------------------------- Cons/Des-truction -----------------------------------------
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
		
		// ----------------------------------------- Entry Methods -----------------------------------------
		/// @entry Deposit some data belonging to the left matrix block
		inline void acceptLeftData(paircalcInputMsg *msg) 		{ (*leftHandler)(msg); }
		/// @entry Deposit some data belonging to the right matrix block
		inline void acceptRightData(paircalcInputMsg *msg) 		{ (*rightHandler)(msg); }
		
		// ----------------------------------------- RDMA Support -----------------------------------------
		/// @entry Senders call this to setup an RDMA link to send the left matrix block
		inline void setupRDMALeft (RDMASetupRequestMsg<RDMApair_GSP_PC> *msg) 	
		{
			#ifndef PC_USE_RDMA
				CkPrintf("InputDataHandler[%d %d %d %d]. Someone called an RDMA method on me when RDMA has not been enabled.\n",
										thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z);
				CkAbort("InputDataHandler aborting because someone called an RDMA method when it has been turned off for me...\n");
			#endif
			#ifdef DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
				CkPrintf("InputHandler[%d,%d,%d,%d] Received an RDMA channel setup request for left data from %d. \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,msg->sender());
			#endif
			/// Get a copy of the transaction token sent in the message
			RDMApair_GSP_PC replyToken = msg->token();
			/// Verify that the setup request has come to the correct place
			CkAssert(replyToken.shouldSendLeft);
			/// Insert receiver identification into the handshake token
			replyToken.pcIndex = thisIndex;
			/// Delegate the actual RDMA setup work to the message handler
			leftHandler->setupRDMA(msg,&replyToken);
		}
		
		
		/// @entry Senders call this to setup an RDMA link to send the right matrix block
		inline void setupRDMARight(RDMASetupRequestMsg<RDMApair_GSP_PC> *msg)
		{
			#ifndef PC_USE_RDMA
				CkPrintf("InputDataHandler[%d %d %d %d]. Someone called an RDMA method on me when RDMA has not been enabled.\n",
										thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z);
				CkAbort("InputDataHandler aborting because someone called an RDMA method when it has been turned off for me...\n");
			#endif
			#ifdef DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
				CkPrintf("InputHandler[%d,%d,%d,%d] Received an RDMA channel setup request for right data from %d. \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,msg->sender());
			#endif
			/// Get a copy of the transaction token sent in the message
			RDMApair_GSP_PC replyToken = msg->token();
			/// Verify that the setup request has come to the correct place
			CkAssert(!replyToken.shouldSendLeft);
			/// Insert receiver identification into the handshake token
			replyToken.pcIndex = thisIndex;
			/// Delegate the actual RDMA setup work to the message handler
			rightHandler->setupRDMA(msg,&replyToken);
		}
		
	private:
		/// A message handler object that actually handles incoming msgs carrying left matrix data
		leftHandlerType  *leftHandler;
		/// A message handler object that actually handles incoming msgs carrying right matrix data
		rightHandlerType *rightHandler;		
};

//} // end namespace paircalc
//} // end namespace cp

/** This include is inside the #define block to avoid redefinition.
 * def.h files dont have include guards and can cause problems for modules with template chares/messages.
 */
#define CK_TEMPLATES_ONLY
	#include "inputDataHandler.def.h"
#undef CK_TEMPLATES_ONLY
#endif // INPUT_DATA_HANDLER_H

