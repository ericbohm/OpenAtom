#include "charm++.h"          ///< This code uses CmiMemcpy and CkCopyMessage, which can be accessed via this include. 
 
#ifndef MESSAGE_DATA_COLLATOR_H
#define MESSAGE_DATA_COLLATOR_H

namespace cp {
namespace paircalc {


/** Class that does the task of buffering data sent in via messages for an expected number of messages before
 * calling on a handler to process them. Assumes that all messages carry the same amount of data.
 *
 * This class is not generic enough to warrant being a template, but it might still be useful to start it as one, 
 * so that it has a chance to becoming abstracted enough with time. 
 *
 * @note: Class msgType should provide the following methods:
 * - int numDataUnits() const  : returns the number of units of dataType in the message (or message size) 
 * - int sender() const        : returns the sender ID or some representation of the senderID as an int
 * - dataType* data()          : returns a pointer to the data held by the message. Does not strictly have to be 
 *                               of dataType as it is only used by CmiMemcpy. 
 *
 * @note: Class triggerType should be a functor that takes three arguments:
 * - void operator() (const dataType *collatedData, const unsigned int numMessages, const unsigned int numUnitsInMsg)
 */
template <class msgType, typename dataType, class triggerType>
class MessageDataCollator 
{
	public:
		/// @warning: Do not use this. This exists only to appease charm migration code which screams for default constructors.
		MessageDataCollator() {}
		/** Collators should be constructed thus. Copying has not been restricted so that collator classes can be passed in as
		 * arguments. Hence, be careful when copying collators to make sure multiple collators dont use the same memory region as
		 * a data buffer and kill you.
		 * 
		 * @warning: This collator class should NEVER delete ptr, it can at best allocate memory if ptr is NULL.
		 */
		MessageDataCollator(triggerType trigger, const unsigned int nExpected, bool manageBuffer=false, const unsigned int offsetStart=0):
			doneTrigger(trigger),
			isBufferMine(manageBuffer),
			dataBuffer(0),
			bufferSizeInUnits(0),
			numBatchesProcessed(0),
			numExpected(nExpected),
			numReceived(0),
			minSender(offsetStart),
			sampleMsg(0)    
		{
			#ifdef DEBUG_CP_PAIRCALC_MESSAGEDATACOLLATOR
				CkPrintf("\nMessageDataCollator instantiated: numExpected = %d and buffer management = %d on(1)/off(0)",numExpected,isBufferMine);
			#endif
		}
		
		/// Destructor
		~MessageDataCollator()
		{
			if (isBufferMine)
				delete [] dataBuffer;
			delete sampleMsg;
		}
		/// Makes this a functor so that it can be fed messages by any entity that just knows the MessageHandler API. 
		void operator() (msgType *msg); 
		/// Get the number of messages this collator has been configured to buffer
		inline unsigned int getNumExpected() const { return numExpected; }
		/// Check the number of messages that have arrived
		inline unsigned int getNumReceived() const { return numReceived; }
		/// Get the number of batches of incoming messages already processed
		inline unsigned int getNumBatchesDone() const { return numBatchesProcessed; }
		/** Get hold of a copy of one of the incoming messages that is then managed by the user
		 *  
		 * Useful in case the msg has info other than the data() field. These can then 
		 * be extracted in a case-specific manner, by the client code directly. 
		 * 
		 * @note: This class makes no guarantees on which message will be held as a sample msg.
		 */
		inline msgType* getSampleMsg() { return reinterpret_cast<msgType*> ( CkCopyMsg(reinterpret_cast<void**>(&sampleMsg)) ); }
		
	private:
		/// This is the functor that is passed in to setup the post-collation operation. Typically a callback
		triggerType doneTrigger;
		/// If buffer is mine, I manage it, viz. reuse the same buffer for every batch of msgs where possible. If not, I allocate a fresh buffer and let the client manage it.
		bool isBufferMine;
		/// A memory region where incoming message data is collated
		dataType* dataBuffer;
		/// Handle to a sample message copied from every batch of incoming messages
		msgType *sampleMsg;
		/// Counters to keep track of the messages
		unsigned int numExpected, numReceived, numBatchesProcessed, bufferSizeInUnits;
		/** Senders should have an ID accessed by msgType::sender(). This stores the minimum value of sender(), 
		 * so that offsets into the memory buffer can be computed
		 *  
		 * @todo: A better way to accomplish that leaves the collator class less rigid might be better
		 */
		 unsigned int minSender;
};



template <class msgType,typename dataType, class triggerType>
void MessageDataCollator<msgType,dataType,triggerType>::operator() (msgType *msg)
{
	/** @todo: Debate if some kind of check is needed to ensure that the message received is the correct iteration etc.
	 * For e.g, if the message defined a uniqID() method that would return a unique tag indicating the group the message
	 * belonged to, then all messages that belonged to this group could be buffered together. uniqID() could return
	 * iteration number if messages from a single iteration are to be grouped, or the sender ID if messages from a single
	 * sender are to be grouped together. Other possible scenarios where it is possible to receive messages from separate
	 * collation groups might occur.
	 */
	/// Get the message size
	int numUnits = msg->numDataUnits();
	/// Get the required buffer size
	int reqdBufferSizeInUnits = numExpected * numUnits;
	/// If buffer exists AND I manage it AND its not big enough for this batch of messages, delete it
	if (dataBuffer && isBufferMine && reqdBufferSizeInUnits > bufferSizeInUnits)
	{
		delete [] dataBuffer;
		dataBuffer = 0;
	}
	/// If the buffer is unallocated ... 
	if (!dataBuffer)
	{
		/// Get some memory to hold the message data
		bufferSizeInUnits =  reqdBufferSizeInUnits;
		dataBuffer = new dataType[bufferSizeInUnits];
		/// Zero out the memory region
		bzero(dataBuffer,bufferSizeInUnits * sizeof(dataType));
	} /// endif
	/// Compute the offset for this message. This will depend on the sender index and thisIndex. 
	int offset = msg->sender() - minSender;
	/// Copy the data from the message into the contiguous data buffer
	CmiMemcpy(&(dataBuffer[offset*numUnits]), msg->data(), numUnits *sizeof(dataType));
	/// Increment the tally of received messages and if the expected number have arrived ...
	if (++numReceived >= numExpected)
	{
		/// Increment the number of message batches that have been processed
		numBatchesProcessed++;
		/// Create a copy of this message for the sample
		if (sampleMsg)
			delete sampleMsg;
		sampleMsg = reinterpret_cast<msgType*> ( CkCopyMsg(reinterpret_cast<void**>(&msg)) );
		/// If the expected number of messages have arrived, call the trigger functor
		doneTrigger(dataBuffer,numReceived,numUnits);
		/// Reset the number received so that we can start collating all over again
		numReceived = 0;
		/// If you are not managing the buffer, then you need to forget about it
		if (!isBufferMine)
			dataBuffer = 0;
	} /// endif
}

} // end namespace paircalc
} // end namespace cp

#endif // MESSAGE_DATA_COLLATOR_H
