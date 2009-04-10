#include "RDMAMessages.decl.h"
// Include cmidirect only if RDMA is enabled
#ifdef ENABLE_RDMA_HANDSHAKES
	#include "cmidirect.h"
#endif
#include <ostream>

#ifndef RDMA_MESSAGES_H
#define RDMA_MESSAGES_H

/// Based on whether RDMA is enabled, the handle type is either the actual handle or just an empty structure
#ifdef ENABLE_RDMA_HANDSHAKES
	typedef struct infiDirectUserHandle rdmaHandleType;
	//PUPmarshallBytes(struct infiDirectUserHandle); ///@todo: Determine if something needs to be done to enable PUPing
#else
	typedef struct {} rdmaHandleType;
#endif


/// A request from a data sender to setup an RDMA link. Initiates the sender-receiver handshake required to setup such a link 
template <class tokenType>
class RDMASetupRequestMsg: public CMessage_RDMASetupRequestMsg<tokenType>
{
	public:
		/// Get a copy of the handshake token the sender created
		inline tokenType token() const 						{ return handshakeToken; }
		/// Get an integer representation of the sender's ID
		inline int sender() const 							{ return senderID; }
		/// Get the PE number of the sender
		inline int senderPE() const 						{ return senderProc; }
		/// Get the number of records of data that will be sent by RDMA 
		inline int numDataUnits() const 					{ return (numUnits); }
		/// Get hold of the callback the sender has configured
		inline const CkCallback& senderCallback() const 	{ return callbk; }
		
		/// Constructors
		RDMASetupRequestMsg(): 
				senderID(-1), senderProc(-1), numUnits(0) 	{}
		RDMASetupRequestMsg(const tokenType tkn, const int sndr, const int sndrPE, const int dataSize,const CkCallback cb):
				handshakeToken(tkn),
				senderID(sndr), senderProc(sndrPE),
				numUnits(dataSize), callbk(cb) 				{}
		
	private:
		int senderID, senderProc, numUnits;
		CkCallback callbk;
		tokenType handshakeToken;

    /// Stream inserter. Note: not a member function
    inline friend std::ostream& operator<< (std::ostream &out, const RDMASetupRequestMsg<tokenType> &obj)
    {
        out<<"RDMA setup request from "<<obj.senderID<<" on proc "<<obj.senderProc
        <<" sending "<<obj.numUnits<<" units. Token data: "<<obj.handshakeToken;
        return out;
    }
};




/// Reply from data receiver to the data sender indicating completion of setup on the receiver side
template <class tokenType> 
class RDMASetupConfirmationMsg: public CMessage_RDMASetupConfirmationMsg<tokenType>
{
	public:
		RDMASetupConfirmationMsg() 							{}
		RDMASetupConfirmationMsg(const tokenType tkn, const rdmaHandleType hndl): handshakeToken(tkn), ourHandle(hndl) {}
		inline tokenType token() const 						{ return handshakeToken; }
		inline rdmaHandleType handle() const 				{ return ourHandle; }
		
	private:
		tokenType handshakeToken;
		rdmaHandleType ourHandle;

    /// Stream inserter. Note: not a member function
    inline friend std::ostream& operator<< (std::ostream &out, const RDMASetupConfirmationMsg<tokenType> &obj)
    {
        out<<"RDMA setup confirmation msg. Token data: "<<obj.handshakeToken;
        return out;
    }
};




/** A (hopefully) tiny token that is unique to every data sender-receiver pair, 
 * and is shared by them during the RDMA setup process. This simply encapsulates 
 * information about the sender and the receiver that each other need.
 * 
 * This can be defined for every unique RDMA pair, and the RDMA setup messages 
 * can be configured with these so that the messages stay somewhat general. This
 * structure holds information relevant to GSpace and PairCalc communication.
 */
class RDMApair_GSP_PC
{
	public:
		RDMApair_GSP_PC(): shouldSendLeft(true), symmetric(true)
        {
            gspIndex.x = gspIndex.y = 0;
            pcIndex.w = pcIndex.x = pcIndex.y = pcIndex.z = 0;
        }
		bool shouldSendLeft, symmetric;
		CkIndex2D gspIndex;
		CkIndex4D pcIndex;

    /// Stream inserter. Note: not a member function 
    inline friend std::ostream& operator<< (std::ostream &out, const RDMApair_GSP_PC &obj)
    {
        out<<"GSpace["<<obj.gspIndex.x<<","<<obj.gspIndex.y
        <<"] - PC["<<obj.pcIndex.w<<","<<obj.pcIndex.x<<","<<obj.pcIndex.y<<","<<obj.pcIndex.z<<","<<obj.symmetric<<"] "
        <<(obj.shouldSendLeft?"Left ":"Right ");
        return out;
    }
};



/** This include is inside the #define block to avoid redefinition.
 * def.h files dont have include guards and can cause problems for modules with template chares/messages.
 */
#define CK_TEMPLATES_ONLY
	#include "RDMAMessages.def.h"
#undef CK_TEMPLATES_ONLY
#endif // RDMA_MESSAGES_H
