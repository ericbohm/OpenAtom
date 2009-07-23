/** \file ckPairCalculator.h
 *
 */

#ifndef _ckPairCalculator_h_
#define _ckPairCalculator_h_
#undef  OLD_COMMLIB
#define OLD_COMMLIB 1
#define _PC_COMMLIB_MULTI_ 0
#include "debug_flags.h"
#include "pairutil.h"
#include "ckmulticast.h"
#include "ckhashtable.h"
#include "PipeBroadcastStrategy.h"
#include "BroadcastStrategy.h"
#include "DirectMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#ifdef CMK_BLUEGENEL
//#include "RectMulticastStrategy.h"
#endif
#include "MultiRingMulticast.h"
#include "NodeMulticast.h"


// Debugging flag for Verbose output
// #define _PAIRCALC_DEBUG_
// #define TEST_ALIGN

// If the machine is capable of RDMA...
#ifdef CMK_DIRECT
    // Enable GSpace-PairCalc RDMA
	#define PC_USE_RDMA
	#ifdef PC_USE_RDMA
		#include "cmidirect.h"
        // Turn RDMA on for the message data collator
		#define COLLATOR_ENABLE_RDMA
        // If RDMA debugging is toggled on, propagate the toggle
        #ifdef DEBUG_CP_PAIRCALC_RDMA
            #define DEBUG_MESSAGEDATACOLLATOR_RDMA
        #endif
	#endif
#endif


#ifdef CMK_BLUEGENEL

#define ALIGN16(x)        (int)((~15)&((x)+15))
#define BUNDLE_USER_EVENT
#define PC_FWD_DGEMM_SPLIT 1
#define PC_BWD_DGEMM_SPLIT 1
// to set split values, use the config parameters: gemmSplitFWk,
// gemmSplitFWm, etc ... 16 for happier align, factor of 6 good for BG/L?

#else

#define PC_FWD_DGEMM_SPLIT 0
#define PC_BWD_DGEMM_SPLIT 0

#endif

// flags to control semantic for matrix contents
#define NORMALPC   0  // standard
#define KEEPORTHO  1  // retain orthoT
#define PSIV       2  // multiply new psiV by retained orthoT

enum redtypes {section=0, machine=1, sparsecontiguous=2};
PUPbytes(redtypes);

#ifdef FORTRANUNDERSCORE
#define ZGEMM zgemm_
#define DGEMM dgemm_
#define DCOPY dcopy_
#define ZTODO ztodo_
#else
#define ZGEMM zgemm
#define DGEMM dgemm
#define DCOPY dcopy
#define ZTODO ztodo
#endif

extern ComlibInstanceHandle mcastInstanceCP;
#define _PAIRCALC_USE_DGEMM_

#ifdef _PAIRCALC_USE_BLAS_
extern "C" complex ZTODO( const int *N,  complex *X, const int *incX, complex *Y, const int *incY);
#endif

#ifdef _PAIRCALC_USE_DGEMM_
// extern "C" void DGEMM(char *,char *, int *,int *, int *,double *,double *,int *, double *,int *,double *,double *,int *);
extern "C" {void DGEMM (char *, char *, int *, int *, int *,double *,double *,
                        int *, double *, int *, double *, double *, int * );}
#endif

#ifdef _PAIRCALC_USE_ZGEMM_
extern "C" void ZGEMM(char *,char *, int *,int *, int *,complex *,complex *,int *,
                       const complex *,int *,complex *,complex *,int *);

extern "C" void DCOPY(int*,double *,int*, double *,int *);
#endif

typedef void (*FuncType) (complex a, complex b);
PUPmarshallBytes(FuncType);


#include "MessageDataCollator.h"
#include "paircalcMessages.h"
/// A shorter name for the namespace
namespace pc = cp::paircalc;
class paircalcInputMsg;
/// Typedefs for the collator types
typedef pc::MessageDataCollator<paircalcInputMsg,double> CollatorType;
struct RDMApair_GSP_PC;
#include "ckPairCalculator.decl.h"




///
class PairCalculator: public CBase_PairCalculator 
{
	public:
		/// @entry (obviously)
		PairCalculator(CProxy_InputDataHandler<CollatorType,CollatorType> inProxy, bool sym, int grainSize, int s, int blkSize, CkArrayID final_callbackid, int final_callback_ep, int callback_ep_tol, int conserveMemory, bool lbpaircalc, redtypes reduce, int orthoGrainSize, bool _AllTiles, bool streambw, bool delaybw, bool gSpaceSum, int gpriority, bool phantomSym, bool useBWBarrier, int _gemmSplitFWk, int _gemmSplitFWm, int _gemmSplitBW, bool expectOrthoT, int instance);
		/// Constructor for migration
		PairCalculator(CkMigrateMessage *);
		/// Destructor (nothing needs to be done?)
		~PairCalculator();
		/// Returns a pointer to the collator that will buffer the left matrix data (only for use by the corresponding InputDataHandler chare)
		inline CollatorType* leftHandler() const { return leftCollator; }
		/// Returns a pointer to the collator that will buffer the right matrix data (only for use by the corresponding InputDataHandler chare)
		inline CollatorType* rightHandler() const { return rightCollator; }
		/// @entry Method to send in the complete block of the left matrix
		void acceptLeftData(paircalcInputMsg *msg); 
		/// @entry Method to send in the complete block of the right matrix
		void acceptRightData(paircalcInputMsg *msg);
        /// NOT an entry method. Called locally from the acceptData* methods to launch the appropriate number-crunching method
        void launchComputations(paircalcInputMsg *aMsg);
		/// Forward path multiply driver. Prepares matrices, calls DGEMM, contributes results to Ortho subTiles and also passes relevant data to phantom PC chares
		void multiplyForward(bool flag_dp);
		/// @entry Simply redirects call to multiplyForward()
		void multiplyForwardRDMA() { multiplyForward(symmetric); }
		/// Piece up a tile and send all the pieces as this PC's contribution to the Ortho chares
		void contributeSubTiles(double *fullOutput);
		/// Contribute orthoGrainSized tiles of data (that are ready?) to the corresponding ortho chares
		void sendTiles(bool flag_dp);
		/// Receive data from ortho chares and copy into matrix
		void collectTile(bool doMatrix1, bool doMatrix2, bool doOrthoT, int orthoX, int orthoY, int orthoGrainSizeX, int orthoGrainSizeY, int numRecdBW, int matrixSize, double *matrix1, double* matrix2);
		/// @entry Initializes the section cookie and the reduction client. Called on startup as the chare world is being created
		void initGRed(initGRedMsg *msg);
		/// @entry During dynamics, each Ortho calls this on the Asymm loop PC instances to send its share of T back to avoid a race condition between Gamma and T.
		void acceptOrthoT(multiplyResultMsg *msg);
		/// @entry Backward path multiplication
		void multiplyResult(multiplyResultMsg *msg);
		/// Dynamics: PsiV Tolerance correction loop called on symmetric instances. Technically, backward path
		void multiplyPsiV();
		/// Multiplies Fpsi by T (from Ortho)
		void bwMultiplyDynOrthoT();
		/// @entry Simply forwards the call to multiplyResult(). @ntodo Dont seem to be any instances in source which call this method. Check.
		void multiplyResultI(multiplyResultMsg *msg);
		/// @entry a debugging tool: a barrier at the end of the backward path before anything is sent over to GSP
		void bwbarrier(CkReductionMsg *msg);
		/// multiplyPsiV() and multiplyResult() call this to perform the matrix multiply math on the backward path. This calls DGEMM routines for the same.
		void bwMultiplyHelper(int size, double *matrix1, double *matrix2, double *amatrix, double *amatrix2, bool unitcoef, int m_in, int n_in, int k_in, int BNAoffset, int BNCoffset, int BTAoffset, int BTCoffset, int orthoX, int orthoY, double beta, int ogx, int ogy);
		/// Called on the normal backward path (non-psiV) to set up the data sends to GSpace
		void bwSendHelper(int orthoX, int orthoY, int sizeX, int sizeY, int ogx, int ogy);
		/// @entry Send the results via multiple reductions as triggered by a prioritized message
		void sendBWResult(sendBWsignalMsg *msg);
		/// @entry
		void sendBWResultDirect(sendBWsignalMsg *msg);
		///
		void sendBWResultColumn(bool other, int startGrain, int endGrain);
		///
		void sendBWResultColumnDirect(bool other, int startGrain, int endGrain);
		/// @entry Initialize the section cookie for each slice of the result
		void initResultSection(initResultMsg *msg);
		/// PUP method
		void pup(PUP::er &);
		///
		void reorder(int *offsetMap, int *revOffsetMap, double *data, double *scratch);
		///
		void dumpMatrixDouble(const char *, double *,int,int, int xstart=0,int ystart=0, int xtra1=0, int xtra2=0);
		///
		void dumpMatrixComplex(const char *, complex *,int,int, int xstart=0, int ystart=0, int iter=0);
		///
		void dgemmSplitBwdM(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, double *B, double *bt, double *C);
		///
		void dgemmSplitFwdStreamMK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *C, int *ldc);
		///
		void dgemmSplitFwdStreamNK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *C, int *ldc);
		///
		void ResumeFromSync();
		/// @entry Something to sync?
		void lbsync() 
		{
			#ifdef _PAIRCALC_DEBUG_
				CkPrintf("[%d,%d,%d,%d] atsyncs\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z);
			#endif
			resumed=false;
			rck=0;
			AtSync();
		};

	private:

        /// Schedules the entry methods that send out the results to GSpace with appropriate priority
        void enqueueBWsend(bool unitcoef, int priority=1);
        /// Cleans up at end of an iteration (fw-bw computation loop); frees mem, resets counters etc
        void cleanupAfterBWPath(); 


		/// A handle to the co-located chare array that handles data input 
		CProxy_InputDataHandler<CollatorType,CollatorType> myMsgHandler;
		/// Data collators for the left and right matrix blocks
		CollatorType *leftCollator, *rightCollator;
		/// Flags indicating if the left and right matrix blocks have been received
		bool isLeftReady, isRightReady;
		int instance; 							/// Instance number of this run in path-integral beads
		int numRecd; 								///< number of messages received
		int numRecdBW; 							///< number of messages received BW
		int numRecdBWOT; 							///< number of messages received BW orthoT
		int numExpected; 							///< number of messages expected all
		int numExpectedX; 						///< number of messages expected x-axis
		int numExpectedY; 						///< number of messages expected y-axis
		int grainSize; 							///< number of states per chare
		int grainSizeX; 							///< number of states per chare x-axis
		int grainSizeY; 							///< number of states per chare y-axis
		int orthoGrainSize; 						///< number of states per ortho tile lower-bound
		int orthoGrainSizeRemX; 					///< sgrainSizeX % orthoGrainSize
		int orthoGrainSizeRemY; 					///< sgrainSizeY % orthoGrainSize
		int blkSize; 								///< number points in gspace plane
		int numStates; 							///< total number of states
		int numPoints; 							///< number of points in this chunk
		int numChunks; 							///< number of blocks the stateplane is divided into

		int streamCaughtR; 						///< number of rows caught so far R stream
		int streamCaughtL; 						///< number of rows caught so far L stream
		
		int numRecLeft; 							///< number of rows so far total left
		int numRecRight; 							///< number of rows so far total right
		
		int gemmSplitFWk; 						///< number of rows in split FW dgemm
		int gemmSplitFWm; 						///< number of columns in split FW dgemm
		int gemmSplitBW; 							///< number of rows in split BW dgemm
		
		
		int *LeftOffsets; 						///< index numbers of caught stream elements
		int *RightOffsets; 						///< index numbers of caught stream elements
		
		int *LeftRev; 							///< reverse index numbers of caught stream elements
		int *RightRev; 							///< reverse index numbers of caught stream elements
		
		double **outTiles; 						///< in output streaming we populate the tiles directly
		
		int *touchedTiles; 						///< tracker to detect when tiles are full
		
		bool symmetric; 							///< if true, one triangle is missing
		int conserveMemory; 						///< Tribool:: -1: dont conserve memory at all; 0: default; +1: conserve memory aggresively and free up matrices when not in use
		bool lbpaircalc;
		bool notOnDiagonal; 						///< being on or off diagonal changes many things
		bool symmetricOnDiagonal; 				///< diagonal symmetric special case
		
		bool phantomSym; 							///< phantoms exist to balance the BW path otherdata work
		
		bool expectOrthoT; 						///< Is true only in asymmetric, dynamics scenario. For PC instances in the asymmetric chare array, orthoT should arrive before end of fwd path 
		bool amPhantom; 							///< consolidate thisIndex.x<thisIndex.y && symmetric && phantomsym
		
		bool useBWBarrier;
		
		bool collectAllTiles; 					///< If true, don't stream compute on tiles in the backward path.
		
		redtypes cpreduce; 						///< which reducer we're using (defunct)
		CkArrayID cb_aid; 						///< bw path callback array ID
		int cb_ep; 								///< bw path callback entry point
		int cb_ep_tol; 							///< bw path callback entry point for psiV tolerance
		bool existsLeft; 							///< inDataLeft allocated
		bool existsRight; 						///< inDataRight allocated
		bool existsOut; 							///< outData allocated
		bool existsNew; 							///< newData allocated
		bool resumed; 							///< have resumed from load balancing
		
		bool PCstreamBWout; 						///< stream output from BW path
		bool PCdelayBWSend; 						///< use priority to delay BW output
		
		bool gSpaceSum; 							///< sum in gspace instead of reduction
		int gpriority; 							///< priority of msg to gspace
		
		complex *mynewData; 						///< results of bw multiply
		complex *othernewData; 					///< results of sym off diagonal multiply,
		                           //! or the C=-1 *inRight* orthoT +c in dynamics
		double *inDataLeft; 						///< the input pair to be transformed
		double *inDataRight; 						///< the input pair to be transformed
        paircalcInputMsg *msgLeft, *msgRight;   ///< Incoming messages with left and right matrix data that are kept around so that we can directly compute on them
		double *outData; 							///< results of fw multiply
		int actionType; 							///< matrix usage control [NORMAL, KEEPORTHO, PSIV]
		
		double *allCaughtLeft; 					///< unordered rows of FW input
		double *allCaughtRight; 					///< unordered rows of FW input
		
		
		double *inResult1; 						///< accumulate ortho or lambda
		double *inResult2; 						///< used in gamma calc (non minimization)
		
		/* to support the simpler section reduction*/
		int rck; 									///< count of received cookies
		CkGroupID mCastGrpIdOrtho; 				///< group id for multicast manager ortho
		int numOrthoCol; 							///< sGrainSizeX/orthoGrainSize
		int numOrthoRow; 							///< sGrainSizeY/orthoGrainSize
		int numOrtho; 							///< number of orthos in our grain = numOrthoCol*numOrthoRow
		
		CkGroupID mCastGrpId; 					///< group id for multicast manager bw
		
		CkSectionInfo *resultCookies; 			///< array of bw path section cookies
		CkSectionInfo *otherResultCookies; 		///< extra array of bw path section cookies for sym off diag, or dynamics
		
		CkCallback *orthoCB; 						///< forward path callbacks
		CkSectionInfo *orthoCookies; 				///< forward path reduction cookie
		int *columnCount; 						///< count of processed rows in BW by column
		int *columnCountOther; 					///< count of processed rows in BW by column
		
		/** Copy the results from outdata1 and outdata2 into the tiles.
		 * Iterate through the source array, look up the destination row in offsetsRow, destination 
		 * col in offsetsCol. This will be the destination row and column for the output if the
		 * output were considered as a single matrix. Use tileSize to map these values into the destination tile.
		 */
		void copyIntoTiles(double *source, double**dest, int sourceRows, int sourceCols, int *offsetsRow, int *offsetsCol, int *touched, int tileSize, int tilesPerRow );
};

//forward declaration
CkReductionMsg *sumMatrixDouble(int nMsg, CkReductionMsg **msgs);
CkReductionMsg *sumBlockGrain(int nMsg, CkReductionMsg **msgs);

void manmult(int numrowsA, int numRowsB, int rowLength, double *A, double *B, double *C, double alpha);
#endif
