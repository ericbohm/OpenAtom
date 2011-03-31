#ifndef CK_PAIR_CALCULATOR_H
#define CK_PAIR_CALCULATOR_H

#undef  OLD_COMMLIB
//#define OLD_COMMLIB 1
#include "debug_flags.h"

#include "ckmulticast.h"
#include "ckhashtable.h"

#include "PipeBroadcastStrategy.h"
#include "BroadcastStrategy.h"
#include "DirectMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#include "MultiRingMulticast.h"
#include "NodeMulticast.h"


// If the machine is capable of RDMA...
#if CMK_DIRECT
    // Enable GSpace-PairCalc RDMA
	//#define PC_USE_RDMA
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


/* If the machine is a BGL, split the GEMMs to intersperse them with CmiNetworkProgress calls
 *
 * To set split values, use the config parameters: gemmSplitFWk, gemmSplitFWm, etc ... 
 * 16 for happier align, factor of 6 good for BG/L?
 */
#ifdef CMK_BLUEGENEL
    #define ALIGN16(x)        (int)((~15)&((x)+15))
    #define BUNDLE_USER_EVENT
    #define PC_FWD_DGEMM_SPLIT 1
    #define PC_BWD_DGEMM_SPLIT 1
#else
    #define PC_FWD_DGEMM_SPLIT 0
    #define PC_BWD_DGEMM_SPLIT 0
#endif


// Define the GEMM macros that paircalc will use to invoke the appropriate matrix multiplys
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


// Flags to control semantic for matrix contents along different paths thro the PairCalcs
#define NORMALPC   0  ///< standard
#define KEEPORTHO  1  ///< retain orthoT
#define PSIV       2  ///< multiply new psiV by retained orthoT

/// The type of data that is input from GSpace and operated on by the PairCalcs
typedef double inputType;

extern ComlibInstanceHandle mcastInstanceCP;

extern "C" {void DGEMM (char *, char *, int *, int *, int *,double *,double *, int *, double *, int *, double *, double *, int * );}

#include "MessageDataCollator.h"
#include "pcMessages.h"
#include "pcConfig.h"
#include "ckPairCalculator.decl.h"




/**
 * PairCalculator (PC) is a 4D chare array that is, at its heart, a glorified wrapper
 * for a bunch of matrix multiplications. It serves the function of managing the
 * complexity of the decomposition and communication logic that is required to
 * scale well. In the process, it introduces some more :).
 *
 * There are two separate PairCalc chare arrays instantiated per OpenAtom simulation
 * instance. Both of these serve GSpace and process data whenever invoked by GSpace.
 * Both these PC arrays share a single Ortho chare array to which some of the multiply
 * work is delegated.
 *
 * The extra complications are for parallelization and the multiplication
 * of the forces and energies.
 *
 * In normal use the calculator is created.  Then forces are sent to it
 * and multiplied in a big dgemm.  Then this result is reduced to the
 * answer matrix and shipped back.  The received left and/or right data is
 * retained for the backward pass which is triggered by the finishPairCalc
 * call.  This carries in another set of matrices for multiplication.
 * The results are again reduced and cast back.  Thus terminating the life
 * cycle of the data in the pair calculator.  As the calculator will be
 * reused again throughout each iteration the calculators themselves are
 * only created once.
 *
 *
 * The paircalculator is a 4 dimensional array.  Those dimensions are:
 *            w: gspace state plane (the second index of the 2D gspace)
 *            x: coordinate offset within plane (a factor of grainsize)
 *            y: coordinate offset within plane (a factor of grainsize)
 *            z: chunk offset within array of nonzero points
 *       So, for an example grainsize of 64 for a 128x128 problem:
 *        numStates/grainsize gives us a 2x2 decomposition.
 *        1st quadrant ranges from [0,0]   to [63,63]    index [w,0,0,0]
 *        2nd quadrant ranges from [0,64]  to [63,127]   index [w,0,64,0]
 *        3rd quadrant ranges from [64,0]  to [127,63]   index [w,64,0,0]
 *        4th quadrant ranges from [64,64] to [127,127]  index [w,64,64,0]
 *
 *       0   64   127
 *     0 _________
 *       |   |   |
 *       | 1 | 2 |
 *    64 ---------
 *       |   |   |
 *       | 3 | 4 |
 *   127 ---------
 *
 *
 *
 * Further complication arises from the fact that each plane is a
 * cross-section of a sphere.  So the actual data is sparse and is
 * represented by a contiguous set of the nonzero elements.  This is
 * commonly referred to as numPoints or size within the calculator.
 *
 * In the dynamics case there are two additional complications. In the
 * backward path of the asymmetric pairCalculator we will receive 2 input
 * matrices, gamma and orthoT.  Where orthoT came from the
 * orthonormalization following the symmetric pairCalculator.  And
 * gamma was produced by a multiplication with orthoT in the Ortho
 * object.
 *
 * If ortho falls out of tolerance then Ortho will signal the GSP that
 * a tolerance update is needed.  We then proceed with the psi
 * calculation as normal.  On receipt of newpsi, Gspace will then
 * react by sending the PC the Psi velocities (PsiV) in the same way
 * (acceptLeft/RightData) that it sends Psi, but with the psiv flag set
 * true.  These will be recorded in the left data array.  We will then
 * multiply with the orthoT we kept from the previous invocation (psi)
 * of the backward path.  We then ship the corrected velocities back
 * to gspace via the acceptnewVpsi reduction.  Same procedure as for
 * acceptnewpsi, just a different entry method.
 *
 * Fourth dimension decomposition is along the axis of the nonzero
 * values in gspace.  Therefore it is fundamentally different from the
 * 2nd and 3rd dimensions which divide the states up into
 * (states/grainsize)^2 pieces.  The fourth dimension divides along
 * the nonzeros of gspace.  A X,0,0,N division will have the entirety
 * of a state, but only a K/Nth (where K is the number of nonzero
 * elements) chunk of the nonzero values.  It can therefore perform
 * the dgemm on that chunk, its multicast degree will be 1, and have a
 * portion of the total solution.  Thereby reducing the PC inbound
 * communication volume substantially.  This comes at the cost of an
 * additional reduction.  The result for the nonzero chunks has to be
 * pasted together to form the result for the entire nonzero.  Then
 * the results are summed together across the planes to produce the
 * complete S or L matrix.  Only the first of those reductions is new.
 *
 * More about this "extra" reduction.  If we consider the Multiply as
 * C = AB.  Where A is nstates x numpoints and B is numpoints x
 * nstates to result in C of nstates x nstates.  The 4th dim
 * decomposition chops only the inner index numpoints.  Thereby
 * resulting in numblocks C submatrices all of size nstates x nstates.
 * Making C from numblock C(i) is just matrix addition.  So for the
 * forward path there is in fact no "extra" reduction or stitching
 * necessary for the 4th dim decomposition.  All the "stitching" is in
 * the statewise decompostion for the forward path.  So the only
 * change for the forward path is in adding the numblock elements to
 * the reduction tree.
 *
 *
 * An important distinction between these methods is that in the
 * absence of a grainsize decomposition the sections for the second
 * reduction to ortho are essentially arbitrary with respect to the
 * PairCalculator decomposition.
 *
 *
 * Similarly, the backward path for a chunk with grainsize==nstates
 * needs input matrices of size nstates X nstates.  Which means that
 * the backward path cannot proceed until all ortho sections broadcast
 * their pieces of the input matrices.  The backward path reduction to
 * gspace becomes richer in that now each chunk is contributing at an
 * offset.  So the acceptnew[psi|lambda|vpsi] methods would all need
 * to paste the input of a contribution at the correct offset.  This
 * recalls the old contiguousreducer logic which did that pasting
 * together in its reduction client.  A functionality that is arguably
 * worth resurrection.  Just in a form that can actually be
 * comprehended by people not from planet brainiac.
 *
 * Which means we should add a distinct parameter for the number of ortho
 * objects.  We'll also need to come up with a way for it to map its
 * grainsize sections onto the chunketized PC chare array.  The
 * constraints on this mapping are only that they should use as many
 * PCs as we can.  The PCs will use the section reduction in their
 * forward path reduction to deposit the S (or lambda) matrix in ortho.
 * Ortho will have to broadcast its T (or lambda, or orthoT and gamma)
 * to the PairCalculator.
 *
 * Which returns us to the bad old days of broadcasting matrices
 * around.  This can be ameliorated slightly by using [nokeep]
 * messages so that we reduce the number of copies to 1 per PE.  But
 * you still have to send numPE*nstates*nstates doubles around (times
 * 2 in the dynamics symmetric case to carry ortho and gamma).
 * Luckily broadcasts turn out to be pretty fast on BG/L.  So this may
 * not be so bad.  The tradeoff against the much larger nonzero
 * multicast is net positive due to the larger constant size of the
 * nonzeros compared to nstates.
 *
 * These communication patterns become more complex in the hybrid case
 * where we have both grainsize and chunksize.  The ortho->PC section
 * mapping could revert to using grainsize, but now has to sum across
 * all the chunks in that grain.
 *
 * If we want independant control over the number of ortho objects
 * then we need to support the overlap issues where grainsize objects
 * do not map nicely onto ortho objects.
 *
 * The reduction out of the backward path of the paircalculator is the
 * one which is made more complicated by 4th dimension decomposition.
 * Here we are sending the transformed Psi, which is necessarily
 * numpoints in size.  Therefore the reduction requires stitching of
 * blocks while summing within a block to assemble the entire g-chare
 * matrix of points.  In practice this is just one big reduction with
 * the userfield used to indicate the chunk number so gspace can lay
 * the results out where they belong.
 *
 * In its current form orthoGrainSize must be an even mod of
 * sGrainSize.  This restriction avoids overlap and boundary issues
 * for tiles between ortho and the calculator.  Thereby avoiding some
 * rather gnarly issues in the setup of multicasts and reductions.
 * This is a fairly minor restriction as long as we do not require
 * nstates % sgrainsize==0 or nstates & orthograinsize.
 *
 * sGrainSize need not be even mod of the number of states.  nstates %
 * sGrainSize = remainder requires some careful handling in the
 * code. Whenever this occurs the multiplies and communications for
 * the border edges which carry the remainder have to cope with
 * asymmetric multiplies and funky remainder logic in communicating
 * the result.
 *
 *
 *
 * Data from GSP travels to the PairCalculators via an InputDataHandler chare 
 * array of the same dimensions as, and bound to, the PairCalculator array. 
 * Appropriate sections/lists of this input handler array for multicasting 
 * the data from GSP to are built in makeLeftTree() and makeRightTree().
 *
 * Each iteration of the GSP-PC-Ortho-PC-GSP loop is started by GSP calling
 * startPairCalcLeft() and startPairCalcRight(). These are simply #defines
 * that turn into the appropriate function: sendLeftData() and sendRightData()
 * or their RDMA equivalents. The input handler chares then collate all the
 * incoming data and then wake their corresponding PC chares once all the data
 * is available. The PCs then do their thing and the result is returned via 
 * the callback set in the create routine (which happens to be Ortho).
 *
 * The backward path is triggered by:
 * finishPairCalcSection(PairCalcID, datasize, data *)
 * Its result is returned via the callback entry point passed in
 * during creation
 *
 * The results of the backward path are returned in a set of section
 * reductions.  The reduction contributes its slice of its matrix of
 * doubles with the offset=thisIndex.z.  The client then returns the
 * sum of each slice to the GSP element that created the section with
 * the offset so it can be copied into the correct place in the points
 * array.
 *
 *
 * NOTE: The magic number 2 appears in 2 contexts.  Either we have
 * twice as many inputs in the non diagonal elements of the matrix.
 * Or in the case of transforming our arrays of complex into arrays of
 * doubles.  This transformation is done so we can use DGEMM instead
 * of ZGEMM.  Motivated by the empirical discovery that BLAS
 * implementors do all their work on DGEMM and leave ZGEMM out in the
 * unoptimized cold.  The latter issue crops up everywhere we have to
 * do allocation or manipulations of the input.  Could arguably be
 * abstracted away by absorbing it into a size field.
 *
 *
 * After the main::main proc 0 phase the array sections used to populate
 * and return data from the paircalculator are called.
 * The forward path section reduction (to and from ortho) is initialized
 * via initOneRedSect() the backward path is initialized via the
 * appropriate makeOneResultSection_X() call.  In each case the call
 * should be made by each GSP or Ortho object.  That way each one has its
 * own proxy and the section tree will only include relevant processors.
 *
 * The chunk decomposition changes the setup substantially.  In chunk
 * decomposition we send a piece of the nonzero points of gspace to
 * the paircalculators.  It is not a multicast.  Each GSP[P,S] will
 * send its ith chunk to PC[P,0,0,i].  Once nstate chunks arrive at a
 * PC the multiply can proceed.
 *
 * In the hybrid case this becomes a multicast of chunks to the
 * appropriate state decomposition destination as before.
 *
 */
class PairCalculator: public CBase_PairCalculator 
{
	public:
		/// @entry (obviously)
		PairCalculator(CProxy_InputDataHandler<CollatorType,CollatorType> inProxy, const pc::pcConfig _cfg);
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
		void multiplyForwardRDMA() { multiplyForward(cfg.isSymmetric); }
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


        /// A private copy of the input configurations
        pc::pcConfig cfg;
		/// A handle to the co-located chare array that handles data input 
		CProxy_InputDataHandler<CollatorType,CollatorType> myMsgHandler;
		/// Data collators for the left and right matrix blocks
		CollatorType *leftCollator, *rightCollator;
		/// Flags indicating if the left and right matrix blocks have been received
		bool isLeftReady, isRightReady;
		int numOrthoCookiesRecvd;                                       ///< counts the number of ortho cookies recieved at startup
		int numRecd; 								///< number of messages received
		int numRecdBW; 							///< number of messages received BW
		int numRecdBWOT; 							///< number of messages received BW orthoT
		int numExpected; 							///< number of messages expected all
		int numExpectedX; 						///< number of messages expected x-axis
		int numExpectedY; 						///< number of messages expected y-axis
		int grainSizeX; 							///< number of states per chare x-axis
		int grainSizeY; 							///< number of states per chare y-axis
		int orthoGrainSizeRemX; 					///< sgrainSizeX % orthoGrainSize
		int orthoGrainSizeRemY; 					///< sgrainSizeY % orthoGrainSize
		int blkSize; 								///< number points in gspace plane
		int numPoints; 							///< number of points in this chunk

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
		
		bool notOnDiagonal; 						///< being on or off diagonal changes many things
		bool symmetricOnDiagonal; 				///< diagonal symmetric special case
		
		bool expectOrthoT; 						///< Is true only in asymmetric, dynamics scenario. For PC instances in the asymmetric chare array, orthoT should arrive before end of fwd path 
		bool amPhantom; 							///< consolidate thisIndex.x<thisIndex.y && cfg.isSymmetric && phantomsym
		
		CkArrayID cb_aid; 						///< bw path callback array ID
		int cb_ep; 								///< bw path callback entry point
		int cb_ep_tol; 							///< bw path callback entry point for psiV tolerance
		bool existsLeft; 							///< inDataLeft allocated
		bool existsRight; 						///< inDataRight allocated
		bool existsOut; 							///< outData allocated
		bool existsNew; 							///< newData allocated
		bool resumed; 							///< have resumed from load balancing
		
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
#endif // CK_PAIR_CALCULATOR_H

