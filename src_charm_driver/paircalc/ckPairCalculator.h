/** \file ckPairCalculator.h
 *
 */

#ifndef _ckPairCalculator_h_
#define _ckPairCalculator_h_
#define _PC_COMMLIB_MULTI_ 0
#include "pairutil.h"
#include "ckmulticast.h"
#include "ckhashtable.h"
#include "PipeBroadcastStrategy.h"
#include "BroadcastStrategy.h"
#include "DirectMulticastStrategy.h"
#include "RingMulticastStrategy.h"
#ifdef CMK_VERSION_BLUEGENE
//#include "RectMulticastStrategy.h"
#endif
#include "MultiRingMulticast.h"
#include "NodeMulticast.h"
#include "../../include/debug_flags.h"

// Debugging flag for Verbose output
// #define _PAIRCALC_DEBUG_
// #define TEST_ALIGN

#ifdef CMK_DIRECT
#define PC_USE_RDMA 1
#ifdef PC_USE_RDMA
#include "cmidirect.h"
#endif
#endif


#ifdef CMK_VERSION_BLUEGENE

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

#include "ckPairCalculator.decl.h"

class RDMAHandle{
 public:
#ifdef PC_USE_RDMA
  struct infiDirectUserHandle handle;
#endif

  ~RDMAHandle(){}
#ifdef PC_USE_RDMA
  RDMAHandle(struct infiDirectUserHandle _handle): handle(_handle){}
  RDMAHandle(){handle.handle=-1;}
#else
  // just to keep the compiler happy
  RDMAHandle(){}
#endif
};

#ifdef PC_USE_RDMA
PUPmarshallBytes(struct infiDirectUserHandle);
#endif

class RDMAHandleMsg : public CMessage_RDMAHandleMsg
{
 public:
#ifdef PC_USE_RDMA
  struct infiDirectUserHandle rhandle;
#endif
  int index;
  int grain;
  bool left;
  bool symmetric;
#ifdef PC_USE_RDMA
  void init(struct infiDirectUserHandle _rhandle, int _index, int _grain,bool _left, bool _symmetric)
    {
      rhandle=_rhandle;
      index=_index;
      grain=_grain;
      left=_left;
      symmetric=_symmetric;
    }
#else
  // who cares?
#endif
  friend class CMessage_RDMAHandleMsg;
};  

class initGRedMsg : public CkMcastBaseMsg, public CMessage_initGRedMsg {
 public:
  CkCallback cb;
  CkGroupID mCastGrpId;
  CkCallback synccb;
  bool lbsync;
  int orthoX;
  int orthoY;
  friend class CMessage_initGRedMsg;
};

class initResultMsg : public CkMcastBaseMsg, public CMessage_initResultMsg {
 public:
  int offset;
  int dest;
  CkGroupID mCastGrpId;
  CkCallback synccb;
  bool lbsync;
  friend class CMessage_initGRedMsg;
};

class sendBWsignalMsg : public CMessage_sendBWsignalMsg{
 public:
  bool otherdata;
};

class sendFWRDMAsignalMsg : public CMessage_sendBWsignalMsg{
 public:
  bool flag_dp;
};

class mySendMsg : public CMessage_mySendMsg {
 public:
  int N;
  complex *data;
  friend class CMessage_mySendMsg;
};

class partialResultMsg : public CMessage_partialResultMsg {
 public:
  complex *result;
  int N;
  int myoffset;
  void init(int _size, int _myoffset, complex *_points)
    {
      N=_size;
      myoffset=_myoffset;
      CmiMemcpy(result,_points,N*sizeof(complex));
    }

  friend class CMessage_partialResultMsg;

#ifdef CMK_VERSION_BLUEGENE
static void* alloc(int msgnum, size_t sz, int *sizes, int pb) {
  int offsets[2];
  offsets[0] = ALIGN16(sz);
  if(sizes==0)
    offsets[1] = offsets[0];
  else
    offsets[1] = offsets[0] + ALIGN16(sizeof(complex)*sizes[0]);
  partialResultMsg *newmsg = (partialResultMsg *) CkAllocMsg(msgnum, offsets[1], pb);
  newmsg->result = (complex *) ((char *)newmsg + offsets[0]);
  return (void *) newmsg;
}

#endif  
};

class priorSumMsg : public CMessage_priorSumMsg {
 public:
  complex *result;
  int N;
  int priority;
  CkCallback cb;

  friend class CMessage_priorSumMsg;

};

class calculatePairsMsg : public CkMcastBaseMsg, public CMessage_calculatePairsMsg {
 public:
  int size;
  int sender;
  complex *points;
  bool fromRow;
  bool flag_dp;
  bool doPsiV;
  int blkSize;
  void init(int _size, int _sender, bool _fromRow, bool _flag_dp, complex *_points , bool _doPsiV, int _blkSize)
    {
      size=_size;
      sender=_sender;
      fromRow=_fromRow;
      flag_dp=_flag_dp;
      doPsiV=_doPsiV;
      blkSize=_blkSize;

      CmiMemcpy(points,_points,size*sizeof(complex));
    }
  friend class CMessage_calculatePairsMsg;

};

class phantomMsg : public CMessage_phantomMsg {
 public:
  int size;
  int numPoints;
  double *points;
  int blkSize;
  int actionType;
  void init(int _size, int _numPoints, bool _flag_dp, double *_points, int _blkSize, int _actionType)
    {
      size=_size;
      numPoints=_numPoints;
      blkSize=_blkSize;
      CmiMemcpy(points,_points,size*sizeof(double));
      actionType=_actionType;
    }

};

class multiplyResultMsg : public CkMcastBaseMsg, public CMessage_multiplyResultMsg {
 public:
  double *matrix1;
  double *matrix2;
  int size;
  int size2;
  int orthoX;
  int orthoY;
  int actionType;
  void init(int _size, int _size2, double *_points1, double *_points2, int _orthoX, int _orthoY, bool _actionType)
    {
      size=_size;
      size2=_size2;
      orthoX=_orthoX;
      orthoY=_orthoY;
      CmiMemcpy(matrix1,_points1,size*sizeof(double));
      CmiMemcpy(matrix2,_points2,size2*sizeof(double));
      actionType=_actionType;
    }
  void init1(int _size, double *_points1, int _orthoX, int _orthoY,int _actionType)
    {
      size=_size;
      size2=0;
      orthoX=_orthoX;
      orthoY=_orthoY;
      CmiMemcpy(matrix1,_points1,size*sizeof(double));
      actionType=_actionType;
      // this field does nothing in minimization
      matrix2=NULL;
    }
#ifdef CMK_VERSION_BLUEGENE
  // if we use our own allocator we can get 16 byte alignment
  // to please BGL
 static  void *alloc(int msgnum, size_t sz, int *sizes, int pb) {
    int offsets[3];
    offsets[0] = ALIGN16(sz);
    if(sizes==0)
      offsets[1] = offsets[0];
    else
      offsets[1] = offsets[0] + ALIGN16(sizeof(double)*sizes[0]);
    if(sizes==0)
      offsets[2] = offsets[0];
    else
      offsets[2] = offsets[1] + ALIGN16(sizeof(double)*sizes[1]);
    multiplyResultMsg *newmsg = (multiplyResultMsg *) CkAllocMsg(msgnum, offsets[2], pb);
    newmsg->matrix1 = (double *) ((char *)newmsg + offsets[0]);
    newmsg->matrix2 = (double *) ((char *)newmsg + offsets[1]);
    return (void *) newmsg;
  }

#endif
  friend class CMessage_multiplyResultMsg;
};

class entireResultMsg : public CMessage_entireResultMsg {
 public:
  double *matrix;
  int size;
  bool symmetric;
  void init(int _size, double *_points, bool _symmetric)
    {
      size=_size;
      symmetric=_symmetric;
      CmiMemcpy(matrix,_points,size*sizeof(double));
    }
  friend class CMessage_entireResultMsg;
};

class entireResultMsg2 : public CMessage_entireResultMsg2 {
 public:
  double *matrix1;
  double *matrix2;
  int size;
  bool symmetric;
  void init(int _size, double *_points1, double *_points2, bool _symmetric)
    {
      size=_size;
      symmetric=_symmetric;
      CmiMemcpy(matrix1,_points1,size*sizeof(double));
      CmiMemcpy(matrix2,_points2,size*sizeof(double));
    }
  friend class CMessage_entireResultMsg2;
};




class PairCalculator: public CBase_PairCalculator {
 public:
  PairCalculator(bool sym, int grainSize, int s, int blkSize, CkCallback cb,  CkArrayID final_callbackid, int final_callback_ep, int callback_ep_tol, int callback_rdma_ep, int conserveMemory, bool lbpaircalc, redtypes reduce, int orthoGrainSize, bool _AllTiles, bool streambw, bool delaybw, int streamFW, bool gSpaceSum, int gpriority, bool phantomSym, bool useBWBarrier, int _gemmSplitFWk, int _gemmSplitFWm, int _gemmSplitBW, bool expectOrthoT);
    
  PairCalculator(CkMigrateMessage *);
  ~PairCalculator();

  void lbsync() {
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d] atsyncs\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z);
#endif
    resumed=false;
    rck=0;
    AtSync();
  };
  void bwbarrier(CkReductionMsg *msg)
    {
      // everyone is done
      delete msg;
      // figure out how to send the results from here sanely
      sendBWsignalMsg *sigmsg;
      if(PCdelayBWSend)
	sigmsg= new (8*sizeof(int)) sendBWsignalMsg;
      else
	sigmsg= new  sendBWsignalMsg;
      //collapse this into 1 flag
      bool unitcoef=true;  //cheap hack for minimzation only
      //collapse this into 1 flag
      if(amPhantom)
	sigmsg->otherdata= true;
      else if(((!phantomSym && symmetric) || !unitcoef) && (thisIndex.x != thisIndex.y))
	sigmsg->otherdata=true;       
      else
	sigmsg->otherdata= false;


      if(PCdelayBWSend)
	{
	  CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(sigmsg) = 1; // just make it slower
					    // than non prioritized
	}
      if(gSpaceSum)
	thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResultDirect(sigmsg);
	else
	  thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);


    }

  void multiplyForwardStream(bool flag_dp);
  void sendTiles(bool flag_dp);
  void collectTile(bool doMatrix1, bool doMatrix2, bool doOrthoT, int orthoX, int orthoY, int orthoGrainSizeX, int orthoGrainSizeY, int numRecdBW, int matrixSize, double *matrix1, double* matrix2);
  void multiplyForward(bool flag_dp);
  void multiplyForwardRDMA( sendFWRDMAsignalMsg *msg){multiplyForward(msg->flag_dp);}
  void contributeSubTiles(double *fullOutput);
  void ResumeFromSync();
  void initGRed(initGRedMsg *msg);
  void acceptPairData(calculatePairsMsg *msg);
  void acceptPhantomData(phantomMsg *msg);
  void acceptOrthoT(multiplyResultMsg *msg);
  void multiplyResult(multiplyResultMsg *msg);
  void multiplyPsiV();
  void bwMultiplyDynOrthoT();
  void receiveRDMASenderNotify(int senderProc, int sender, bool fromRow, int size, int totalsize);
  void multiplyResultI(multiplyResultMsg *msg);
  void bwMultiplyHelper(int size, double *matrix1, double *matrix2, double *amatrix, double *amatrix2, bool unitcoef, int m_in, int n_in, int k_in, int BNAoffset, int BNCoffset, int BTAoffset, int BTCoffset, int orthoX, int orthoY, double beta, int ogx, int ogy);
  void bwSendHelper(int orthoX, int orthoY, int sizeX, int sizeY, int ogx, int ogy);
  void sendBWResult(sendBWsignalMsg *msg);
  void sendBWResultDirect(sendBWsignalMsg *msg);
  void sendBWResultColumn(bool other, int startGrain, int endGrain);
  void sendBWResultColumnDirect(bool other, int startGrain, int endGrain);
  void initResultSection(initResultMsg *msg);
  void pup(PUP::er &);
  void reorder(int *offsetMap, int *revOffsetMap, double *data, double *scratch);
  void dumpMatrixDouble(const char *, double *,int,int, int xstart=0,int ystart=0, int xtra1=0, int xtra2=0);
  void dumpMatrixComplex(const char *, complex *,int,int, int xstart=0, int ystart=0, int iter=0);
  void dgemmSplitBwdM(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, double *B, double *bt, double *C);
  void dgemmSplitFwdStreamMK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *C, int *ldc);
  void dgemmSplitFwdStreamNK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *C, int *ldc);
  static void Wrapper_To_CallBack(void* pt2Object){
    // explicitly cast to a pointer to CLA_Matrix
    PairCalculator* mySelf = (PairCalculator*) pt2Object;
    mySelf->checkRDMADone();

  }
  void checkRDMADone()
  {
#ifdef _PAIRCALC_DEBUG_RDMA_
    CkPrintf("pairCalc[%d %d %d %d %d] got rdma, count=%d numExpected=%d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numRecd, numExpected);
#endif
    bool streamready=false;
    numRecd++;
    bool flag_dp=symmetric;
    bool doPsiV=false;
    if(symmetricOnDiagonal) //only left
      streamready=((streamCaughtL>=streamFW) && (streamFW>0));
    else 
      {
	streamready=
	  // left or right has streamFW 
	  ((streamCaughtL>=streamFW)||(streamCaughtR>=streamFW)) 
	  // total left and total right have at least stream FW
	  && ((numRecLeft>=streamFW) && (numRecRight>=streamFW) 
	      && (streamFW>0));
	//      CkPrintf("[%d,%d,%d,%d,%d] streamFW %d streamReady %d streamCL %d streamCR %d RL %d RR %d TR %d NE %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric, streamFW, streamready , streamCaughtL, streamCaughtR, numRecLeft, numRecRight, numRecd, numExpected);

      }

    // call member funcs to determine if we're ready to multiply
    if(streamready || ((streamFW>0) && (numRecd == numExpected) && (!doPsiV) ))
      {
	multiplyForwardStream(flag_dp);	
	// not yet supported for dynamic psiV
      }
    else if (numRecd == numExpected)
      {
	if(!doPsiV)
	  {  //normal behavior
	    actionType=0;
	    if(!expectOrthoT || numRecdBWOT==numOrtho)
	      {
		//		multiplyForward(flag_dp);	
		// needs to become a regular charm message
		  sendFWRDMAsignalMsg *sigmsg=new (8*sizeof(int)) sendFWRDMAsignalMsg;
		  // needs to be prioritized and go through the usual scheduler
		  CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
		  *(int*)CkPriorityPtr(sigmsg) = 300000000; // lambda prio
		  sigmsg->flag_dp=flag_dp;
		  thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).multiplyForwardRDMA(sigmsg);	
	      }
	    else
	      {
		if(expectOrthoT)
		  CkPrintf("GAMMA BEAT ORTHOT, holding\n");
	      }
	    if(expectOrthoT && numRecdBWOT==numOrtho)
	      { // we must also multiply orthoT by Fpsi
		bwMultiplyDynOrthoT();
		//	      CkPrintf("[%d,%d,%d,%d,%d] cleanup numRecdBWOT now %d \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numRecdBWOT);
		numRecdBWOT=0;
	      }
	  }
	else
	  {
	    // tolerance correction psiV
	    multiplyPsiV();
	  }
      }
    else
      {
	//      CkPrintf("[%d,%d,%d,%d,%d] no fwd yet numRecd %d numExpected %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric, numRecd, numExpected);
      }
  }
    


 private:
  RDMAHandle *RDMAHandlesRight;
  RDMAHandle *RDMAHandlesLeft;
  int numRecd;               //! number of messages received
  int numRecdBW;               //! number of messages received BW
  int numRecdBWOT;               //! number of messages received BW orthoT
  int numExpected;           //! number of messages expected all
  int numExpectedX;           //! number of messages expected x-axis
  int numExpectedY;           //! number of messages expected y-axis
  int grainSize;             //! number of states per chare
  int grainSizeX;             //! number of states per chare x-axis
  int grainSizeY;             //! number of states per chare y-axis
  int orthoGrainSize;        //! number of states per ortho tile lower-bound
  int orthoGrainSizeRemX;    //! sgrainSizeX % orthoGrainSize
  int orthoGrainSizeRemY;    //! sgrainSizeY % orthoGrainSize
  int blkSize;               //! number points in gspace plane
  int numStates;             //! total number of states
  int numPoints;             //! number of points in this chunk
  int numChunks;             //! number of blocks the stateplane is
			     //divided into

  int streamFW;              //! number of rows to accumulate before
			     //computing
  
  int streamCaughtR;          //! number of rows caught so far R stream
  int streamCaughtL;          //! number of rows caught so far L stream
  
  int numRecLeft;           //! number of rows so far total left
  int numRecRight;          //! number of rows so far total right

  int gemmSplitFWk;        //! number of rows in split FW dgemm
  int gemmSplitFWm;        //! number of columns in split FW dgemm
  int gemmSplitBW;        //! number of rows in split BW dgemm


  int *LeftOffsets;           //! index numbers of caught stream elements
  int *RightOffsets;           //! index numbers of caught stream elements

  int *LeftRev;           //! reverse index numbers of caught stream elements
  int *RightRev;           //! reverse index numbers of caught stream elements

  double **outTiles;         //! in output streaming we populate the
			     //! tiles directly

  int *touchedTiles;         //! tracker to detect when tiles are full
  
  bool symmetric;            //! if true, one triangle is missing
  int conserveMemory;       //! free up matrices when not in use
  bool lbpaircalc;
  bool notOnDiagonal;              //! being on or off diagonal changes many things
  bool symmetricOnDiagonal;     //! diagonal symmetric special case

  bool phantomSym;           //! phantoms exist to balance the BW path
			     //otherdata work

  bool expectOrthoT;         //! orthoT should arrive before end of
			     //  fwd path
  bool amPhantom;            //! consolidate thisIndex.x<thisIndex.y && symmetric && phantomsym
  
  bool useBWBarrier;
  
  bool collectAllTiles;      //! If true, don't stream compute on tiles in the backward path.
   
  redtypes cpreduce;         //! which reducer we're using (defunct)
  CkArrayID cb_aid;          //! bw path callback array ID 
  int cb_ep;                 //! bw path callback entry point 
  int rdma_ep;               //! rdma setup callback entry point 
  int cb_ep_tol;             //! bw path callback entry point for psiV tolerance
  bool existsLeft;           //! inDataLeft allocated 
  bool existsRight;          //! inDataRight allocated 
  bool existsOut;            //! outData allocated
  bool existsNew;            //! newData allocated
  bool resumed;              //! have resumed from load balancing

  bool PCstreamBWout;        //! stream output from BW path       
  bool PCdelayBWSend;        //! use priority to delay BW output 

  bool gSpaceSum;            //! sum in gspace instead of reduction
  int gpriority;            //! priority of msg to gspace

  complex *mynewData;        //! results of bw multiply
  complex *othernewData;     //! results of sym off diagonal multiply,
                             //! or the C=-1 *inRight* orthoT +c in dynamics
  double *inDataLeft;        //! the input pair to be transformed
  double *inDataRight;       //! the input pair to be transformed
  double *outData;           //! results of fw multiply
  int actionType;            //! matrix usage control [NORMAL, KEEPORTHO, PSIV]

  double *allCaughtLeft;     //! unordered rows of FW input
  double *allCaughtRight;    //! unordered rows of FW input
  

  double *inResult1;         //! accumulate ortho or lambda
  double *inResult2;         //! used in gamma calc (non minization)

  /* to support the simpler section reduction*/
  int rck;                   //! count of received cookies
  CkGroupID mCastGrpIdOrtho;  //! group id for multicast manager ortho
  int numOrthoCol;            //! sGrainSizeX/orthoGrainSize
  int numOrthoRow;            //! sGrainSizeY/orthoGrainSize
  int numOrtho;               //! number of orthos in our grain = numOrthoCol*numOrthoRow

  CkGroupID mCastGrpId;      //! group id for multicast manager bw

  CkSectionInfo *resultCookies;  //! array of bw path section cookies
  CkSectionInfo *otherResultCookies;  //! extra array of bw path
                                      //! section cookies
                                      //! for sym off diag, or dynamics

  CkCallback *orthoCB;             //! forward path callbacks
  CkSectionInfo *orthoCookies;      //! forward path reduction cookie 
  int *columnCount;                  //! count of processed rows in BW
				    // by column 
  int *columnCountOther;             //! count of processed rows in BW
				    //by column
// copy the results from outdata1 and outdata2 into the tiles
/**
 * Iterate through the source array, look up the destination row in
 * offsetsRow, destination col in offsetsCol
 *
 * This will be the destination row and column for the output if the
 * output were considered as a single matrix.
 *
 * Use tileSize to map these values into the destination tile.
 *
 */
  void copyIntoTiles(double *source, double**dest, int sourceRows, int sourceCols, int *offsetsRow, int *offsetsCol, int *touched, int tileSize, int tilesPerRow );

};

//forward declaration
CkReductionMsg *sumMatrixDouble(int nMsg, CkReductionMsg **msgs);
CkReductionMsg *sumBlockGrain(int nMsg, CkReductionMsg **msgs);

void manmult(int numrowsA, int numRowsB, int rowLength, double *A, double *B, double *C, double alpha);
#endif
