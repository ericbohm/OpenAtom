/** \file ckPairCalculator.h
 *
 */
//#define _PAIRCALC_NAN_CHECK_
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
#include "MultiRingMulticast.h"
#include "NodeMulticast.h"
#include "debug_flags.h"
// Flag to use sparse reduction or regular reduction

// Debugging flag for Verbose output
//#define _PAIRCALC_DEBUG_
#ifdef CMK_VERSION_BLUEGENE
#define ALIGN16(x)        (int)((~15)&((x)+15))
//#define TEST_ALIGN
#define BUNDLE_USER_EVENT  
#define PC_FWD_DGEMM_SPLIT 16   //multiple of 6 for BG/L?  use 16 for happier align 
#define PC_BWD_DGEMM_SPLIT 16
#else
#define PC_FWD_DGEMM_SPLIT 0
#define PC_BWD_DGEMM_SPLIT 0
#endif

//flags to control semantic for matrix contents
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

//extern "C" void DGEMM(char *,char *, int *,int *, int *,double *,double *,int *, double *,int *,double *,double *,int *);
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
      memcpy(points,_points,size*sizeof(complex));
    }
  friend class CMessage_calculatePairsMsg;

};

class phantomMsg : public CMessage_phantomMsg {
 public:
  int size;
  int numPoints;
  double *points;
  int blkSize;
  void init(int _size, int _numPoints, bool _flag_dp, double *_points, int _blkSize)
    {
      size=_size;
      numPoints=_numPoints;
      blkSize=_blkSize;
      memcpy(points,_points,size*sizeof(double));
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
      memcpy(matrix1,_points1,size*sizeof(double));
      memcpy(matrix2,_points2,size2*sizeof(double));
      actionType=_actionType;
    }
  void init1(int _size, double *_points1, int _orthoX, int _orthoY,int _actionType)
    {
      size=_size;
      size2=0;
      orthoX=_orthoX;
      orthoY=_orthoY;
      memcpy(matrix1,_points1,size*sizeof(double));
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
      memcpy(matrix,_points,size*sizeof(double));
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
      memcpy(matrix1,_points1,size*sizeof(double));
      memcpy(matrix2,_points2,size*sizeof(double));
    }
  friend class CMessage_entireResultMsg2;
};




class PairCalculator: public CBase_PairCalculator {
 public:
  PairCalculator(bool sym, int grainSize, int s, int blkSize, CkCallback cb,  CkArrayID final_callbackid, int final_callback_ep, int callback_ep_tol, bool conserveMemory, bool lbpaircalc, redtypes reduce, int orthoGrainSize, bool _AllTiles, bool streambw, bool delaybw, int streamFW, bool gSpaceSum, int gpriority, bool phantomSym, bool useBWBarrier, int _gemmSplitFWk, int _gemmSplitFWm, int _gemmSplitBW);
    
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
  void multiplyForward(bool);
  void contributeSubTiles(double *fullOutput);
  void ResumeFromSync();
  void initGRed(initGRedMsg *msg);
  void acceptPairData(calculatePairsMsg *msg);
  void acceptPhantomData(phantomMsg *msg);
  void sendBWResult(sendBWsignalMsg *msg);
  void sendBWResultDirect(sendBWsignalMsg *msg);
  void sendBWResultColumn(bool other, int startGrain, int endGrain);
  void sendBWResultColumnDirect(bool other, int startGrain, int endGrain);
  void multiplyResult(multiplyResultMsg *msg);
  void multiplyPsiV();
  void multiplyResultI(multiplyResultMsg *msg);
  void initResultSection(initResultMsg *msg);
  void pup(PUP::er &);
  void dgemmSplitBwdM(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, double *B, double *bt, double *C);
  void dgemmSplitFwdStreamMK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *bt, double *C, int *ldc);
  void dgemmSplitFwdStreamNK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *bt, double *C, int *ldc);
 
  void dumpMatrixDouble(const char *, double *,int,int, int xstart=0,int ystart=0 );
  void dumpMatrixComplex(const char *, complex *,int,int, int xstart=0, int ystart=0);

 private:

  int numRecd;               //! number of messages received
  int numRecdBW;               //! number of messages received BW
  int numExpected;           //! number of messages expected
  int grainSize;             //! number of states per chare
  int orthoGrainSize;        //! number of states per ortho tile
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
  bool conserveMemory;       //! free up matrices when not in use
  bool lbpaircalc;
  bool phantomSym;           //! phantoms exist to balance the BW path
			     //otherdata work

  bool amPhantom;            //! consolidate thisIndex.x<thisIndex.y && symmetric && phantomsym
  
  bool useBWBarrier;
  
  bool collectAllTiles;      //! If true, don't stream compute on tiles in the backward path.
   
  redtypes cpreduce;         //! which reducer we're using (defunct)
  CkArrayID cb_aid;          //! bw path callback array ID 
  int cb_ep;                 //! bw path callback entry point 
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


#endif
