#include "pcMessages.decl.h"
#include "ckcomplex.h"

#ifndef PC_MESSAGES_H
#define PC_MESSAGES_H

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
    CkIndex4D sndr;
    void init(CkIndex4D _sndr, int _size, int _myoffset, complex *_points)
    {
      sndr = _sndr;
      N=_size;
      myoffset=_myoffset;
      CmiMemcpy(result,_points,N*sizeof(complex));
    }

    friend class CMessage_partialResultMsg;

};



class priorSumMsg : public CMessage_priorSumMsg {
  public:
    complex *result;
    int N;
    int priority;
    CkCallback cb;

    friend class CMessage_priorSumMsg;

};



/** The new message for sending input data to the PairCalculator. Implements API as required by MessageDataCollator viz.
 * the sender(), numDataUnits() and data() methods.
 */
class paircalcInputMsg: public CkMcastBaseMsg, public CMessage_paircalcInputMsg
{
  friend class CMessage_paircalcInputMsg;

  public:
  /// Provide the API required by the MessageDataCollator
  /// An incomplete constructor just for the collator
  paircalcInputMsg(const int _sender, const int _nCols, const int _nRows=1)
    : nCols(_nCols), nRows(_nRows), senderID(_sender), doPsiV(false) {}
  /// The number of rows in the data array stored in this msg
  inline int numRows() const { return nRows; }
  /// The number of columns of data units in the data array stored in this msg
  inline int numCols() const { return nCols; }
  /// An integer representation of the sender's ID
  inline int sender() const  { return senderID; }
  /// A pointer to the message data. No checks on pointer validity. Use with a pinch of salt
  inline inputType* data()   { return points; }
  /// Constructor used to create actual GSpace to PC messages
  paircalcInputMsg(int _size, int _sender, bool _fromRow, bool _flag_dp, inputType *_points , bool _doPsiV, int _blkSize, int _nRows=1)
  {
    nCols    =_size;
    nRows    =_nRows;
    senderID =_sender;
    fromRow  =_fromRow;
    flag_dp  =_flag_dp;
    doPsiV   =_doPsiV;
    blkSize  =_blkSize;
    CmiMemcpy(points,_points,nRows*nCols*sizeof(inputType));
  }
  /// @todo: Message data, should slowly be hidden from the world. The sender and end user could become friends
  inputType *points;
  bool fromRow, flag_dp, doPsiV;
  int blkSize; ///< @todo: blkSize is used in paircalc only for dumping data files in the backward path. Also, it cannot be retreived when usng RDMA. Is there an alternative?

  private:
  int senderID, nCols, nRows;
};


/*
   class phantomMsg : public CMessage_phantomMsg {
   public:
   int size;
   int numPoints;
   double *points;
   int blkSize;
   int actionType;
   void init(int _size, int _numPoints, bool _flag_dp, double *_points, int _blkSize, int _actionType)
   {
   CkPrintf("WARNING: WARNING: phantomMsg has been created\n");
   size=_size;
   numPoints=_numPoints;
   blkSize=_blkSize;
   CmiMemcpy(points,_points,size*sizeof(double));
   actionType=_actionType;
   }

   };

 */

class multiplyResultMsg : public CkMcastBaseMsg, public CMessage_multiplyResultMsg {
  public:
    internalType *matrix1;
    internalType *matrix2;
    int size;
    int size2;
    int orthoX;
    int orthoY;
    int actionType;
    void init(int _size, int _size2, internalType *_points1, internalType *_points2, int _orthoX, int _orthoY, bool _actionType)
    {
      size=_size;
      size2=_size2;
      orthoX=_orthoX;
      orthoY=_orthoY;
      CmiMemcpy(matrix1,_points1,size*sizeof(internalType));
      CmiMemcpy(matrix2,_points2,size2*sizeof(internalType));
      actionType=_actionType;
    }
    void init1(int _size, internalType *_points1, int _orthoX, int _orthoY,int _actionType)
    {
      size=_size;
      size2=0;
      orthoX=_orthoX;
      orthoY=_orthoY;
      CmiMemcpy(matrix1,_points1,size*sizeof(internalType));
      actionType=_actionType;
      // this field does nothing in minimization
      matrix2=NULL;
    }
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

#endif // PC_MESSAGES_H

