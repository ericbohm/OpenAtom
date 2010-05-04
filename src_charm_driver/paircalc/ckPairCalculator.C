#include "ckPairCalculator.h"
#include "utility/matrix2file.h"
#include <sstream> 

ComlibInstanceHandle mcastInstanceCP;
ComlibInstanceHandle mcastInstanceACP;

CkReduction::reducerType sumMatrixDoubleType;


void registersumMatrixDouble(void)
{
  sumMatrixDoubleType=CkReduction::addReducer(sumMatrixDouble);
}
/*CkReduction::reducerType sumFastDoubleType;
CkReductionMsg *sumFastDouble(int nMsg, CkReductionMsg **msgs);*/
void fastAdd (double *a, double *b, int nelem);


// sum together matrices of doubles
// possibly faster than CkReduction::sum_double due to minimizing copies
// and calling CmiNetworkProgress
inline CkReductionMsg *sumMatrixDouble(int nMsg, CkReductionMsg **msgs)
{
  double *ret=(double *)msgs[0]->getData();

  //  CkAssert ((unsigned int) ret % 8 == 0);
#ifdef CMK_BLUEGENEL
  //      __alignx(16,ret);
#endif
  int size0=msgs[0]->getSize();
  int size=size0/sizeof(double);

  double *inmatrix;
  //  int progcount=0;
  if(nMsg>3) // switch loops and unroll
    {
      int i=1;
      // idea here is to have only 1 store for 4 loads
#ifdef CMK_BLUEGENEL
#pragma unroll(10)
      // how much doth XLC sucketh?
#endif
      for(int d=0;d<size;d++)
	{
	  for(i=1; i<nMsg-3;i+=3)
	    ret[d]+= ((double *) msgs[i]->getData())[d] + ((double *) msgs[i+1]->getData())[d] + ((double *) msgs[i+2]->getData())[d];
	  for(; i<nMsg;i++)
	    {
	      ret[d]+=((double *) msgs[i]->getData())[d];
	    }
	}
    }
  else
  for(int i=1; i<nMsg;i++)
    {

      inmatrix=(double *) msgs[i]->getData();
#ifdef CMK_BLUEGENEL
      //      __alignx(16,inmatrix);
#pragma disjoint(*ret,*inmatrix)
#pragma unroll(16)
#endif
	for(int d=0;d<size;d++)
	  ret[d]+=inmatrix[d];
    }
  //  CmiNetworkProgress();
  return CkReductionMsg::buildNew(size*sizeof(double),ret);
}



// A functor to simply delegate a gemm to either zgemm or dgemm based on how its instantiated
class gemmDelegator
{
    public:
        /// Construct a gemm delegate functor
        gemmDelegator(bool useZgemm = false): useComplex(useZgemm)
        {
            CkAssert(sizeof(complex)/sizeof(double) == 2);
        }

        /// Use it just like a GEMM call
        inline void operator() (char *opA, char *opB, int *m, int *n, int *k, double *alpha, complex *A, int *lda, complex *B, int *ldb, double *beta, complex *C, int *ldc)
        {
            if (useComplex)
                //ZGEMM(opA, opB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
                CkAbort("zgemm is not yet available");
            else
            {
                DGEMM(opA, opB, m, n, k, alpha, reinterpret_cast<double*>(A), lda, reinterpret_cast<double*>(B), ldb, beta, reinterpret_cast<double*>(C), ldc);
            }
        }

    private:
        bool useComplex;
};



PairCalculator::PairCalculator(CkMigrateMessage *m) { }

PairCalculator::PairCalculator(CProxy_InputDataHandler<CollatorType,CollatorType> inProxy, const pc::pcConfig _cfg): cfg(_cfg)
{
#ifdef _PAIRCALC_DEBUG_PLACE_
  CkPrintf("{%d} [PAIRCALC] [%d,%d,%d,%d,%d] inited on pe %d \n", _instance,thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,_sym, CkMyPe());
#endif


  int remainder        = cfg.numStates%cfg.grainSize;
  grainSizeX=(cfg.numStates- thisIndex.x == cfg.grainSize+remainder) ? cfg.grainSize+remainder: cfg.grainSize;
  grainSizeY=(cfg.numStates- thisIndex.y == cfg.grainSize+remainder) ? cfg.grainSize+remainder: cfg.grainSize;

  this->numPoints = -1;
  this->cb_aid         = cfg.gSpaceAID;
  this->cb_ep          = cfg.gSpaceEP;
  this->cb_ep_tol      = cfg.PsiVEP;
  orthoGrainSizeRemX=grainSizeX%cfg.orthoGrainSize;
  orthoGrainSizeRemY=grainSizeY%cfg.orthoGrainSize;
  gemmSplitFWk         = cfg.gemmSplitFWk;
  gemmSplitFWk         = cfg.gemmSplitFWm;
  gemmSplitBW          = cfg.gemmSplitBW;
  existsLeft=false;
  existsRight=false;
  existsOut=false;
  existsNew=false;
  numRecd = 0;
  numRecdBW = 0;
  numRecdBWOT = 0;
  numRecRight = 0;
  numRecLeft = 0;
  streamCaughtR=0;
  streamCaughtL=0;
  expectOrthoT= (cfg.isDynamics && !cfg.isSymmetric);
  amPhantom=(cfg.arePhantomsOn && (thisIndex.y<thisIndex.x) && cfg.isSymmetric ) ? true : false;
  /*  if(amPhantom)
    { // ye old switcheroo for the phantoms
      numExpectedX=grainSizeY;
      grainSizeY=grainSizeX;
      grainSizeX=numExpectedX;
    }
  */
  numExpectedX = grainSizeX;
  numExpectedY = grainSizeY;
  numExpected = numExpectedX + numExpectedY;

  notOnDiagonal= (thisIndex.x!=thisIndex.y) ? true: false;
  symmetricOnDiagonal=(cfg.isSymmetric && thisIndex.x==thisIndex.y) ? true: false;
  // If I lie on the chare array diagonal, I expect to get only left matrix data
  if(symmetricOnDiagonal)
    numExpected=numExpectedX;
  // If I am a phantom chare, I expect to get only right matrix data
  if(amPhantom)
  {
      ///@todo: This is a hack to ensure that phantoms which get matrix blocks
      // which should actually be *their* left matrix blocks but arrive as right
      // data dont choke when running remaindery grainSizes. Read commit log for info
      numExpected = numExpectedX;
      numExpectedX= numExpectedY;
      numExpectedY= numExpected;
  }
  resumed=true;

  touchedTiles=NULL;
  msgLeft = msgRight = 0;
  inDataLeft = NULL;
  inDataRight = NULL;
  allCaughtLeft=NULL;
  allCaughtRight=NULL;
  outData = NULL;
  mynewData= NULL;
  othernewData= NULL;
  inResult1=NULL;
  inResult2=NULL;
  usesAtSync=true;
  if(cfg.isLBon)
      setMigratable(true);
  else
      setMigratable(false);
  resultCookies=NULL;
  otherResultCookies=NULL;
  // TODO: technically we can make fewer of these if cfg.grainSize>grainSizeX || grainSizeY but you'll only save a few bytes
  resultCookies=new CkSectionInfo[numExpectedX];
  numOrthoCol=grainSizeX/cfg.orthoGrainSize;
  numOrthoRow=grainSizeY/cfg.orthoGrainSize;
  numOrtho=numOrthoCol*numOrthoRow;
  orthoCookies=new CkSectionInfo[numOrtho*2];
  orthoCB=new CkCallback[numOrtho*2];
  if(cfg.isBWstreaming && !cfg.areBWTilesCollected)
    {
      columnCount= new int[numOrthoCol];
      columnCountOther= new int[numOrthoCol];
      bzero(columnCount, sizeof(int) * numOrthoCol);
      bzero(columnCountOther, sizeof(int) * numOrthoCol);
    }
  else
    {
      columnCountOther=NULL;
      columnCount=NULL;
    }
  if(notOnDiagonal)
    // we don't actually use these in the asymmetric minimization case
    // but we make them anyway
    otherResultCookies=new CkSectionInfo[numExpectedY];

	/** -------- Setup the forward path input message handling -------- **/
	/// Set isDataReady flags to false only for those inputs this chare will be getting
    if (amPhantom)
    {
        /// Phantoms will get only a right matrix input from their mirror chares
        isLeftReady  = true;
        isRightReady = false;
    }
    else
    {
        /// Non-phantoms will always get at least a left matrix input
        isLeftReady = false;
        /// Chares on the chare array diagonals will get ONLY a left matrix (which will also be the right). 
        if (symmetricOnDiagonal)
            isRightReady = true;
        else
            isRightReady = false;
    }
	/// Setup the callbacks 
    CkArrayIndex4D myIndex(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z);
    CkCallback leftTrigger (CkIndex_PairCalculator::acceptLeftData (NULL),myIndex,thisProxy,true);
    CkCallback rightTrigger(CkIndex_PairCalculator::acceptRightData(NULL),myIndex,thisProxy,true);
	/// Create a string that holds the chare ID and pass it to the message handlers
	std::ostringstream idStream;
	idStream<<"["<<thisIndex.w<<","<<thisIndex.x<<","<<thisIndex.y<<","<<thisIndex.z<<","<<cfg.isSymmetric<<"]";
	/// Create the message handlers for the left and right input matrix blocks
	leftCollator = new CollatorType (idStream.str()+" LeftHandler" , leftTrigger, numExpectedX,(cfg.conserveMemory<=0),thisIndex.x);
	rightCollator= new CollatorType (idStream.str()+" RightHandler",rightTrigger, numExpectedY,(cfg.conserveMemory<=0),thisIndex.y);
	#ifdef DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
		CkPrintf("[%d,%d,%d,%d,%d] My left and right data collators: %p %p\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,leftCollator,rightCollator);
	#endif
	/// This is the first point during execution when I can supply my InputDataHandler with pointers to the msg handlers, hence
	/// it is (now) safe to insert the [w,x,y,z]th element of the InputDataHandler chare array (as it will immediately clamor 
	/// for access to these handlers)
	myMsgHandler = inProxy;
	#ifdef DEBUG_CP_PAIRCALC_CREATION
		CkPrintf("[%d,%d,%d,%d,%d] Inserting my InputDataHandler\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
	#endif
	myMsgHandler(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).insert(thisProxy);
	myMsgHandler.doneInserting();
}



void
PairCalculator::pup(PUP::er &p)
{
  ArrayElement4D::pup(p);
  p|numRecd;
  p|numRecdBW;
  p|numRecdBWOT;
  p|numExpected;
  p|numExpectedX;
  p|numExpectedY;
  p|grainSizeX;
  p|grainSizeY;
  p|numPoints;
  p|symmetricOnDiagonal;
  p|notOnDiagonal;
  p|cb_aid;
  p|cb_ep;
  p|cb_ep_tol;
  p|existsLeft;
  p|existsRight;
  p|existsOut;
  p|existsNew;
  p|mCastGrpId;
  p|mCastGrpIdOrtho;
  p|resumed;
  p|rck;
  p|actionType;
  p|orthoGrainSizeRemX;
  p|orthoGrainSizeRemY;
  p|expectOrthoT;
  p|amPhantom;
  p|numOrthoCol;
  p|numOrthoRow;
  p|numOrtho;
  if (p.isUnpacking())
  {
      mynewData=NULL;
      othernewData=NULL;
      resultCookies=new CkSectionInfo[cfg.grainSize];
      if(notOnDiagonal)
	otherResultCookies= new CkSectionInfo[cfg.grainSize];
      else
	otherResultCookies=NULL;
      if(existsOut)
	  outData= new double[grainSizeX*grainSizeY];
      else
	  outData=NULL;
      /// @todo: Fix this to allocate or grab a msgLeft and msgRight. inDataLeft/Right is no longer allocated directly
      if(existsLeft)
	  inDataLeft = reinterpret_cast<double*> ( new inputType[numExpectedX*numPoints] );
      else
	  inDataLeft=NULL;
      if(existsRight)
	  inDataRight = reinterpret_cast<double*> ( new inputType[numExpectedY*numPoints] );
      else
	  inDataRight=NULL;
      orthoCookies= new CkSectionInfo[numOrtho*2];
      orthoCB= new CkCallback[numOrtho*2];
      if(cfg.isBWstreaming && !cfg.areBWTilesCollected)
	{
	  columnCount= new int[numOrthoCol];
	  columnCountOther= new int[numOrthoCol];

	}

  }
  int i;
  for (i=0; i<cfg.grainSize; i++) p|resultCookies[i];
  if(notOnDiagonal)
    for (i=0; i<cfg.grainSize; i++) p|otherResultCookies[i];
  for (int i=0; i<cfg.grainSize; i++)
    CmiAssert(resultCookies[i].get_redNo() > 0);
  if(existsOut)
    p(outData, cfg.grainSize * cfg.grainSize);
    /** @todo: Fix this to pack msgLeft and msgRight directly. inDataLeft/Right is no longer allocated directly
      * msgLeft and msgRight will already be packed and available as msgs. They just wont be packed together with the rest of paircalc. 
      * Is there any way we could simply hand the msgs back to charm and ask it to deliver them to me after I am done migrating.?
      * Or am I talking crap?
      */
  if(existsLeft)
    p(inDataLeft, numExpectedX * numPoints * 2);
  if(existsRight)
    p(inDataRight, numExpectedY* numPoints * 2);
  PUParray(p,orthoCookies,numOrtho*2);
  PUParray(p,orthoCB,numOrtho*2);
  if(cfg.isBWstreaming && !cfg.areBWTilesCollected)
    {
      PUParray(p,columnCount,numOrthoCol);
      PUParray(p,columnCountOther,numOrthoCol);
    }
#ifdef _PAIRCALC_DEBUG_
  if (p.isUnpacking())
    {
      CkPrintf("[%d,%d,%d,%d,%d] pup unpacking on %d resumed=%d memory %d\n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,cfg.isSymmetric,CkMyPe(),resumed, CmiMemoryUsage());
      CkPrintf("[%d,%d,%d,%d,%d] pupped : %d,%d,%d,%d,%d %d %d %d %d  %d %d cb cb_aid %d %d %d cb_lb inDataLeft inDataRight outData  %d \n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numRecd, numExpected, cfg.grainSize, cfg.numStates, cfg.numChunks, numPoints, cfg.isSymmetric, cfg.conserveMemory, cfg.isLBon, cfg.reduce, cb_ep, existsLeft, existsRight,  resumed);

    }
  else
    CkPrintf("[%d,%d,%d,%d,%d] pup called on %d\n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,cfg.isSymmetric,CkMyPe());
#endif


}

PairCalculator::~PairCalculator()
{

#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] destructs on [%d]\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, cfg.isSymmetric, CkMyPe());
#endif
  if(outData!=NULL)
    delete [] outData;
  if(msgLeft)
      delete msgLeft;
  if(msgRight)
      delete msgRight;
  if(mynewData!=NULL)
    delete [] mynewData;
  if(othernewData!=NULL)
    delete [] othernewData;
  // redundant paranoia
  outData=NULL;
  inDataRight=NULL;
  inDataLeft=NULL;
  existsLeft=false;
  existsRight=false;
  existsNew=false;
  existsOut=false;
  CkAbort("paircalc pup needs fw streaming work");
  if(orthoCookies!=NULL)
    delete [] orthoCookies;
  if(orthoCB!=NULL)
    delete [] orthoCB;
  if(resultCookies!=NULL)
    delete [] resultCookies;
  if(notOnDiagonal && otherResultCookies !=NULL)
    delete [] otherResultCookies;
}



void PairCalculator::initGRed(initGRedMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================

  int maxorthostateindex=(cfg.numStates/cfg.orthoGrainSize-1)*cfg.orthoGrainSize;
  //  int orthoIndexX=(msg->orthoX*cfg.orthoGrainSize-thisIndex.x)/cfg.orthoGrainSize;
  //  int orthoIndexY=(msg->orthoY*cfg.orthoGrainSize-thisIndex.y)/cfg.orthoGrainSize;
  int orthoIndexX=msg->orthoX*cfg.orthoGrainSize;
  orthoIndexX= (orthoIndexX>maxorthostateindex) ? maxorthostateindex : orthoIndexX;
  int orthoIndexY=msg->orthoY*cfg.orthoGrainSize;
  orthoIndexY= (orthoIndexY>maxorthostateindex) ? maxorthostateindex : orthoIndexY;
  orthoIndexX=(orthoIndexX-thisIndex.x)/cfg.orthoGrainSize;
  orthoIndexY=(orthoIndexY-thisIndex.y)/cfg.orthoGrainSize;

  int orthoIndex=orthoIndexX*numOrthoCol+orthoIndexY;

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] initGRed ox %d oy %d oindex %d oxindex %d oyindex %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, cfg.isSymmetric,msg->orthoX, msg->orthoY,orthoIndex, orthoIndexX, orthoIndexY);
#endif
 // numOrtho here is numOrtho per sGrain
  CkAssert(orthoIndex<numOrtho*2);
  CkGetSectionInfo(orthoCookies[orthoIndex],msg);
  orthoCB[orthoIndex]=msg->cb;
  mCastGrpIdOrtho=msg->mCastGrpId;
  /*  cfg.reduce=section;
  if(msg->lbsync)
  {
      int foo=1;
      contribute(sizeof(int), &foo , CkReduction::sum_int, msg->synccb);
  }

  */

  /// @note: numRecd here is just used as some counter during the init phase. Not related to its usual purpose
  if(!cfg.isSymmetric && ++numRecd==numOrtho)
  {
    struct s_array { //(sendArray, bcastArray)
      int ep; //Entry point to call
      CkGroupID id; //Array ID to call it on
      CkArrayIndexStruct idx; //Index to send to (if any)
    } array, *ap;


      contribute(sizeof(int), &numRecd , CkReduction::sum_int, cfg.uponSetupCompletion, cfg.instanceIndex);
      numRecd=0;
  }

  //  do not delete nokeep msg
}



void PairCalculator::initResultSection(initResultMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================

  CkAssert(msg->offset<cfg.grainSize);
  if(msg->dest == thisIndex.x && thisIndex.x != thisIndex.y)
  {
      CkGetSectionInfo(otherResultCookies[msg->offset],msg);

#ifdef _PAIRCALC_DEBUG_SPROXY_
      CkPrintf("[%d,%d,%d,%d,%d] other initResultSection for dest %d offset %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, cfg.isSymmetric,msg->dest, msg->offset);
#endif
  }
  else
  {
      CkGetSectionInfo(resultCookies[msg->offset],msg);

#ifdef _PAIRCALC_DEBUG_SPROXY_
      CkPrintf("[%d,%d,%d,%d,%d] initResultSection for dest %d offset %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, cfg.isSymmetric,msg->dest, msg->offset);
#endif
  }
  rck++;

  mCastGrpId=msg->mCastGrpId;
  //to force synchronize in lb resumption
  if(msg->lbsync && rck==cfg.grainSize)
  {
      contribute(sizeof(int), &rck, CkReduction::sum_int, msg->synccb);
  }
  // do not delete nokeep msg
}

void PairCalculator::ResumeFromSync() {
  resumed=true;
  // we own no proxies so we have nothing to reset
  if(!resumed){

#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] resumes from sync\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, cfg.isSymmetric);
#endif
  }
}



void PairCalculator::acceptLeftData(paircalcInputMsg *msg) 
{
    complex *data = msg->data();
    const int numRows  = msg->numRows();
    const int numCols  = msg->numCols();
	/// Assert that data is a valid pointer
	CkAssert(data != NULL);
	/// Assert that numRows is as expected
	CkAssert(numRows == numExpectedX);
	/// Check data validity
	#ifdef _NAN_CHECK_
		for(int i=0; i < numRows*numCols; i++)
			CkAssert( finite(data[i].re) && finite(data[i].im) );
	#endif
    /// Once the basic checks have passed, and if we're debugging print status info
	#ifdef _PAIRCALC_DEBUG_
		CkPrintf("[%d,%d,%d,%d,%d] Received left matrix block of size %d x %d at %p\n",
                                thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numRows,numCols,data);
	#endif

	/// Set member data pertinent to the left block
    msgLeft      = msg;
	inDataLeft   = reinterpret_cast<double*> (data);
	existsLeft   = true;
	numRecd     += numRows;
	numPoints    = numCols;
	isLeftReady  = true;

	/// If all data is ready 
	if (isLeftReady && isRightReady)
    {
        // Obtain a sample incoming msg from the collator and extract relevant data from it
        paircalcInputMsg *sampleMsg = leftCollator->getSampleMsg();
        if (sampleMsg)
        {
            msgLeft->doPsiV = sampleMsg->doPsiV;
            msgLeft->blkSize= sampleMsg->blkSize;
            msgLeft->flag_dp= sampleMsg->flag_dp;
        }
        else 
        {
            // If RDMA is enabled, there will be no sample msgs available for any non-PsiV loop. Hence doPsiV is false
            #ifdef PC_USE_RDMA
                msgLeft->doPsiV = false;
            // If RDMA is NOT enabled, we should have obtained a sample message. Something must be wrong
            #else
                std::stringstream dbgStr;
                dbgStr<<"["<<thisIndex.w<<","<<thisIndex.x<<","<<thisIndex.y<<","<<thisIndex.z<<","<<cfg.isSymmetric<<"]"
                    <<" My collator was not able to give me a sample message. Aborting...";
                CkAbort(dbgStr.str().c_str());
            #endif
        }
        delete sampleMsg;
        /// Trigger the computation 
        launchComputations(msgLeft);
    }
}



void PairCalculator::acceptRightData(paircalcInputMsg *msg) 
{
    complex *data = msg->data();
    const int numRows  = msg->numRows();
    const int numCols  = msg->numCols();
	/// Assert that data is a valid pointer
	CkAssert(data != NULL);
	/// Assert that numRows is as expected
	CkAssert(numRows == numExpectedY);
	/// Check data validity
	#ifdef _NAN_CHECK_
		for(int i=0; i < numRows*numCols; i++)
			CkAssert( finite(data[i].re) && finite(data[i].im) );
	#endif
    /// Once the basic checks have passed, and if we're debugging print status info
	#ifdef _PAIRCALC_DEBUG_
		CkPrintf("[%d,%d,%d,%d,%d] Received right matrix block of size %d x %d at %p\n",
                                thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numRows,numCols,data);
	#endif

	/// Set member data pertinent to the right block
    msgRight     = msg;
	inDataRight  = reinterpret_cast<double*> (data);
	existsRight  = true;
	numRecd     += numRows;
	numPoints    = numCols;
	isRightReady = true;

	/// If all data is ready 
	if (isLeftReady && isRightReady)
    {
        /// Phantom chares already have correct info in the incoming msg. Only non-phantoms have to be handled
        if (!amPhantom)
        {
            // Obtain a sample incoming msg from the collator and extract relevant data from it
            paircalcInputMsg *sampleMsg = rightCollator->getSampleMsg();
            if (sampleMsg)
            {
                msgRight->doPsiV = sampleMsg->doPsiV;
                msgRight->blkSize= sampleMsg->blkSize;
                msgRight->flag_dp= sampleMsg->flag_dp;
            }
            else
            {
                // If RDMA is enabled, there will be no sample msgs available for any non-PsiV loop. Hence doPsiV is false
                #ifdef PC_USE_RDMA
                    msgRight->doPsiV = false;
                // If RDMA is NOT enabled, we should have obtained a sample message. Something must be wrong
                #else
                    std::stringstream dbgStr;
                    dbgStr<<"["<<thisIndex.w<<","<<thisIndex.x<<","<<thisIndex.y<<","<<thisIndex.z<<","<<cfg.isSymmetric<<"]"
                        <<" My collator was not able to give me a sample message. Aborting...";
                    CkAbort(dbgStr.str().c_str());
                #endif
            }
            delete sampleMsg;
        }
        /// Trigger the computation 
        launchComputations(msgRight);
    }
}



void PairCalculator::launchComputations(paircalcInputMsg *aMsg)
{
    #ifdef _PAIRCALC_DEBUG_
        CkPrintf("[%d,%d,%d,%d,%d] Going to launch computations... numRecd = %d numExpected = %d numExpectedX = %d, numExpectedY = %d\n",
                                thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,
                                numRecd, numExpected, numExpectedX, numExpectedY);
    #endif
    /// Ensure that we're really ready to launch computations
    CkAssert(numRecd == numExpected);
    blkSize = aMsg->blkSize;
    bool isForwardPathPending = false; ///< Flag to make sure we cleanup only if FW path starts
    
    // If this is not a PsiV loop, trigger the forward path for just the non-phantom chares
    if(!aMsg->doPsiV)
    {
        // This iteration is a normal loop. Hence normal behavior
        actionType = NORMALPC;

        // Start the forward path substep timer
        #ifdef _CP_SUBSTEP_TIMING_
            if(cfg.forwardTimerID > 0)
            {
                double pstart=CmiWallTimer();
                contribute(sizeof(double),&pstart,CkReduction::min_double, cfg.beginTimerCB , cfg.forwardTimerID);
            }
        #endif

        if (!amPhantom)
        {
            /** expectOrthoT is false in any scenario other than asymmetric, dynamics.
            * numRecdBWOT is equal to numOrtho only when it is asymm, dynamics and T has been 
            * received completely (from Ortho). So this condition, invokes multiplyForward() on 
            * all cases except (asymm, dynamics when T has not been received)
            * 
            * In that exception scenario, we dont do anything now. Later, when all of T is received, 
            * both multiplyForward() and bwMultiplyDynOrthoT() are called. Look for these calls in acceptOrthoT().
            */
            if(!expectOrthoT || numRecdBWOT==numOrtho)
                thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).multiplyForward(aMsg->flag_dp);
            else
            {
                isForwardPathPending = true;
                CkPrintf("[%d,%d,%d,%d,%d] Gamma beat OrthoT. Waiting for T to arrive before proceeding with forward path\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
            }
        }
        else
        {
            /// Do nothing for the phantom chare, non-psiv loops. Computation will be triggered only in the backward path
            // Just stop the forward path substep timer for the phantom chares
            #ifdef _CP_SUBSTEP_TIMING_
               if(cfg.forwardTimerID > 0)
               {
                   double pstart=CmiWallTimer();
                   contribute(sizeof(double),&pstart,CkReduction::max_double, cfg.endTimerCB , cfg.forwardTimerID);
               }
            #endif
        }
    }
    // else, if this is a PsiV loop (there is no forward path, only backward path computations to update PsiV)
    else
    {
        // This is a PsiV loop. Hence behave accordingly
        actionType = PSIV;
        thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).multiplyPsiV();
        #ifdef PC_USE_RDMA
            // Let the collators know that they should now expect the next batch of (non-PsiV) data via RDMA
            leftCollator->expectNext();
            rightCollator->expectNext();
        #endif
    }
    
    /// Reset the counters and flags for the next iteration

    /** If asymm,dyn and T has not yet arrived completely, fw path has not yet been triggered. 
     * In this case numRecd will be reset in acceptOrthoT() after the fw path has been triggered. 
     * numRecd is reset here for all other cases.
     */
    if (!isForwardPathPending)
        numRecd = 0;
    /// All non-phantoms should expect left matrix data again
    if (!amPhantom)
        isLeftReady = false;
    /// All non(symm, on-diagonal) chares should expect right matrix data again
    if (!symmetricOnDiagonal)
        isRightReady = false;
}




void PairCalculator::reorder(int * offsetMap, int *revOffsetMap, double *data, double *scratch)
{

  int actualPoints= numPoints*2;
  // rather ugly, but each iteration correctly places 2 rows
  CkAssert(numExpectedY==numExpectedX);
  //TODO HANDLE BORDER grainSizeX|Y funky cases

  int datasize=actualPoints*sizeof(double);
  for(int off=0;off<numExpected;off++)
    {
      if(off!=offsetMap[off])
	{  // gotta move
	  int want = revOffsetMap[off];
	  int found = offsetMap[off];
	  int currentOffset = off * actualPoints;
	  int wantOffset = want * actualPoints;

	  CmiMemcpy(scratch, &(data[currentOffset]), datasize);
	  CmiMemcpy(&(data[currentOffset]), &(data[wantOffset]), datasize);
	  if(want==found) //simple exchange
	    {
	      CmiMemcpy(&(data[wantOffset]),scratch, datasize);
	    }
	  else  // three card shuffle
	    {
	      int foundOffset = found * actualPoints;
	      CmiMemcpy(&(data[wantOffset]), &(data[foundOffset]), datasize);
	      CmiMemcpy(&(data[foundOffset]),scratch, datasize);
	      // 1 more entry is changed
	      offsetMap[want]=offsetMap[found];
	      revOffsetMap[offsetMap[found]]=want;

	    }
	  revOffsetMap[off]=off;
	  offsetMap[off]=off;
	  offsetMap[found]=found;
	  revOffsetMap[found]=found;
	}
    }
  bzero(offsetMap, numExpected*sizeof(int));
  bzero(revOffsetMap, numExpected*sizeof(int));
}

/** @todo: The only use of flag_dp has been commented out. Check if this argument is still needed, else weed it out. */
void
PairCalculator::sendTiles(bool flag_dp)
{
  // TODO: Changes necessary for remainder logic
  // border tiles aren't of uniform size!
  // For now we aren't supporting fw streaming with remainder
  // fw streaming produces incorrect results under as yet unisolated conditions
  CkAssert(orthoGrainSizeRemX==0 && orthoGrainSizeRemY==0);
  int tilesq = cfg.orthoGrainSize * cfg.orthoGrainSize;
  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpIdOrtho).ckLocalBranch();
  //  for(int orthoX=0; orthoX<numOrtho; orthoX++)
  //    for(int orthoY=0; orthoY<numOrtho; orthoY++)

  int progcounter=0;
  for(int orthoIndex=0;orthoIndex<numOrtho;orthoIndex++)
      {
	// copy into submatrix, contribute
	// we need to stride by cfg.grainSize and copy by orthoGrainSize
	//	int orthoIndex=orthoX*numOrthoCol+orthoY;

	if(touchedTiles[orthoIndex]==tilesq)
	  {
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	    if(cfg.isSymmetric)
	      CkPrintf("[%d,%d,%d,%d,%d]: contributes %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, orthoIndex);
#endif
	    //	    if (flag_dp) {
	    //	      for (int i = 0; i < cfg.orthoGrainSize * cfg.orthoGrainSize; i++)
	    //		outTiles[orthoIndex][i] *= 2.0;
	    //	    }
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	    int orthoX=orthoIndex/numOrthoCol;
	    int orthoY=orthoIndex%numOrthoCol;
	    char filename[80];
	    snprintf(filename,80,"fwoutTile_%d_%d:",orthoX,orthoY);
	    dumpMatrixDouble(filename, outTiles[orthoIndex], cfg.orthoGrainSize, cfg.orthoGrainSize,thisIndex.x+orthoX*cfg.orthoGrainSize, thisIndex.y+orthoY*cfg.orthoGrainSize);
#endif

#ifdef _NAN_CHECK_
	    for(int i=0; i<cfg.orthoGrainSize * cfg.orthoGrainSize; i++)
	      CkAssert(finite(outTiles[orthoIndex][i]));
#endif

	    mcastGrp->contribute(cfg.orthoGrainSize * cfg.orthoGrainSize*sizeof(double), outTiles[orthoIndex], sumMatrixDoubleType, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
	    //mcastGrp->contribute(cfg.orthoGrainSize*orthoGrainSize*sizeof(double), outTiles[orthoIndex], CkReduction::sum_double, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
	    touchedTiles[orthoIndex]=0;
	    if(++progcounter>8)
	      {progcounter=0;CmiNetworkProgress();}
	  }
	else if(touchedTiles[orthoIndex]>tilesq)
	  {
	    CkPrintf("tile orthoIndex %d has %d vs %d\n",orthoIndex, touchedTiles[orthoIndex], tilesq);
	    CkAbort("invalid large number of tiles touched");
	  }
	else
	  {
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	      CkPrintf("[%d,%d,%d,%d,%d]: %i not ready with %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, orthoIndex, touchedTiles[orthoIndex]);
#endif
	  }
      }
}


/**
 * Forward path multiply. Essentially does a simple GEMM on the input matrices.
 *
 * We get a left and right matrix each of dimension [numExpectedX/Y, numPoints]
 * We multiply them to get the S matrix of dimension [numExpectedX, numExpectedY]
 *
 * GEMM performs: C = alpha * op(A).op(B) + beta * C
 * where the matrix dimensions are:
 *  op(A) = [m,k]
 *  op(B) = [k,n]
 *     C  = [m,n]
 * with
 *      m = numExpectedX
 *      k = numPoints
 *      n = numExpectedY
 *
 * [numExpectedX, numExpectedY] = [numExpectedX, numPoints] X [numPoints, numExpectedY]
 * To make this work, we transpose the first matrix (A).
 *
 * In C++ it appears to be:         [ydima, ydimb] = [ydima, xdima] X [xdimb, ydimb]
 * which would be wrong, this works because we're using fortran BLAS,
 * which has a transposed perspective (column major).
 * So the actual multiplication is: [xdima, xdimb] = [xdima, ydima] X [ydimb, xdimb]
 *
 * Since xdima == xdimb == numExpected == cfg.grainSize, this gives us the solution matrix we want in one step.
 * In the border case it gives numExpectedX x numExpectedY which is what we want on the borders.
 */
void PairCalculator::multiplyForward(bool flag_dp)
{
    #ifdef _PAIRCALC_DEBUG_
        CkPrintf("[%d,%d,%d,%d,%d] PairCalculator::multiplyForward() Starting forward path computations.\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
    #endif

    // Allocate space for the fw path results if needed
    if(!existsOut)
    {
        CkAssert(outData==NULL);
        outData = reinterpret_cast<double*> (new inputType[grainSizeX * grainSizeY]);
        bzero(outData, sizeof(inputType)* grainSizeX * grainSizeY);
        existsOut=true;
        #ifdef _PAIRCALC_DEBUG_
            CkPrintf("[%d,%d,%d,%d,%d] Allocated outData %d * %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,grainSizeX, grainSizeY);
        #endif
    }

    // Configure the inputs to the GEMM describing the matrix dimensions and operations
    char transformT = 'T';           // Transpose matrix A
    char transform  = 'N';           // Retain matrix B as it is
    int m_in        = numExpectedY;  // Rows of op(A)    = Rows of C
    int k_in        = numPoints;     // Columns of op(A) = Rows of op(B)
    int n_in        = numExpectedX;  // Columns of op(B) = Columns of C
    double alpha    = double(1.0);   // Scale B.A by this scalar factor
    double beta     = double(0.0);   // Scale initial value of C by this factor

    // Get handles to the input and output matrices
    inputType *matrixC = reinterpret_cast<inputType*> (outData);
    inputType *matrixB = msgLeft->data();
    inputType *matrixA;
    if(!symmetricOnDiagonal)
        matrixA = msgRight->data();
    else
    {
        // Symm PC chares on the array diagonal only get a left matrix. For these B serves as A too
        matrixA = matrixB;
        // Redundant, as numExpectedX == numExpectedY (except for the border chares?)
        m_in    = numExpectedX;
    }
    #ifdef TEST_ALIGN
        CkAssert((unsigned int)matrixA%16==0);
        CkAssert((unsigned int)matrixB%16==0);
        CkAssert((unsigned int)matrixC%16==0);
    #endif

    // If we're sending the work to dgemm
    if (!cfg.useComplexMath)
    {
        // Treat each complex as 2 doubles
        k_in *= 2;
        // Double packing (possible only in symm PC for real input) entails a scaling factor for psi
        if (flag_dp)
            alpha = 2.0;
    }


    // Create a gemm delegator which will call the appropriate gemm
    gemmDelegator myGEMM;

    // with dgemm splitting
    #if PC_FWD_DGEMM_SPLIT > 0
        // Results of each smaller gemm must be accumulated in C, not overwritten
        double betap = 1.0;
        // Determine the num of rows of output (C) to compute in each smaller gemm
        int Ksplit_m = gemmSplitFWk;
        int Ksplit   = (k_in > Ksplit_m) ? Ksplit_m : k_in;
        // Calculate the number of gemms that will be required to compute the whole output
        int Krem     = k_in % Ksplit;
        int Kloop    = k_in/Ksplit-1;

        #ifdef PRINT_DGEMM_PARAMS
            CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, Ksplit, alpha, beta, k_in, k_in, m_in);
        #endif

        // Invoke the first split gemm, but with beta=0 so that outData is overwritten (and bracket it in projections)
        #ifndef CMK_OPTIMIZE
            double StartTime=CmiWallTimer();
        #endif
        myGEMM(&transformT, &transform, &m_in, &n_in, &Ksplit, &alpha, matrixA , &k_in, matrixB, &k_in, &beta, matrixC, &m_in);
        CmiNetworkProgress();
        #ifndef CMK_OPTIMIZE
            traceUserBracketEvent(210, StartTime, CmiWallTimer());
        #endif

        // Call each split gemm until all the output is accumulated (beta=1)
        for(int i=1;i<=Kloop;i++)
        {
            int off     = i * Ksplit;
            int KsplitU = (i==Kloop) ? Ksplit+Krem : Ksplit;
            #ifdef TEST_ALIGN
                CkAssert((unsigned int)&(matrixA[off])%16==0);
                CkAssert((unsigned int)&(matrixB[off]) %16==0);
                CkAssert((unsigned int)matrixC%16==0);
            #endif
            #ifdef PRINT_DGEMM_PARAMS
                CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, KsplitU, alpha, beta, k_in, k_in, m_in);
            #endif

            #ifndef CMK_OPTIMIZE
                StartTime=CmiWallTimer();
            #endif
            myGEMM(&transformT, &transform, &m_in, &n_in, &KsplitU, &alpha, &matrixA[off], &k_in, &matrixB[off], &k_in, &betap, matrixC, &m_in);
            CmiNetworkProgress();
            #ifndef CMK_OPTIMIZE
                traceUserBracketEvent(210, StartTime, CmiWallTimer());
            #endif
        } //endfor

    // without dgemm splitting
    #else
        #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
            dumpMatrixDouble("fwlmdata", matrixB, numExpectedX, numPoints*2, thisIndex.x, 0);
            dumpMatrixDouble("fwrmdata", matrixA, numExpectedY, numPoints*2, thisIndex.y, 0);
        #endif
        #ifdef PRINT_DGEMM_PARAMS
            CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, k_in, k_in, m_in);
        #endif

        // Invoke the DGEMM (and bracket it in projections)
        #ifndef CMK_OPTIMIZE
            double StartTime=CmiWallTimer();
        #endif
        myGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, matrixA, &k_in, matrixB, &k_in, &beta, matrixC, &m_in);
        #ifndef CMK_OPTIMIZE
            traceUserBracketEvent(210, StartTime, CmiWallTimer());
        #endif

    #endif  // gemm splitting


    #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
        dumpMatrixDouble("fwgmodata",matrixC,grainSizeX, grainSizeY,thisIndex.x, thisIndex.y);
    #endif
    #ifndef CMK_OPTIMIZE
        StartTime=CmiWallTimer();
    #endif

    // Do the slicing and dicing, and contribute the appropriate bits to the redn that reaches Ortho
    contributeSubTiles(matrixC);
    #ifdef _CP_SUBSTEP_TIMING_
        if(cfg.forwardTimerID > 0)
        {
            double pstart=CmiWallTimer();
            contribute(sizeof(double),&pstart,CkReduction::max_double, cfg.endTimerCB , cfg.forwardTimerID);
        }
    #endif
    #ifndef CMK_OPTIMIZE
        traceUserBracketEvent(220, StartTime, CmiWallTimer());
    #endif

    // Mirror our input data to the phantoms so that they can use it in the bw path
    if(cfg.arePhantomsOn && cfg.isSymmetric && notOnDiagonal)
    {
        CkAssert(existsRight);
        paircalcInputMsg *msg2phantom = new (numExpectedY*numPoints, 8*sizeof(int)) paircalcInputMsg(numPoints,0,false,flag_dp,msgRight->data(),false,blkSize,numExpectedY);
        bool prioPhan=false;
        if(prioPhan)
        {
            CkSetQueueing(msg2phantom, CK_QUEUEING_IFIFO);
            *(int*)CkPriorityPtr(msg2phantom) = 1; // just make it slower than non prioritized
        }
        thisProxy(thisIndex.w,thisIndex.y, thisIndex.x,thisIndex.z).acceptRightData(msg2phantom);
    }

    /** If this is an asymmetric loop, dynamics case AND Ortho has already sent T,
     * call bwMultiplyDynOrthoT() as we must also multiply orthoT by Fpsi
	 *
	 * @note: This if condition originally lived in acceptPairData(). Has been shoveled here
	 * to reduce branching over there.
	 */
    if(expectOrthoT && numRecdBWOT==numOrtho)
        thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).bwMultiplyDynOrthoT();
}




void PairCalculator::contributeSubTiles(inputType *fullOutput)
{
    #ifdef _PAIRCALC_DEBUG_
        CkPrintf("[%d,%d,%d,%d,%d]: contributeSubTiles \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric);
    #endif
    /**
     * Done:
     * necessary changes:
     * 1: border tiles are not uniform size
     * 2: offset calculation must handle non-uniformity
     * solutions, use sGrainSize for all row skips
     * use orthoGrainSizeX or orthoGrainSizeY as needed for row and column iteration
     */

    CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpIdOrtho).ckLocalBranch();
    #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
        dumpMatrixDouble("fullOutput", fullOutput, grainSizeX, grainSizeY, thisIndex.x, thisIndex.y);
    #endif
    inputType *outTile;
    bool reuseTile = false;
    bool borderX   = orthoGrainSizeRemX != 0;
    bool borderY   = orthoGrainSizeRemY != 0;
    bool borderXY  = borderX && borderY;

    // remainder logic happens if you are a border sGrain
    if(! borderX && ! borderY)
    { // only do once to cut down on new/delete
        reuseTile = true;
        outTile  = new inputType[cfg.orthoGrainSize * cfg.orthoGrainSize];
        bzero(outTile,sizeof(inputType) * cfg.orthoGrainSize * cfg.orthoGrainSize);
    }
    // forward multiply ldc
    int bigGindex=grainSizeY;

    for(int orthoX = 0; orthoX < numOrthoCol ; orthoX++)
    {
        // advance tilestart to new column
        // only the size of tiles in last column is affected
        int orthoGrainSizeX=(orthoX == numOrthoCol-1) ? cfg.orthoGrainSize + orthoGrainSizeRemX : cfg.orthoGrainSize;
        int orthoXoff=cfg.orthoGrainSize*bigGindex*orthoX;
        for(int orthoY=0; orthoY<numOrthoRow; orthoY++)
        {
            int orthoYoff=orthoY*cfg.orthoGrainSize;
            int tileStart=orthoYoff+orthoXoff;
            // only the last row is affected
            int orthoGrainSizeY=(orthoY==numOrthoRow-1) ? cfg.orthoGrainSize+orthoGrainSizeRemY : cfg.orthoGrainSize;
            int tileSize=orthoGrainSizeX*orthoGrainSizeY;
            int bigOindex=orthoGrainSizeY;
            int ocopySize=bigOindex*sizeof(inputType);
            int orthoIndex=orthoX*numOrthoCol+orthoY;
            if(! reuseTile)
            {
                outTile = new inputType[tileSize];
                bzero(outTile,sizeof(inputType)*tileSize);
            }
            // copy into submatrix, contribute
            // we need to stride by cfg.grainSize and copy by orthoGrainSize
            CkAssert(orthoIndex<numOrtho);
            for(int ystart=0, itileStart=tileStart; ystart<tileSize; ystart+=bigOindex, itileStart+=bigGindex)
                CmiMemcpy(&(outTile[ystart]),&(fullOutput[itileStart]),ocopySize);

            #ifdef _PAIRCALC_DEBUG_PARANOID_FW_
                char filename[80];
                snprintf(filename,80,"fwoutTile_%d_%d:",orthoX,orthoY);
                dumpMatrixDouble(filename, outTile, orthoGrainSizeX, orthoGrainSizeY,thisIndex.x+orthoX*cfg.orthoGrainSize, thisIndex.y+orthoY*cfg.orthoGrainSize);
            #endif

            mcastGrp->contribute(tileSize*sizeof(inputType), outTile, sumMatrixDoubleType, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
            if(! reuseTile)
                delete [] outTile;
        }
    }
    if(reuseTile)
        delete [] outTile;
}




/** Basically calls collectTile to accept incoming ortho data from an Ortho chare
 * @note This method is used only during dynamics in the asymm loop. Talk to EB for more info
 */
void PairCalculator::acceptOrthoT(multiplyResultMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
  //collect orthoT from matrix2
  CkAssert(expectOrthoT);

  numRecdBWOT++;
  //  CkPrintf("[%d,%d,%d,%d,%d] acceptOrthoT, numRecdBWOT now %d \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numRecdBWOT);
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size;i++)
      CkAssert(finite(msg->matrix1[i]));
#endif

  int size=msg->size;
  int size2=msg->size2;
  double *matrix1=msg->matrix1;
  double *matrix2=msg->matrix2;
#ifdef TEST_ALIGN
  CkAssert((unsigned int) msg->matrix2 %16 ==0);
#endif
  int maxorthostateindex=(cfg.numStates/cfg.orthoGrainSize-1)*cfg.orthoGrainSize;
  // find our tile indices within this sGrain
  int orthoX=msg->orthoX*cfg.orthoGrainSize;

  int orthoY=msg->orthoY*cfg.orthoGrainSize;
  ///? @todo Document this after talking with EB. Shouldnt it be an error if orthoX/Y > maxorthostateindex?
  orthoX= (orthoX>maxorthostateindex) ? maxorthostateindex : orthoX;
  orthoY= (orthoY>maxorthostateindex) ? maxorthostateindex : orthoY;
  orthoX=(orthoX-thisIndex.x)/cfg.orthoGrainSize;
  orthoY=(orthoY-thisIndex.y)/cfg.orthoGrainSize;

  int orthoGrainSizeY=(orthoY==numOrthoRow-1) ? cfg.orthoGrainSize+orthoGrainSizeRemY : cfg.orthoGrainSize;

  int orthoGrainSizeX=(orthoX == numOrthoCol-1) ? cfg.orthoGrainSize + orthoGrainSizeRemX : cfg.orthoGrainSize;
  int matrixSize=grainSizeX*grainSizeY;

  collectTile(false, true, false,orthoX, orthoY, orthoGrainSizeX, orthoGrainSizeY, numRecdBWOT, matrixSize, matrix2, matrix1);
  if ((numRecdBWOT==numOrtho) && (numRecd == numExpected)) 
    { // forward path beat orthoT
      CkPrintf("[%d,%d,%d,%d,%d] Just received all of OrthoT. Triggering forward path multiply which has been pending since Gamma arrival\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
      actionType  = NORMALPC;
      bool myfalse= false;
      thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).multiplyForward(myfalse);
      /** @note: multiplyForward() already checks for numRecdBWOT and calls bwMultiplyDynOrthoT().
       * There is no need to call it again here. It appears, this would have been a correctness bug 
       * that would have showed up if gamma did beat OrthoT. Verify.
       */
      //thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).bwMultiplyDynOrthoT();
      numRecd     = 0;
    }
}




inline void PairCalculator::enqueueBWsend(bool unitcoef, int priority)
{
    /// Create a signal message
    sendBWsignalMsg *sigmsg;
    if(cfg.shouldDelayBWsend)
    {
        sigmsg = new (8*sizeof(int)) sendBWsignalMsg;
        CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
        *(int*)CkPriorityPtr(sigmsg) = priority; // Just make it slower (default value is 1)
    }
    else
        sigmsg = new sendBWsignalMsg;

    /// Collapse this into 1 flag
    if(amPhantom)
        sigmsg->otherdata= true;
    else if(((!cfg.arePhantomsOn && cfg.isSymmetric) || !unitcoef) && notOnDiagonal)
        sigmsg->otherdata=true;
    else
        sigmsg->otherdata= false;

    /// Either reduce the results or sum them direct in GSpace
    if(cfg.isOutputReduced)
        thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);
    else
        thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResultDirect(sigmsg);
}




//PairCalculator::multiplyResult(int size, double *matrix1, double *matrix2)
void PairCalculator::multiplyResultI(multiplyResultMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
    multiplyResult(msg);
}


/**
 * Tolerance correction PsiV Backward path multiplication
 * This is the same as the regular one matrix backward path with the following exceptions:
 * - inDataLeft and inDataRight contain PsiV
 * - outData contains the orthoT from the previous (standard) backward path invocation
 */
void PairCalculator::multiplyPsiV()
{
	#ifdef DEBUG_CP_PAIRCALC_PSIV
		CkPrintf("[%d,%d,%d,%d,%d] In multiplyPsiV\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
	#endif
        
    // If I am a non-phantom chare in the symmetric instance, send a copy of my data to my mirror phantom chare
    if(!amPhantom && cfg.arePhantomsOn && cfg.isSymmetric && notOnDiagonal) 
    {
        CkAssert(existsRight);
        paircalcInputMsg *msg2phantom = new (numExpectedY*numPoints, 8*sizeof(int)) paircalcInputMsg(numPoints,0,false,true,(complex*)inDataRight,true,blkSize,numExpectedY);
        bool prioPhan=false;
        if(prioPhan)
        {
            CkSetQueueing(msg2phantom, CK_QUEUEING_IFIFO);
            *(int*)CkPriorityPtr(msg2phantom) = 1; // just make it slower than non prioritized
        }
        thisProxy(thisIndex.w,thisIndex.y, thisIndex.x,thisIndex.z).acceptRightData(msg2phantom);
    }

    /// We do not need to go through the all of multiplyresult for PsiV. All we really need is the setup for multiplyHelper
    // call helper function to do the math
    int  size=grainSizeX*grainSizeY;
    bool unitcoef=true;
    // TODO: figure out relationship between n_in k_in and grainSizeX grainSizeY
    int m_in=numPoints*2;   // rows of op(A)==rows C
    int n_in=grainSizeY;     // columns of op(B)==columns C
    int k_in=grainSizeX;     // columns op(A) == rows op(B)
    /*
    if(amPhantom)
    {
        n_in=grainSizeX;
        k_in=grainSizeY;
    }
    */
    double beta(0.0);
    int orthoX=0;
    int orthoY=0;
    //BTransform=T offsets for C and A matrices
    int BTCoffset=0;
    int BTAoffset=0;
    //BTransform=N offsets for C and A matrices
    int BNCoffset=0;
    int BNAoffset=0;
    actionType=PSIV;
    bwMultiplyHelper(size, outData, NULL, outData, NULL,  unitcoef, m_in, n_in, k_in, BNAoffset, BNCoffset, BTAoffset, BTCoffset, orthoX, orthoY, beta, grainSizeX, grainSizeY);
    /// Schedule the entry methods that will send the bw results out
    enqueueBWsend(unitcoef);
}




/**
 * Backward path multiplication
 */
void PairCalculator::multiplyResult(multiplyResultMsg *msg)
{
    //============================================================================
    // Do not delete msg. Its a nokeep.
    //============================================================================
    
    #ifdef _PAIRCALC_DEBUG_
        CkPrintf("[%d,%d,%d,%d,%d]: MultiplyResult from orthoX %d orthoY %d size %d numRecdBW %d actionType %d amPhantom %d notOnDiagonal %d cfg.arePhantomsOn %d symmetricOnDiagonal %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, msg->orthoX, msg->orthoY, msg->size, numRecdBW, msg->actionType, amPhantom, notOnDiagonal, cfg.arePhantomsOn, symmetricOnDiagonal);
    #endif
    #ifdef _CP_SUBSTEP_TIMING_
        if(cfg.backwardTimerID>0)
            if(numRecdBW==0)
            {
                double pstart=CmiWallTimer();
                contribute(sizeof(double),&pstart,CkReduction::min_double, cfg.beginTimerCB , cfg.backwardTimerID);
            }
    #endif

    /// Increment the number of tiles received in the backward path
    numRecdBW++;

    #ifdef _NAN_CHECK_
        for(int i=0;i<msg->size;i++)
            CkAssert(finite(msg->matrix1[i]));
    #endif

    int size=msg->size;
    int size2=msg->size2;
    double *matrix1=msg->matrix1;
    double *matrix2=msg->matrix2;

    #ifdef TEST_ALIGN
        if((unsigned int) msg->matrix1 %16 !=0)
        {
            CkPrintf("msg->matrix1 is %p\n",msg->matrix1);
        }
        CkAssert((unsigned int) msg->matrix1 %16 ==0);
        CkAssert((unsigned int) msg->matrix2 %16 ==0);
    #endif

    actionType=msg->actionType;
    bool unitcoef = false;

    int maxorthostateindex=(cfg.numStates/cfg.orthoGrainSize-1)*cfg.orthoGrainSize;
    /// Find our tile indices within this sGrain
    int orthoX=msg->orthoX*cfg.orthoGrainSize;

    int orthoY=msg->orthoY*cfg.orthoGrainSize;
    orthoX= (orthoX>maxorthostateindex) ? maxorthostateindex : orthoX;
    orthoY= (orthoY>maxorthostateindex) ? maxorthostateindex : orthoY;
    if(amPhantom)
    {
        orthoX=(orthoX-thisIndex.y)/cfg.orthoGrainSize;
        orthoY=(orthoY-thisIndex.x)/cfg.orthoGrainSize;
        int swap=orthoY;
        orthoY=orthoX;
        orthoX=swap;
        //      CkPrintf("[%d,%d,%d,%d,%d]: phantom MultiplyResult with size %d numRecdBW %d actionType %d numPoints %d orthoX %d orthoY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, msg->size, numRecdBW, msg->actionType, numPoints, orthoX, orthoY);
    }
    else
    {
        orthoX=(orthoX-thisIndex.x)/cfg.orthoGrainSize;
        orthoY=(orthoY-thisIndex.y)/cfg.orthoGrainSize;
    }
    
    int orthoGrainSizeY=(orthoY==numOrthoRow-1) ? cfg.orthoGrainSize+orthoGrainSizeRemY : cfg.orthoGrainSize;
    int orthoGrainSizeX=(orthoX == numOrthoCol-1) ? cfg.orthoGrainSize + orthoGrainSizeRemX : cfg.orthoGrainSize;
    //  CkPrintf("orthoGrainSizeX %d orthoGrainSizeY %d orthoX %d orthoY %d e1 %.10g\n",orthoGrainSizeX, orthoGrainSizeY, orthoX, orthoY, msg->matrix1[0]);
    //  CkPrintf("orthoGrainSizeX*orthoGrainSizeY = %d msg->size %d\n",orthoGrainSizeY*orthoGrainSizeX, msg->size);
    if(matrix2==NULL||size2<1)
    {
        unitcoef = true;
    }

    /// If I am a phantom chare and have not yet received the right matrix data from my non-phantom mirror
    if(amPhantom && !existsRight)
    { //our mirror data is delayed
        cfg.areBWTilesCollected=true;
        cfg.isBWstreaming=false;
        CkPrintf("[%d,%d,%d,%d,%d] Warning! phantom got bw before fw, forcing tile collection\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
    }

    int matrixSize=grainSizeX*grainSizeY;

    /// @note: ASSUMING TMATRIX IS REAL (LOSS OF GENERALITY)

    double *amatrix=NULL;
    double *amatrix2=matrix2;  ///< @note: may be overridden later

    #ifdef _PAIRCALC_DEBUG_PARANOID_BW_
        CkPrintf("orthoGrainSizeX %d orthoGrainSizeY %d orthoX %d orthoY %d e1 %.10g\n",orthoGrainSizeX, orthoGrainSizeY, orthoX, orthoY, msg->matrix1[0]);
        if(cfg.grainSize==cfg.orthoGrainSize)
        {
            dumpMatrixDouble("bwm1idata",msg->matrix1,grainSizeX,grainSizeY,0,0,orthoX, orthoY);
            if(!unitcoef)
            { // CG non minimization case
                dumpMatrixDouble("bwm2idata",msg->matrix2,grainSizeX,grainSizeY);
            }
        }
        else
        {          
            dumpMatrixDouble("bwm1idata",msg->matrix1,orthoGrainSizeX,orthoGrainSizeY,0,0,orthoX, orthoY);
            if(!unitcoef)
            { // CG non minimization case
                dumpMatrixDouble("bwm2idata",msg->matrix2,orthoGrainSizeX,orthoGrainSizeY,0,0,orthoX, orthoY);
            }
        }
    #endif

    // default DGEMM for non streaming comp case
    int m_in=numPoints*2;   // rows of op(A)==rows C
    // TODO: is n_in grainSizeX or is k_in?
    int n_in=grainSizeY;     // columns of op(B)==columns C
    int k_in=grainSizeX;     // columns op(A) == rows op(B)
    double beta(0.0);

    /// expectOrthoT is true only in asymm, dynamics. Hence only for the asymm, off-diagonal chares in dynamics
    if(expectOrthoT && !notOnDiagonal)
        beta=1.0;  // need to subtract fpsi*orthoT
    // default these to 0, will be set for streaming comp if !collectAllTiles
    //BTransform=T offsets for C and A matrices
    int BTCoffset=0;
    int BTAoffset=0;
    //BTransform=N offsets for C and A matrices
    int BNCoffset=0;
    int BNAoffset=0;
    
    if(cfg.numStates==cfg.grainSize)// all at once no malloc
    {
        amatrix=matrix1;  // index is 0 in this case, so this is silly
    }

    /// If the grain size for paircalc and ortho are the same
    if (cfg.orthoGrainSize==cfg.grainSize)
    { // you were sent the correct section only
        amatrix=matrix1;
        // the other tiles were already collected for PSIV
        numRecdBW=numOrtho;
    }
    /// else, if this is a PsiV loop
    else if(actionType==PSIV)
    {
        amatrix=matrix1;
        // the other tiles were already collected for PSIV
        numRecdBW=numOrtho;
    }
    /// else, if PC is collecting all the tiles from Ortho 
    else if (cfg.areBWTilesCollected)
    {
        collectTile(true, !unitcoef, false,orthoX, orthoY, orthoGrainSizeX, orthoGrainSizeY, numRecdBW, matrixSize, matrix1, matrix2);
        amatrix = inResult1;
        amatrix2 = inResult2;
    } //else !collectAllTiles
    else
    {  // settings for streaming computation on each tile
        /* For Reference to collect tiles we do this
         * if(cfg.isSymmetric && notOnDiagonal) //swap the non diagonals
         * {  // this is trickier to do than one might expect
         * // because the tiles can have funny shapes
         * bigGindex=grainSizeX;
         * bigOindex=orthoGrainSizeX;
         * orthoXoff=cfg.orthoGrainSize*bigGindex*orthoY;
         * orthoYoff=orthoX*cfg.orthoGrainSize;
         * 
         * }
         * // which means we set up our multiply for oY lines of oX size
         * // embedded in gY lines of gX size
         */
        
        amatrix=matrix1;
        if(!unitcoef)
            amatrix2=matrix2;
        // fix n_in and k_in
        n_in=orthoGrainSizeY;
        k_in=orthoGrainSizeX;
        if(cfg.isSymmetric && notOnDiagonal) //swap the non diagonals
        {
            if(!amPhantom)
            {
                int swap=orthoX;
                orthoX=orthoY;
                orthoY=swap;
            }
            /* int swap=orthoX;
             * orthoX=orthoY;
             * orthoY=swap;
             * if(amPhantom)
             * {
             *      n_in=orthoGrainSizeX;
             *      k_in=orthoGrainSizeY;
             * }
             * k_in=(orthoX==numOrthoRow-1) ? cfg.orthoGrainSize+orthoGrainSizeRemX : cfg.orthoGrainSize;
             * n_in=(orthoY == numOrthoCol-1) ? cfg.orthoGrainSize + orthoGrainSizeRemY : cfg.orthoGrainSize;
             */
        }
        
        // skip to the rows which apply to this ortho
        BTCoffset=orthoY * m_in * cfg.orthoGrainSize;
        BTAoffset=orthoX * m_in * cfg.orthoGrainSize;
        BNCoffset=orthoX * m_in * cfg.orthoGrainSize;
        BNAoffset=orthoY * m_in * cfg.orthoGrainSize;

        /*
        if(symmetricOnDiagonal)
        {
            BNCoffset=orthoY * m_in * cfg.orthoGrainSize;
            BNAoffset=orthoX * m_in * cfg.orthoGrainSize;
        }
        */
        beta=1.0;  // need to sum over tiles within orthoY columns
    }
    
    /// If we've received all the tiles, or if we're streaming the computations
    if(cfg.orthoGrainSize==cfg.grainSize || numRecdBW==numOrtho || !cfg.areBWTilesCollected || actionType==PSIV)
    { // have all the input we need
        // call helper function to do the math
        if(actionType!=PSIV && !cfg.areBWTilesCollected && n_in*k_in>size)
        {
            CkPrintf("[%d,%d,%d,%d,%d] Warning! your n_in %d and k_in %d is larger than matrix1->size %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric, n_in,k_in,size);
            n_in=cfg.orthoGrainSize;
            k_in=cfg.orthoGrainSize;
        }
        bwMultiplyHelper(size, matrix1, matrix2, amatrix, amatrix2,  unitcoef, m_in, n_in, k_in, BNAoffset, BNCoffset, BTAoffset, BTCoffset, orthoX, orthoY, beta, orthoGrainSizeX, orthoGrainSizeY);
    }

    //#define SKIP_PARTIAL_SEND
    #ifndef SKIP_PARTIAL_SEND
        /// If we're stream computing without any barriers, simply send whatever results are ready now
        if(cfg.isBWstreaming && !cfg.areBWTilesCollected && !cfg.isBWbarriered && actionType!=PSIV)
        { // send results which are complete and not yet sent
            /*
            if(cfg.isSymmetric && notOnDiagonal &&!amPhantom) //swap the non diagonals
            {
                k_in=(orthoX==numOrthoRow-1) ? cfg.orthoGrainSize+orthoGrainSizeRemX : cfg.orthoGrainSize;
                n_in=(orthoY == numOrthoCol-1) ? cfg.orthoGrainSize + orthoGrainSizeRemY : cfg.orthoGrainSize;
            }
            */
            bwSendHelper( orthoX, orthoY, k_in, n_in, orthoGrainSizeX, orthoGrainSizeY);
            // bwSendHelper( orthoX, orthoY, k_in, n_in, k_in, n_in);

            /// If we're done with all the paircalc work in this loop (iteration), then cleanup
            if(numRecdBW==numOrtho)
                cleanupAfterBWPath();
        }
    #else // dump them all after multiply complete
        /// If SKIP_PARTIAL_SEND is defined, and we're streaming, then we send only if we've received all the tiles
        if((cfg.isBWstreaming && !cfg.areBWTilesCollected && !cfg.isBWbarriered && actionType!=PSIV) && numRecdBW==numOrtho)
            enqueueBWsend(unitcoef);
        else
        {
            CkPrintf("[%d,%d,%d,%d,%d] not sending cfg.isBWstreaming %d cfg.areBWTilesCollected %d cfg.isBWbarriered %d actionType %d  numRecdBW %d numOrtho%d \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,cfg.isBWstreaming,cfg.areBWTilesCollected,cfg.isBWbarriered,actionType, numRecdBW,numOrtho);
        }
    #endif

    /** If all the computations are done, and we're either collecting tiles or are barriered, then we need to send our results to GSpace now
     * @note: We do not enter this if we're streaming without barriers, as in that case data will be sent out as it is ready
     */
    /** @note: cfg.orthoGrainSize == cfg.grainSize implies this PC chare will get only one msg in the bw path. 
     * This is covered by the check numRecdBW==numOrtho, hence is just for paranoia.
     *
     * @note: It appears actionType = PSIV will not enter multiplyResult. Needs to be verified.
     */
    if(   ((!cfg.isBWstreaming && cfg.areBWTilesCollected) && (cfg.orthoGrainSize==cfg.grainSize || numRecdBW==numOrtho))
       || (cfg.isBWbarriered && (cfg.orthoGrainSize==cfg.grainSize || numRecdBW==numOrtho))
       || actionType==PSIV
      )
    { // clean up
        /// If we're barriered on the backward path
        if(cfg.isBWbarriered)
        {
            int wehaveours=1;
            contribute(sizeof(int),&wehaveours,CkReduction::sum_int,CkCallback(CkIndex_PairCalculator::bwbarrier(NULL),thisProxy));
        }
        else
            enqueueBWsend(unitcoef);

        /// If we're not streaming, delete inResult*
        if(cfg.conserveMemory>=0 && (cfg.areBWTilesCollected || !unitcoef))
        {
            if(inResult2!=NULL)
                delete [] inResult2;
            if(inResult1!=NULL)
                delete [] inResult1;
            inResult1=NULL;
            inResult2=NULL;
        }
    }
}



/** Resets counters, zeroes out data structures, frees mem when appropriate. 
 * Readies other components (collators) for the next iteration. Ends the backward path substep
 * timing. Paircalc behavior under different cfg.conserveMemory settings are defined solely in this 
 * method.
 */
void PairCalculator::cleanupAfterBWPath()
{
    #ifdef _PAIRCALC_DEBUG_
        CkPrintf("[%d,%d,%d,%d,%d] Cleaning up at end of BW path\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric);
    #endif

    /// Reset the backward path tile counter
    numRecdBW = 0;

    /// For all mem usage modes
    {
        /// If we're stream computing, zero out these containers
        if (cfg.isBWstreaming)
        {
            bzero(columnCount, sizeof(int) * numOrthoCol);
            bzero(columnCountOther, sizeof(int) * numOrthoCol);
        }

        /// Phantom chares will get a fresh right matrix every time in an incoming msg. Hence, delete it irrespective of cfg.conserveMemory or RDMA settings
        if (amPhantom)
        {
            delete msgRight;
            msgRight    = 0;
            inDataRight = 0;
            existsRight = false;
        }
    }

    /// For only the speed-optimized (high-mem) mode
    if (cfg.conserveMemory == -1)
    {
        if (!amPhantom)
        {
            bzero(mynewData,numPoints*numExpectedY* sizeof(complex));
            if (othernewData) ///< @todo: is this redundant? isnt othernewData always allocated?
                bzero(othernewData,numPoints*numExpectedX * sizeof(complex));
	    }
        else
        {
            if (othernewData) ///< @todo: is this redundant? isnt othernewData always allocated?
                bzero(othernewData,numPoints*numExpectedY * sizeof(complex));
        }
    }


    /// For the normal and strict low-mem modes
    if (cfg.conserveMemory >= 0)
    {
        /// These deletes lived in sendBWResult() and sendBWResultDirect() 
        if (!amPhantom)
        {
            if (mynewData)
                delete [] mynewData;
            mynewData    = 0;
            existsNew    = false;
        }
        if (othernewData)
            delete [] othernewData;
        othernewData = 0;
    }
    

    /// For the strictest of the low-mem footprint modes, 
    if (cfg.conserveMemory == 1)
    {
        /** Delete the input matrices and they'll get reallocated on the next pass. 
         * Only do this if we dont use RDMA because we dont want to setup an RDMA channel every iteration
         */
        #ifndef PC_USE_RDMA
            delete msgLeft;
            msgLeft       = 0;
            inDataLeft    = NULL;
            existsLeft    = false;
            /// If this chare received a right matrix message and it still exists
            if(msgRight && (!cfg.isSymmetric || (cfg.isSymmetric && notOnDiagonal)) )
            {
                delete msgRight;
                msgRight    = 0;
                inDataRight = 0;
                existsRight = false;  
            }
        #endif

        /// Delete the storage for the output matrix too
        if(outData!=NULL && actionType!=KEEPORTHO)
        {
            delete [] outData;
            outData = NULL;
            existsOut=false;
        }
    }


    #ifdef PC_USE_RDMA
        // Unless a PsiV step is next, let the collators know that they should now expect the next batch of data via RDMA. If a PsiV step is next, then PsiV data will come in via traditional messages. expectNext() will be called later at the end of multiplyPsiV().
        if ( !(cfg.isSymmetric && actionType == KEEPORTHO) )
        {
            leftCollator->expectNext();
            rightCollator->expectNext();
        }
        else
        {
            #ifdef DEBUG_CP_PAIRCALC_PSIV
                CkPrintf("[%d,%d,%d,%d,%d]: Am NOT notifying the message handlers to expectNext() as a PsiV step is next (actionType=%d). Data should be arriving in messages. \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,  numRecdBW, actionType);
            #endif
        }
    #endif


    #ifdef _CP_SUBSTEP_TIMING_
        /// End the timing for this sub-step.
        if(cfg.backwardTimerID > 0)
        {
            double pstart=CmiWallTimer();
            contribute(sizeof(double),&pstart,CkReduction::max_double, cfg.endTimerCB , cfg.backwardTimerID);
        }
    #endif
}




void PairCalculator::collectTile(bool doMatrix1, bool doMatrix2, bool doOrthoT, int orthoX, int orthoY, int orthoGrainSizeX, int orthoGrainSizeY, int numRecdBW, int matrixSize, double *matrix1, double* matrix2)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: collectTile aggregating numRecdBW %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numRecdBW);
#endif


  //do strided copy
  // stride amatrix by cfg.grainSize,
  // stride matrix1 by orthoGrainSize
  // copy out the orthograinSize section
  // into amatrix[goffset]
  // goffset+=grainSize
  if(!doOrthoT && doMatrix1 && numRecdBW==1 && inResult1==NULL) //alloc on first receipt
    {
      CkAssert(inResult1==NULL);
      inResult1 = new double[matrixSize];
      bzero(inResult1,sizeof(double)*matrixSize);
    }

  // advance tilestart to new row and column
  int bigGindex=grainSizeY;
  int bigOindex=orthoGrainSizeY;

  int orthoXoff=cfg.orthoGrainSize*bigGindex*orthoX;
  int orthoYoff=orthoY*cfg.orthoGrainSize;

  if(cfg.isSymmetric && notOnDiagonal) //swap the non diagonals
    {
      if(amPhantom)
	{
	  bigGindex=grainSizeY;
	  bigOindex=orthoGrainSizeY;
	  orthoXoff=cfg.orthoGrainSize*bigGindex*orthoX;
	  orthoYoff=orthoY*cfg.orthoGrainSize;
	}
      else{
	bigGindex=grainSizeX;
	bigOindex=orthoGrainSizeX;
	orthoXoff=cfg.orthoGrainSize*bigGindex*orthoY;
	orthoYoff=orthoX*cfg.orthoGrainSize;
      }
    }
  int tileStart=orthoYoff+orthoXoff;


  int tileSize=orthoGrainSizeX*orthoGrainSizeY;

  int ocopySize=bigOindex*sizeof(double);

  int savetileStart=tileStart;

  if(doMatrix2){ // CG non minimization case have GAMMA
    if(inResult2==NULL && numRecdBW==1) //alloc on first receipt
      {
	inResult2 = new double[matrixSize];
	bzero(inResult2,sizeof(double)*matrixSize);
      }
	for(int i=0; i<tileSize; i+=bigOindex,tileStart+=bigGindex)
	  CmiMemcpy(&(inResult2[tileStart]),&(matrix2[i]),ocopySize);

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
    char filename[80];
    snprintf(filename,80,"bwinResult2_%d_%d:",orthoX,orthoY);
    dumpMatrixDouble(filename, matrix2, orthoGrainSizeX, orthoGrainSizeY,orthoX*cfg.orthoGrainSize, orthoY*cfg.orthoGrainSize);
#endif

  }
  double *dest= (doOrthoT) ? outData : inResult1;

  tileStart=savetileStart;

  if(doMatrix1)
    {
#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
  CkPrintf("[%d,%d,%d,%d,%d]: collectTile copying tile bigOindex %d bigGindex %d tileSize %d tileStart %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, bigOindex, bigGindex, tileSize,tileStart);
#endif

      for(int i=0; i< tileSize; i+=bigOindex,tileStart+=bigGindex)
	CmiMemcpy(&(dest[tileStart]),&(matrix1[i]),ocopySize);
#ifdef _PAIRCALC_DEBUG_PARANOID_BW_

      char filename[80];
      if(doOrthoT)
	snprintf(filename,80,"bworthoT_%d_%d:",orthoX,orthoY);
      else
	snprintf(filename,80,"bwinResult1_%d_%d:",orthoX,orthoY);
      dumpMatrixDouble(filename, matrix1, orthoGrainSizeX, orthoGrainSizeY,orthoX*cfg.orthoGrainSize, orthoY*cfg.orthoGrainSize);
#endif
    }


}



void PairCalculator::bwbarrier(CkReductionMsg *msg)
{
      // everyone is done
      delete msg;
      bool unitcoef=true;  //cheap hack for minimzation only
      enqueueBWsend(unitcoef);
}



// PRE: amatrix contains the entire matrix for the multiplication.
//      matrix contains just what was sent in the trigger message.
void PairCalculator::bwMultiplyHelper(int size, double *matrix1, double *matrix2, double *amatrix, double *amatrix2, bool unitcoef, int m_in, int n_in, int k_in, int BNAoffset, int BNCoffset, int BTAoffset, int BTCoffset, int orthoX, int orthoY, double beta, int orthoGrainSizeX, int orthoGrainSizeY)
{


#ifdef _PAIRCALC_DEBUG_

  if(numRecdBW==numOrtho|| cfg.orthoGrainSize==cfg.grainSize || cfg.areBWTilesCollected)
    {
      streamCaughtR++;
    }
  CkPrintf("[%d,%d,%d,%d,%d]: bwMultiplyHelper with size %d numRecdBW %d actionType %d orthoX %d orthoY %d orthoGrainSizeX %d orthoGrainSizeY %d BTCoffset %d BNCoffset %d m_in %d n_in %d k_in %d iter %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, size, numRecdBW, actionType, orthoX, orthoY,orthoGrainSizeX, orthoGrainSizeY, BTCoffset, BNCoffset, m_in, n_in, k_in, streamCaughtR);
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
  if(cfg.orthoGrainSize==cfg.grainSize || cfg.areBWTilesCollected)
    {
      dumpMatrixDouble("bwm1cidata",amatrix,grainSizeX,grainSizeY,0,0,0,streamCaughtR);
      if(!unitcoef)
	{ // CG non minimization case
	  dumpMatrixDouble("bwm2cidata",amatrix2,grainSizeX, grainSizeY,0,0,0,streamCaughtR);
	}
    }
#endif
  int  matrixSize=grainSizeX*cfg.grainSize;
  if(cfg.isSymmetric && actionType==KEEPORTHO) // there will be a psiV step following
    {

      if(outData==NULL)
	{ //reuse outData to hold onto ortho
	  CkAssert(!existsOut);
	  outData=new double[matrixSize];
	  bzero(outData,sizeof(double)*matrixSize);
	  existsOut=true;
	}
      //keep the orthoT we just received in matrix1
      if(!cfg.areBWTilesCollected && cfg.orthoGrainSize!=cfg.grainSize)
	{ // copy this tile
	  collectTile(true, false, true,orthoX, orthoY, orthoGrainSizeX, orthoGrainSizeY, numRecdBW, matrixSize, matrix1, matrix2);
	}
      else
	{
	  CmiMemcpy(outData, amatrix, size*sizeof(double));
	}
      // it is safe to reuse this memory
      // normal backward path has no use for outData
      // forward path won't be called again until after we're done with outData.
    }
  if(!amPhantom && mynewData==NULL)
    {
      CkAssert(numPoints>0);
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] Allocated mynewData %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numExpectedY,numPoints);
#endif

      mynewData = new complex[numPoints*numExpectedY];
      bzero(mynewData,numPoints*numExpectedY* sizeof(complex));
      existsNew=true;
      if(!amPhantom && ((cfg.isSymmetric || !unitcoef) && notOnDiagonal)){
	if(othernewData==NULL)
	  {
#ifdef _PAIRCALC_DEBUG_
	    CkPrintf("[%d,%d,%d,%d,%d] Allocated other %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numExpectedX,numPoints);
#endif
	    othernewData = new complex[numPoints*numExpectedX];
	    bzero(othernewData,numPoints*numExpectedX * sizeof(complex));
	  }
      }
    }
  else if(amPhantom)
    { // phantoms live in a bizarre reverso world
      if(othernewData==NULL)
	{
#ifdef _PAIRCALC_DEBUG_
	    CkPrintf("[%d,%d,%d,%d,%d] Allocated other %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numExpectedY,numPoints);
#endif

	  othernewData = new complex[numPoints*numExpectedY];
	  bzero(othernewData,numPoints*numExpectedY * sizeof(complex));
	}
  }


  double *mynewDatad= reinterpret_cast <double *> (mynewData);
  double alpha(1.0);
  char transform='N';
  char transformT='T';

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_

  int chunksize=blkSize/cfg.numChunks;
  int ystart=chunksize*thisIndex.z;
  if(!amPhantom)
    dumpMatrixDouble("bwmlodata",inDataLeft,numExpectedX,numPoints*2,thisIndex.x,0);
  if(!unitcoef||amPhantom){ // CG non minimization case
    dumpMatrixDouble("bwmrodata",inDataRight,numExpectedY,numPoints*2,thisIndex.y,0);
  }
#endif
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

  //first multiply to apply the T or L matrix
  if(!amPhantom)
    {
      int lk_in=k_in;
      int ln_in=n_in;
      if(symmetricOnDiagonal && orthoX!=orthoY && !cfg.areBWTilesCollected && actionType!=PSIV )
	{
	  // here is where we need some more juggling
	  lk_in=n_in;
	  ln_in=k_in;
	}


#if PC_BWD_DGEMM_SPLIT > 0
      if(cfg.isSymmetric)
	{dgemmSplitBwdM(m_in, ln_in, lk_in, &transform, &transform, &alpha,
			&(inDataLeft[BNAoffset]),  amatrix, &beta,
			&(mynewDatad[BNCoffset]));}
      else
	{dgemmSplitBwdM(m_in, n_in, k_in, &transform, &transformT, &alpha,
			&(inDataLeft[BTAoffset]),  amatrix, &beta,
			&(mynewDatad[BTCoffset]));}
#else // no split

#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(inDataLeft[BNAoffset] )%16==0);
      CkAssert((unsigned int) amatrix %16==0);
      CkAssert((unsigned int)&(mynewDatad[BNCoffset] )%16==0);
#endif

      if(cfg.isSymmetric) {
#ifdef PRINT_DGEMM_PARAMS
	CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d BNAoffset %d BNCoffset %d\n", transform, transform, m_in, n_in, k_in, alpha, beta, m_in, n_in, m_in, BNAoffset, BNCoffset);
#endif
	//note, your choice of k_in or n_in only matters on the
	//border column, elsewhere k_in==n_in
	//orig no valgrind invalid read 8 on A

	DGEMM(&transform, &transform, &m_in, &ln_in, &lk_in, &alpha,
	      &(inDataLeft[BNAoffset]), &m_in, amatrix, &lk_in, &beta,
	      &(mynewDatad[BNCoffset]), &m_in);
      }
      else {
#ifdef PRINT_DGEMM_PARAMS
	CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transformT, m_in, n_in, k_in, alpha, beta, m_in, n_in, m_in);
#endif
	DGEMM(&transform, &transformT, &m_in, &n_in, &k_in, &alpha,
	      &(inDataLeft[BTAoffset]), &m_in,  amatrix, &n_in, &beta,
	      &(mynewDatad[BTCoffset]), &m_in);
      }
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(230, StartTime, CmiWallTimer());
#endif

#endif // end of else clause for split

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
      char snark[80];
      snprintf(snark,80,"bwgmodata_%d_%d:",orthoX,orthoY);
      if(cfg.isSymmetric)
	dumpMatrixComplex(snark,mynewData,numExpectedY,numPoints,0,ystart,streamCaughtR);
      else
	dumpMatrixComplex(snark,mynewData,numExpectedY,numPoints,0,ystart,streamCaughtR);
#endif
    }// end of !amPhantom

#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  // Multiply to compensate for the missing triangle in symmetric case
  if((amPhantom || (!cfg.arePhantomsOn && cfg.isSymmetric && notOnDiagonal)) && existsRight)
    {
      double *othernewDatad= reinterpret_cast <double *> (othernewData);
      if(amPhantom)
	{
	  int swap= n_in;
	  n_in=k_in;
	  k_in=swap;

	}
#if PC_BWD_DGEMM_SPLIT > 0
      dgemmSplitBwdM(m_in, k_in, n_in, &transform, &transformT, &alpha, &(inDataRight[BTAoffset]),  amatrix, &beta, &(othernewDatad[BTCoffset]));
#else  // no split

      CmiNetworkProgress();
#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(inDataRight[BTAoffset] )%16==0);
      CkAssert((unsigned int) amatrix %16==0);
      CkAssert((unsigned int)&(othernewDatad[BTCoffset] )%16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
      CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n",
	       transform, transformT, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif
      DGEMM(&transform, &transformT, &m_in, &k_in, &n_in, &alpha,
	    &(inDataRight[BTAoffset]), &m_in,  amatrix, &k_in, &beta,
	    &(othernewDatad[BTCoffset]), &m_in);

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(250, StartTime, CmiWallTimer());
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
      char snark[80];
      snprintf(snark,80,"bwgomodata_%d_%d:",orthoX,orthoY);
      if(amPhantom)
	dumpMatrixComplex(snark,othernewData,numExpectedY,numPoints,0,ystart,streamCaughtR);
      else
	dumpMatrixComplex(snark,othernewData,numExpectedX,numPoints,0,ystart,streamCaughtR);
#endif

#endif  // end of split
    }



  if(!unitcoef){ // CG non minimization case  GAMMA

    CkAssert(cfg.areBWTilesCollected || cfg.orthoGrainSize==cfg.grainSize);

    // this code is only correct in the non streaming case
    // output modified by subtracting an application of orthoT

    // C = alpha*A*B + beta*C
    // C= -1 * inRight * orthoT + C
    double *othernewDatad;

    alpha=-1.0;  //comes in with a minus sign
    if(notOnDiagonal){
      // setting beta to 0 here means you cannot stream process the
      // application of orthoT by tile, you can only process them once
      // the completed column has been accumulated in inResult2
      beta=0.0; // new contribution off-diagonal
      if(!cfg.areBWTilesCollected && cfg.orthoGrainSize!=cfg.grainSize)
       	beta=1.0; // need to accumulate
      othernewDatad= reinterpret_cast <double *> (othernewData);
    }else{
      beta=1.0; //subtract contribution from existing on diagonal
      othernewDatad=mynewDatad;
    }//endif

    // Funny thing here, this logic works unchanged for remainder case.
    // off diagonals use the usual funny size othernewData
    // diagonals use a MxM newData
    CmiNetworkProgress();

#ifdef PRINT_DGEMM_PARAMS
    CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transform, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif

#if PC_BWD_DGEMM_SPLIT > 0

    dgemmSplitBwdM(m_in, k_in, n_in, &transform, &transform, &alpha,
		   &(inDataRight[BNAoffset]),  amatrix, &beta,
		   &(othernewDatad[BNCoffset]));

#else // no split
#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif
    DGEMM(&transform, &transform, &m_in, &k_in, &n_in, &alpha,
	  &(inDataRight[BNAoffset]), &m_in, amatrix2, &n_in, &beta,
	  &(othernewDatad[BNCoffset]), &m_in);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(240, StartTime, CmiWallTimer());
#endif
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
    if(notOnDiagonal)  // exists right matrix
      if(amPhantom)
	dumpMatrixComplex("bwg2modata",othernewData,numExpectedY,numPoints,0,ystart,streamCaughtR);
      else
	dumpMatrixComplex("bwg2modata",othernewData,numExpectedX,numPoints,0,ystart,streamCaughtR);
#endif

  } // end  CG case

#ifdef _PAIRCALC_VALID_OUT_
  CkPrintf("[PAIRCALC] [%d,%d,%d,%d,%d] backward gemm out %.10g %.10g %.10g %.10g \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,cfg.isSymmetric, mynewDatad[0],mynewDatad[1],mynewData[numPoints*grainSizeX-1].re,mynewData[numPoints*grainSizeX-1].im);
#endif

}



/* when streaming and dynamics and ograin<sgrain
 * You cannot stream process both gamma and orthoT at the same time.
 * Due to the different tranpose the tile ox=1 oy=0 and ox=0 oy=1
 * operates on one axis of the gamma operation and a different axes
 * for the orthoT.  So the gamma output rows you want to use in orthoT
 * may not be ready yet.
 *
 * To avoid this race condition we send orthoT ahead of time (end of
 * Symm S->T) and multiply it with Fpsi at the end of the Asymm
 * forward path.
 *
 * Then when gamma is ready for the backward path we do the gamma
 * multiply, streaming if desired, and simply perform the subtraction
 * in place by using a nonzero beta in the gamma multiply.
 */
void PairCalculator::bwMultiplyDynOrthoT()
{

    /// Ensure that we're the correct paircalc instance
    CkAssert(expectOrthoT);
    /// Ensure that we really have received all the OrthoT pieces
    CkAssert( numRecdBWOT == numOrtho);
    /// Reset the counter
    numRecdBWOT = 0;

    // output modified by subtracting an application of orthoT
    // C = alpha*A*B + beta*C
    // C= -1 * inRight * orthoT + C

  if(othernewData==NULL)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] Allocated other %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numExpectedX,numPoints);
#endif
      othernewData = new complex[numPoints*numExpectedX];
      bzero(othernewData,numPoints*numExpectedX * sizeof(complex));
      existsNew=true;
    }
  if(mynewData==NULL)
    {
      CkAssert(numPoints>0);
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] Allocated mynewData %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,cfg.isSymmetric,numExpectedY,numPoints);
#endif

      mynewData = new complex[numPoints*numExpectedY];
      bzero(mynewData,numPoints*numExpectedY* sizeof(complex));
      existsNew=true;
    }

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    double alpha=-1.0;  //set it up with a minus sign
    double beta=0.0; // new contribution off-diagonal
    double *othernewDatad= reinterpret_cast <double *> (mynewData);
    if(notOnDiagonal)
      othernewDatad= reinterpret_cast <double *> (othernewData);

    CmiNetworkProgress();

#ifdef PRINT_DGEMM_PARAMS
    CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transform, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif
    char transform='N';
    int m_in=numPoints*2;
    int n_in=grainSizeY;     // columns of op(B)==columns C
    int k_in=grainSizeX;     // columns op(A) == rows op(B)

#if PC_BWD_DGEMM_SPLIT > 0

    dgemmSplitBwdM(m_in, k_in, n_in, &transform, &transform, &alpha,
		   inDataRight, inResult2, &beta,
		   othernewDatad);
#else
    DGEMM(&transform, &transform, &m_in, &k_in, &n_in, &alpha, inDataRight,
    	  &m_in, inResult2, &n_in, &beta, othernewDatad, &m_in);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(240, StartTime, CmiWallTimer());
#endif

#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
    int chunksize=blkSize/cfg.numChunks;
    int ystart=chunksize*thisIndex.z;
    if(notOnDiagonal)
      dumpMatrixComplex("bwg2modata",othernewData,numExpectedX,numPoints,0,ystart,streamCaughtR);
#endif

}




void PairCalculator::bwSendHelper(int orthoX, int orthoY, int sizeX, int sizeY, int orthoGrainSizeX, int orthoGrainSizeY)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper with  numRecdBW %d actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,  numRecdBW, actionType);
#endif


  // Check to see if this column is  complete, if so send it.

  // symmetric uses BNCoffset on newData for diagonal.
  // symmetric also BTCoffset on otherNewData for off diagonal.
  // asymmetric uses BTC offset on newData

  // determine the orthograinsizeXY for the orthoX,orthoY tile
  // do the calculation based on numExpectedX numExpectedY
  // to account for last row col remainder expansion in PC


  int maxorthostateindex=(cfg.numStates/cfg.orthoGrainSize-1)*cfg.orthoGrainSize;

  int orthoXgrain=orthoX*cfg.orthoGrainSize;
  int orthoYgrain=orthoY*cfg.orthoGrainSize;

  orthoXgrain= (orthoXgrain>maxorthostateindex) ? maxorthostateindex : orthoXgrain;
  orthoYgrain= (orthoYgrain>maxorthostateindex) ? maxorthostateindex : orthoYgrain;

  //  CkAssert(orthoGrainSizeY==sizeY);
  //  CkAssert(orthoGrainSizeX==sizeX);
  if(cfg.isSymmetric)
    {
      if(!amPhantom)
	{
	  /*	  if(notOnDiagonal)
	    {
	      int size=sizeX;
	      sizeX=sizeY;
	      sizeY=size;
	      }*/
	  int index=orthoY;
	  if(notOnDiagonal)
	    index=orthoX;
	  if(++columnCount[index]==numOrthoCol ) // BNC
	    {
#ifdef _PAIRCALC_DEBUG_
	      CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper !amPhantom orthoXgrain %d sizeX %d orthoYgrain %d sizeY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,  orthoXgrain, sizeX, orthoYgrain, sizeY);
#endif

	      int   startGrain=orthoYgrain;
	      int   endGrain=startGrain+sizeY;
	      if(notOnDiagonal)
		{
		  startGrain=orthoXgrain;
		  endGrain=startGrain+sizeY;
		}
	      // send orthoX in newData
	      if(cfg.isOutputReduced)
              sendBWResultColumn(false, startGrain, endGrain);
          else
              sendBWResultColumnDirect(false, startGrain, endGrain);
	    }
	  CkAssert(columnCount[index]<=numOrthoCol);
	}

    }
  else //asymm
    {
      if(++columnCount[orthoY]==numOrthoCol) // BTC
	{
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper asymm orthoYgrain %d sizeY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,  orthoYgrain, sizeY);
#endif

	  //send orthoY in newdata
	  if(cfg.isOutputReduced)
          sendBWResultColumn(false, orthoYgrain, orthoYgrain+sizeY);
      else
          sendBWResultColumnDirect(false, orthoYgrain, orthoYgrain+sizeY);
	}
      CkAssert(columnCount[orthoY]<=numOrthoCol);
    }
  // send the othernewData for symm off diag (including phantom) and
  // asymm off diag dynamics
  if((amPhantom || (!cfg.arePhantomsOn && (othernewData!=NULL)&& notOnDiagonal))&& existsRight)
    {
      int index =orthoX;
      if(cfg.isSymmetric)
      	index=orthoY;
      if(++columnCountOther[index]==numOrthoCol) // BTC
	{
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper amPhantom | other orthoXgrain %d sizeX %d orthoYgrain %d sizeY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,  orthoXgrain, sizeX, orthoYgrain, sizeY);
#endif


	    int startGrain=orthoXgrain;
	    int endGrain=sizeX+startGrain;

	    //int startGrain=orthoYgrain;
	    //int endGrain=sizeY+startGrain;
	    if(cfg.isSymmetric)
	      {
		startGrain=orthoYgrain;
		 endGrain=sizeX+startGrain;
		}
	    if(amPhantom)
	      {
		startGrain=orthoYgrain;
		//startGrain=orthoYgrain;
		//		endGrain=startGrain+orthoGrainSizeY;
		//				endGrain=startGrain+sizeY;
		endGrain=startGrain+sizeY;
	      }
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper amPhantom | other startGrain %d endGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,  startGrain, endGrain );
#endif
	  //send orthoY in otherNewData
	  if(cfg.isOutputReduced)
          sendBWResultColumn(true, startGrain, endGrain);
      else
          sendBWResultColumnDirect(true, startGrain, endGrain);
	}
      CkAssert(columnCountOther[index]<=numOrthoCol);
    }

  // this could be refined to track an array of completed columns
  // and send them in some grouping scheme
  CkAssert(numRecdBW<=numOrtho);
}




void
PairCalculator::sendBWResultColumnDirect(bool otherdata, int startGrain, int endGrain  )
{
#ifdef _PAIRCALC_DEBUG_CONTRIB_
  CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumnDirect with otherdata %d actionType %d startGrain %d endGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, otherdata, actionType, startGrain, endGrain);
#endif


  int cp_entry=cb_ep;
  if(actionType==PSIV)
    {
      cp_entry= cb_ep_tol;
    }
  CkAssert(endGrain<=cfg.numStates);
  int numToSend=endGrain-startGrain;
  int permuter=(int) ((float) thisIndex.z/ (float) cfg.numChunks) * (float) numToSend;
  if(otherdata){  // we have this othernewdata issue for the symmetric case
    // and the asymmetric dynamic case
    // for the off diagonal elements
    CkAssert(othernewData!=NULL);
    int index=thisIndex.x;
    if(amPhantom)
      index=thisIndex.y;

    for(int j=startGrain;j<endGrain;j++)
      {
	// shift the send order proportional to chunk index
	int jPermuted=(j+permuter>endGrain) ? (j+permuter-numToSend) : (j+permuter) ;

	complex *computed=&(othernewData[jPermuted*numPoints]);
	CkCallback mycb(cp_entry, CkArrayIndex2D(jPermuted+ index ,thisIndex.w), cb_aid);
	partialResultMsg *msg=new (numPoints, 8*sizeof(int) ) partialResultMsg;
	if(cfg.resultMsgPriority)
	  {
	    *((int*)CkPriorityPtr(msg)) = cfg.resultMsgPriority;
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  }
	msg->init(thisIndex,numPoints, thisIndex.z, computed);

#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d]: sending partial other of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numPoints,jPermuted,index+jPermuted,thisIndex.w);
#endif

#ifdef _NAN_CHECK_
	for(int i=0;i<msg->N ;i++)
	  {
	    if(!finite(msg->result[i].re)|| !finite(msg->result[i].im))
	      {
		CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, i);
	      }
	    CkAssert(finite(msg->result[i].re));
	    CkAssert(finite(msg->result[i].im));
	  }
#endif
	mycb.send(msg);
      }
  }
  else
    {
      CkAssert(!amPhantom);
      CkAssert(mynewData!=NULL);
      for(int j=startGrain;j<endGrain;j++) //mynewdata
	{
	int jPermuted=(j+permuter>endGrain) ? (j+permuter-numToSend) : (j+permuter) ;
	complex *computed=&(mynewData[jPermuted*numPoints]);
#ifdef _NAN_CHECK_
	for(int i=0;i<numPoints ;i++)
	  {
	    if(!finite(computed[i].re) || !finite(computed[i].im))
	      {
		fprintf(stderr,"[%d,%d,%d,%d,%d]: sendBWResultColumnDirect nan in computed at %d %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, jPermuted,i);
		CkAssert(finite(computed[i].re));
		CkAssert(finite(computed[i].im));
	      }
	  }
#endif
	int index=thisIndex.y;
	CkCallback mycb(cp_entry, CkArrayIndex2D(jPermuted+index ,thisIndex.w), cb_aid);
	partialResultMsg *msg=new (numPoints, 8*sizeof(int) ) partialResultMsg;
	if(cfg.resultMsgPriority)
	  {
	    *((int*)CkPriorityPtr(msg)) = cfg.resultMsgPriority;
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  }
	msg->init(thisIndex,numPoints, thisIndex.z, computed);
	/*
	  msg->N=numPoints;
	  msg->myoffset = thisIndex.z; // chunkth
	  CmiMemcpy(msg->result,mynewData+jPermuted*numPoints,msg->N*sizeof(complex));
	*/
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d]:sending partial of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numPoints,jPermuted,index+jPermuted,thisIndex.w);
#endif

#ifdef _NAN_CHECK_
	for(int i=0;i<msg->N ;i++)
	  {
	    if(!finite(msg->result[i].re)|| !finite(msg->result[i].im))
	      {
		CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, i);
	      }
	    CkAssert(finite(msg->result[i].re));
	    CkAssert(finite(msg->result[i].im));
	  }
#endif
	mycb.send(msg);
	}
    }
}



/**
 * Send the result for this column
 */

void PairCalculator::sendBWResultColumn(bool otherdata, int startGrain, int endGrain  )
{


#ifdef _PAIRCALC_DEBUG_CONTRIB_
  CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumn with actionType %d startGrain %d sendGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, actionType, startGrain, endGrain);
#endif

  // Now we have results in mynewData and if(otherData) othernewData
  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  int cp_entry=cb_ep;
  if(actionType==PSIV)
    {
      cp_entry= cb_ep_tol;
    }
  if(otherdata){  // we have this othernewdata issue for the symmetric case
    // and the asymmetric dynamic case
    // for the off diagonal elements
    CkAssert(othernewData!=NULL);
    int index=thisIndex.x;
    if(amPhantom)
      index=thisIndex.y;
    for(int j=startGrain;j<endGrain;j++)
      {
	//this callback creation could be obviated by keeping an
	//array of callbacks, not clearly worth doing
	complex *computed=&(othernewData[j*numPoints]);
#ifndef CMK_OPTIMIZE
	double StartTime=CmiWallTimer();
#endif

	CkCallback mycb(cp_entry, CkArrayIndex2D(j+index ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	CkPrintf("[%d,%d,%d,%d,%d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
#endif

	int outOffset=thisIndex.z;
	mcastGrp->contribute(numPoints*sizeof(complex), computed, sumMatrixDoubleType, otherResultCookies[j], mycb, outOffset);

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
	//	if((j-startGrain) % 8)
	  CmiNetworkProgress();

      }
  }
  else
    {
      CkAssert(mynewData!=NULL);
      int index=thisIndex.y;
      for(int j=startGrain;j<endGrain;j++) //mynewdata
	{

	complex *computed=&(mynewData[j*numPoints]);
#ifndef CMK_OPTIMIZE
	double StartTime=CmiWallTimer();
#endif
	CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);

#ifdef _PAIRCALC_DEBUG_CONTRIB_
	  CkPrintf("[%d,%d,%d,%d,%d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
#endif

	  int outOffset=thisIndex.z;
	  mcastGrp->contribute(numPoints*sizeof(complex), computed, sumMatrixDoubleType, resultCookies[j], mycb, outOffset);

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif

	  //	  if((j-startGrain) % 8)
	    CmiNetworkProgress();
	}
    }
}



void PairCalculator::sendBWResultDirect(sendBWsignalMsg *msg)
{
	#ifdef _PAIRCALC_DEBUG_
		CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultDirect with actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, actionType);
	#endif
	// Now we have results in mynewData and if(cfg.isSymmetric) othernewData
	bool otherdata=msg->otherdata;
	delete msg;
	/// Determine the callback to use
	int cp_entry=cb_ep;
	if(actionType==PSIV)
		cp_entry= cb_ep_tol;

	/// Off-diagonal PCs in the symm (and asymm, in dynamics) instance, have to handle othernewdata
	if(otherdata)
	{
		CkAssert(othernewData!=NULL);
		int size=grainSizeX;
		int index=thisIndex.x;
		if(amPhantom)
		{
			index=thisIndex.y;
			size=grainSizeY;
		}
		for(int j=0;j<size;j++)
		{
			complex *computed=&(othernewData[j*numPoints]);
			//this callback creation could be obviated by keeping an
			//array of callbacks, not clearly worth doing
			CkCallback mycb(cp_entry, CkArrayIndex2D(j+ index ,thisIndex.w), cb_aid);
			partialResultMsg *omsg= new (numPoints, 8*sizeof(int) ) partialResultMsg;
			if(cfg.resultMsgPriority)
			{
				*((int*)CkPriorityPtr(omsg)) = cfg.resultMsgPriority;
				CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
			}
			omsg->init(thisIndex,numPoints, thisIndex.z, computed);
			#ifdef _PAIRCALC_DEBUG_
				CkPrintf("[%d,%d,%d,%d,%d]:sending other partial of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numPoints,j,index+j,thisIndex.w);
			#endif
			#ifdef _NAN_CHECK_
				for(int i=0;i<omsg->N ;i++)
				{
					if(!finite(omsg->result[i].re) || !finite(omsg->result[i].im))
						CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, i);
					CkAssert(finite(omsg->result[i].re));
					CkAssert(finite(omsg->result[i].im));
				}
			#endif
			mycb.send(omsg);
		}
	}
	
	///
	if(!amPhantom)
	{
		CkAssert(mynewData!=NULL);
		int outsize=grainSizeY;
		int index=thisIndex.y;
		for(int j=0;j<outsize;j++) //mynewdata
		{
			complex *computed=&(mynewData[j*numPoints]);
			CkCallback mycb(cp_entry, CkArrayIndex2D(j+index ,thisIndex.w), cb_aid);
			partialResultMsg *omsg= new (numPoints, 8*sizeof(int) ) partialResultMsg;
			if(cfg.resultMsgPriority)
			{
				*((int*)CkPriorityPtr(omsg)) = cfg.resultMsgPriority;
				CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
			}
			omsg->init(thisIndex,numPoints, thisIndex.z, computed);
			#ifdef _PAIRCALC_DEBUG_
				CkPrintf("[%d,%d,%d,%d,%d]:sending partial of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,numPoints,j,index+j,thisIndex.w);
			#endif
			#ifdef _NAN_CHECK_
				for(int i=0;i<omsg->N ;i++)
				{
					if(!finite(omsg->result[i].re) || !finite(omsg->result[i].im))
						CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, i);
					CkAssert(finite(omsg->result[i].re));
					CkAssert(finite(omsg->result[i].im));
				}
			#endif
			mycb.send(omsg);
		}
	}

    /// If we're done with all the paircalc work in this loop (iteration), then cleanup
    if (numRecdBW == numOrtho || actionType == PSIV)
        cleanupAfterBWPath();
}




// entry method to allow us to delay this outbound communication
// to minimize brain dead BG/L interference we have a signal to prioritize this
void
PairCalculator::sendBWResult(sendBWsignalMsg *msg)
{

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: sendBWResult with actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, actionType);
#endif
  bool otherdata=msg->otherdata;
  delete msg;
  // Now we have results in mynewData and if(cfg.isSymmetric) othernewData
  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  int cp_entry=cb_ep;
  if(actionType==PSIV)
    {
      cp_entry= cb_ep_tol;
    }
  /*
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
  */
  if(otherdata){  // we have this othernewdata issue for the symmetric case
    // and the asymmetric dynamic case
    // for the off diagonal elements
    CkAssert(othernewData!=NULL);
    int size=grainSizeX;
    int index=thisIndex.x;
    if(amPhantom)
      {
	index=thisIndex.y;
	size=grainSizeY;
      }
    for(int j=0;j<size;j++)
      {
	//this callback creation could be obviated by keeping an
	//array of callbacks, not clearly worth doing
	complex *computed=&(othernewData[j*numPoints]);
	CkCallback mycb(cp_entry, CkArrayIndex2D(j+index ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	CkPrintf("[%d,%d,%d,%d,%d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
#endif
	/*
#ifndef CMK_OPTIMIZE
	StartTime=CmiWallTimer();
#endif
	*/
	int outOffset=thisIndex.z;
	mcastGrp->contribute(numPoints*sizeof(complex),computed, sumMatrixDoubleType, otherResultCookies[j], mycb, outOffset);
	/*
#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
	*/

      }
  }
  if(!amPhantom)
    {
      int outsize=grainSizeY;
      int index=thisIndex.y;

      CkAssert(mynewData!=NULL);
      for(int j=0;j<outsize;j++) //mynewdata
	{
	  complex *computed=&(mynewData[j*numPoints]);
	  CkCallback mycb(cp_entry, CkArrayIndex2D(j+index,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	  CkPrintf("[%d,%d,%d,%d,%d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, cfg.isSymmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
#endif
	  /*
	    #ifndef CMK_OPTIMIZE
	    StartTime=CmiWallTimer();
	    #endif
	  */
	  int outOffset=thisIndex.z;
	  mcastGrp->contribute(numPoints*sizeof(complex), computed, sumMatrixDoubleType, resultCookies[j], mycb, outOffset);
	  /*
	    #ifndef CMK_OPTIMIZE
	    traceUserBracketEvent(220, StartTime, CmiWallTimer());
	    #endif
	  */
	}
    }

    /// If we're done with all the paircalc work in this loop (iteration), then cleanup
    if (numRecdBW == numOrtho || actionType == PSIV)
        cleanupAfterBWPath();
}



void PairCalculator::dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim, int xstart, int ystart, int xtra1, int xtra2)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename, fmt, thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,
	  xtra1, xtra2, cfg.isSymmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",i+xstart,j+ystart,matrix[i*ydim+j]);
  fclose(loutfile);
}

void PairCalculator::dumpMatrixComplex(const char *infilename, complex *matrix, int xdim, int ydim, int xstart, int ystart, int iter)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename, fmt, thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, iter ,cfg.isSymmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g %.12g \n",i+xstart,j+ystart,matrix[i*ydim+j].re, matrix[i*ydim+j].im);
  fclose(loutfile);
}

void PairCalculator::copyIntoTiles(double *source, double**dest, int sourceRows, int sourceCols, int *offsetsRow, int *offsetsCol, int *touched, int tileSize, int tilesPerRow )

{

  // foreach element
  for(int j=0;j<sourceRows;j++) // x is inner index
    for(int i=0;i<sourceCols;i++)
      {
	int x=offsetsRow[i];
	int y=offsetsCol[j];
	int sourceOffset=j*sourceCols+i;
	double value=source[sourceOffset];
	int tilex=x/tileSize;
	int tiley=y/tileSize*tilesPerRow;
	int tilei=x%tileSize;
	int tilej=y%tileSize;
	touched[tiley+tilex]++;
#ifdef _NAN_CHECK_
	CkAssert(finite(value));
	CkAssert(touched[tiley+tilex]<=tileSize*tileSize);
	CkAssert(tilex+tiley<numOrtho);
#endif
	dest[tiley+tilex][tilej*tileSize+tilei]=value;
	//	if(cfg.isSymmetric)
	//	  CkPrintf(" j %d i %d, x %d y %d copy %g into tilex %d tiley %d, offset %d touched %d\n", j,i, x, y, value, tilex, tiley, tilej*tileSize+tilei, touched[tiley+tilex]);

      }
}

void PairCalculator::dgemmSplitFwdStreamMK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *C, int *ldc)
{
#ifndef CMK_OPTIMIZE
        double StartTime=CmiWallTimer();
#endif
        double betap = 1.0;
        int Ksplit_m =  gemmSplitFWk;
        int Ksplit   = ( (k > Ksplit_m) ? Ksplit_m : k);
        int Krem     = (k % Ksplit);
        int Kloop    = k/Ksplit-1;
        int Msplit_m = gemmSplitFWm;
        int Msplit   = ( (m > Msplit_m) ? Msplit_m : m);
        int Mrem     = (m % Msplit);
        int Mloop    = m/Msplit;

        for(int ms=1;ms<=Mloop;ms++)
        {
          int moff    = (ms-1)*k*Msplit;
          int moffc   = (ms-1)*Msplit;
          int MsplitU = (ms==Mloop ? Msplit+Mrem : Msplit);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
          StartTime=CmiWallTimer();
#endif
#endif
#ifdef TEST_ALIGN
          CkAssert((unsigned int) A[moff] % 16==0);
          CkAssert((unsigned int) B % 16==0);
          CkAssert((unsigned int) C[moffc] % 16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", *transT, *trans, MsplitU, n, Ksplit, *alpha, betap, *lda, *ldb, *ldc);
#endif
          DGEMM(transT, trans, &MsplitU, &n, &Ksplit, alpha, &A[moff], lda, B, ldb, &betap, &C[moffc], ldc);

#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
          traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif
          CmiNetworkProgress();
          for(int ks=1;ks<=Kloop;ks++)
          {
            int koff    = ks*Ksplit;
            int KsplitU = (ks==Kloop ? Ksplit+Krem : Ksplit);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
            StartTime=CmiWallTimer();
#endif
#endif
#ifdef TEST_ALIGN
            CkAssert((unsigned int) A[koff+moff] % 16==0);
            CkAssert((unsigned int) B[koff] % 16==0);
            CkAssert((unsigned int) C[moffc] % 16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", *transT, *trans, MsplitU, n, KsplitU, *alpha, betap, *lda, *ldb, *ldc);
#endif
            DGEMM(transT, trans, &MsplitU, &n, &KsplitU, alpha, &A[koff+moff], lda, &B[koff], ldb, &betap, &C[moffc], ldc);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE

            traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif
          CmiNetworkProgress();
          }//endfor
        }//endfor

#ifdef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
        traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif
}

void PairCalculator::dgemmSplitFwdStreamNK(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, int *lda, double *B, int *ldb, double *C, int *ldc)
{
#ifndef CMK_OPTIMIZE
        double StartTime=CmiWallTimer();
#endif
        double betap = 1.0;
        int Ksplit_m =  gemmSplitFWk;
        int Ksplit   = ( (k > Ksplit_m) ? Ksplit_m : k);
        int Krem     = (k % Ksplit);
        int Kloop    = k/Ksplit-1;
        int Nsplit_m = gemmSplitFWm;
        int Nsplit   = ( (n > Nsplit_m) ? Nsplit_m : n);
        int Nrem     = (n % Nsplit);
        int Nloop    = n/Nsplit;
#ifndef CMK_OPTIMIZE
            StartTime=CmiWallTimer();
#endif
        for(int ns=1;ns<=Nloop;ns++)
        {
          int noff    = (ns-1)*k*Nsplit;
          int noffc   = (ns-1)*(*ldc)*Nsplit;
          int NsplitU = (ns==Nloop ? Nsplit+Nrem : Nsplit);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
          StartTime=CmiWallTimer();
#endif
#endif
#ifdef TEST_ALIGN
          CkAssert((unsigned int) A % 16==0);
          CkAssert((unsigned int) B[noff] % 16==0);
          CkAssert((unsigned int) C[noffc] % 16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", *transT, *trans, m, NsplitU, Ksplit, *alpha, betap, *lda, *ldb, *ldc);
#endif
            DGEMM(transT, trans, &m, &NsplitU, &Ksplit, alpha, A, lda, &B[noff], ldb, &betap, &C[noffc], ldc);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
            traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif
            CmiNetworkProgress();
            for(int ks=1;ks<=Kloop;ks++){
              int koff    = ks*Ksplit;
              int KsplitU = (ks==Kloop ? Ksplit+Krem : Ksplit);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
              StartTime=CmiWallTimer();
#endif
#endif
#ifdef TEST_ALIGN
              CkAssert((unsigned int) A[koff] % 16==0);
              CkAssert((unsigned int) B[koff+noff] % 16==0);
              CkAssert((unsigned int) C[noffc] % 16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", *transT, *trans, m, NsplitU, KsplitU, *alpha, betap, *lda, *ldb, *ldc);
#endif
              DGEMM(transT, trans, &m, &NsplitU, &KsplitU, alpha, &A[koff], lda, &B[koff+noff], ldb, &betap, &C[noffc], ldc);

#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
              traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif
              CmiNetworkProgress();
            }//endfor
          }//endfor

#ifdef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
          traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif
}

void PairCalculator::dgemmSplitBwdM(int m, int n, int k, char *trans, char *transT, double *alpha, double *A, double *B, double *bt, double *C)
{
#ifndef CMK_OPTIMIZE
        double StartTime=CmiWallTimer();
#endif
        int Msplit_m = gemmSplitBW;
        int Msplit   = ( (m > Msplit_m) ? Msplit_m : m);
        int Mrem     = (m % Msplit);
        int Mloop    = m/Msplit-1;

#ifdef TEST_ALIGN
        CkAssert((unsigned int) A % 16==0);
        CkAssert((unsigned int) B % 16==0);
        // CkAssert((unsigned int) C % 16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", *trans, *transT, Msplit, n, k, *alpha, *bt, m, k, m);
#endif
  int ldb=k;
  if(*transT=='T')
    ldb=n;
        DGEMM(trans, transT, &Msplit, &n, &k, alpha, A, &m, B, &ldb, bt, C, &m);
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
        traceUserBracketEvent(230, StartTime, CmiWallTimer());
#endif
#endif

        CmiNetworkProgress();
        for(int i=1;i<=Mloop;i++)
        {
          int off = i*Msplit;
          if(i==Mloop) { Msplit+=Mrem; }
#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
        StartTime=CmiWallTimer();
#endif
#endif
#ifdef TEST_ALIGN
        CkAssert((unsigned int) A[off] % 16==0);
        CkAssert((unsigned int) B % 16==0);
        // CkAssert((unsigned int) C[off] % 16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", *trans, *transT, Msplit, n, k, *alpha, *bt, m, k, m);
#endif
        DGEMM(trans, transT, &Msplit, &n, &k, alpha, &A[off], &m, B, &ldb, bt, &C[off], &m);

#ifndef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
        traceUserBracketEvent(230, StartTime, CmiWallTimer());
#endif
#endif
        CmiNetworkProgress();
      } //endfor
#ifdef BUNDLE_USER_EVENT
#ifndef CMK_OPTIMIZE
        traceUserBracketEvent(230, StartTime, CmiWallTimer());
#endif
#endif
}

void manmult(int numRowsA, int numRowsB, int rowLength, double *A, double *B, double *C, double alpha)
{
  // foreach row multiply all members
  for( int i=0; i<numRowsA; i++)
    {
      double *athrow=&(A[i*rowLength]);
      for( int j=0; j<numRowsB; j++)
	{
	  double *bthrow=&(B[j*rowLength]);
	  for( int e=0; e<rowLength; e++)
	    {
	      //	      C[i*numRowsB +j]+=alpha*athrow[e]*bthrow[e];
	      C[j*numRowsA +i]+=alpha*athrow[e]*bthrow[e];
	    }
	}
    }
}


#include "pcMessages.def.h"
#include "InputDataHandler.h"
#include "inputDataHandler.def.h"
#include "ckPairCalculator.def.h"

