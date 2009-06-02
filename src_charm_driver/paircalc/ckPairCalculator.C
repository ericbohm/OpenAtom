//**************************************************************************
 /** \file ckPairCalculator.C
 * This is a matrix multiply library with extra frills to communicate the
 * results back to gspace or the calling ortho char as directed by the
 * callback.
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
 * NOTE: The magic number 2 appears in 2 contexts.  Either we have
 * twice as many inputs in the non diagonal elements of the matrix.
 * Or in the case of transforming our arrays of complex into arrays of
 * doubles.  This transformation is done so we can use DGEMM instead
 * of ZGEMM.  Motivated by the empirical discovery that BLAS
 * implementors do all their work on DGEMM and leave ZGEMM out in the
 * unoptimized cold.  The latter issue crops up everywhere we have to
 * do allocation or manipulations of the input.  Could arguably be
 * abstracted away by absorbing it into a size field.
*/
//**************************************************************************


#include "ckPairCalculator.h"
#include "pairCalculator.h"
#include <sstream> 

#ifdef CMK_BLUEGENEL
//#include "builtins.h"
// void __alignx(int alignment,  const void *address);
// void _alignx(int alignment,  const void *address);
// void alignx(int alignment,  const void *address);
#endif
//#define PRINT_DGEMM_PARAMS

ComlibInstanceHandle mcastInstanceCP;
ComlibInstanceHandle mcastInstanceACP;

extern CkVec <PairCalcID> UpairCalcID1;
extern CkVec <PairCalcID> UpairCalcID2;

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


PairCalculator::PairCalculator(CkMigrateMessage *m) { }

PairCalculator::PairCalculator(CProxy_InputDataHandler<CollatorType,CollatorType> inProxy, bool _sym, int _grainSize, int _s, int _numChunks, CkCallback _cb, CkArrayID _cb_aid, int _cb_ep, int _cb_ep_tol, int _conserveMemory, bool _lbpaircalc,  redtypes _cpreduce, int _orthoGrainSize, bool _collectTiles, bool _PCstreamBWout, bool _PCdelayBWSend, bool _gSpaceSum, int _gpriority, bool _phantomSym, bool _useBWBarrier, int _gemmSplitFWk, int _gemmSplitFWm, int _gemmSplitBW, bool _expectOrthoT, int _instance)
{
#ifdef _PAIRCALC_DEBUG_PLACE_
  CkPrintf("{%d} [PAIRCALC] [%d,%d,%d,%d,%d] inited on pe %d \n", _instance,thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,_sym, CkMyPe());
#endif


  this->conserveMemory=_conserveMemory;
  this->symmetric = _sym;
  this->grainSize = _grainSize;
  this->numStates = _s;
  instance=_instance;
  int remainder=numStates%grainSize;
  grainSizeX=(numStates- thisIndex.x == grainSize+remainder) ? grainSize+remainder: grainSize;
  grainSizeY=(numStates- thisIndex.y == grainSize+remainder) ? grainSize+remainder: grainSize;

  this->numChunks = _numChunks;
  this->numPoints = -1;
  this->cb_aid = _cb_aid;
  this->cb_ep = _cb_ep;
  this->cb_ep_tol = _cb_ep_tol;
  useBWBarrier=_useBWBarrier;
  PCstreamBWout=_PCstreamBWout;
  PCdelayBWSend=_PCdelayBWSend;
  orthoGrainSize=_orthoGrainSize;
  orthoGrainSizeRemX=grainSizeX%orthoGrainSize;
  orthoGrainSizeRemY=grainSizeY%orthoGrainSize;
  collectAllTiles=_collectTiles;
  gSpaceSum=_gSpaceSum;
  gpriority=_gpriority;
  phantomSym=_phantomSym;
  gemmSplitFWk=_gemmSplitFWk;
  gemmSplitFWm=_gemmSplitFWm;
  gemmSplitBW=_gemmSplitBW;
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
  lbpaircalc=_lbpaircalc;
  expectOrthoT=_expectOrthoT;
  amPhantom=(phantomSym && (thisIndex.y<thisIndex.x) && symmetric ) ? true : false;
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
  symmetricOnDiagonal=(symmetric && thisIndex.x==thisIndex.y) ? true: false;
  // If I lie on the chare array diagonal, I expect to get only left matrix data
  if(symmetricOnDiagonal)
    numExpected=numExpectedX;
  // If I am a phantom chare, I expect to get only right matrix data
  if(amPhantom)
      numExpected = numExpectedY;
  cpreduce=_cpreduce;
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
  if(lbpaircalc)
      setMigratable(true);
  else
      setMigratable(false);
  resultCookies=NULL;
  otherResultCookies=NULL;
  // TODO: technically we can make fewer of these if grainSize>grainSizeX || grainSizeY but you'll only save a few bytes
  resultCookies=new CkSectionInfo[numExpectedX];
  numOrthoCol=grainSizeX/orthoGrainSize;
  numOrthoRow=grainSizeY/orthoGrainSize;
  numOrtho=numOrthoCol*numOrthoRow;
  orthoCookies=new CkSectionInfo[numOrtho*2];
  orthoCB=new CkCallback[numOrtho*2];
  if(PCstreamBWout && !collectAllTiles)
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
	idStream<<"["<<thisIndex.w<<","<<thisIndex.x<<","<<thisIndex.y<<","<<thisIndex.z<<","<<symmetric<<"]";
	/// Create the message handlers for the left and right input matrix blocks
	leftCollator = new CollatorType (idStream.str()+" LeftHandler" , leftTrigger, numExpectedX,(conserveMemory<=0),thisIndex.x);
	rightCollator= new CollatorType (idStream.str()+" RightHandler",rightTrigger, numExpectedY,(conserveMemory<=0),thisIndex.y);
	#ifdef DEBUG_CP_PAIRCALC_INPUTDATAHANDLER
		CkPrintf("[%d,%d,%d,%d,%d] My left and right data collators: %p %p\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,leftCollator,rightCollator);
	#endif
	/// This is the first point during execution when I can supply my InputDataHandler with pointers to the msg handlers, hence
	/// it is (now) safe to insert the [w,x,y,z]th element of the InputDataHandler chare array (as it will immediately clamor 
	/// for access to these handlers)
	myMsgHandler = inProxy;
	#ifdef DEBUG_CP_PAIRCALC_CREATION
		CkPrintf("[%d,%d,%d,%d,%d] Inserting my InputDataHandler\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
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
  p|grainSize;
  p|grainSizeX;
  p|grainSizeY;
  p|numStates;
  p|numChunks;
  p|numPoints;
  p|symmetric;
  p|symmetricOnDiagonal;
  p|notOnDiagonal;
  p|conserveMemory;
  p|lbpaircalc;
  p|cpreduce;
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
  p|orthoGrainSize;
  p|orthoGrainSizeRemX;
  p|orthoGrainSizeRemY;
  p|PCstreamBWout;
  p|PCdelayBWSend;
  p|gSpaceSum;
  p|gpriority;
  p|expectOrthoT;
  p|phantomSym;
  p|amPhantom;
  p|useBWBarrier;
  p|numOrthoCol;
  p|numOrthoRow;
  p|numOrtho;
  if (p.isUnpacking())
  {
      mynewData=NULL;
      othernewData=NULL;
      resultCookies=new CkSectionInfo[grainSize];
      if(notOnDiagonal)
	otherResultCookies= new CkSectionInfo[grainSize];
      else
	otherResultCookies=NULL;
      if(existsOut)
	  outData= new double[grainSizeX*grainSizeY];
      else
	  outData=NULL;
      /// @todo: Fix this to allocate or grab a msgLeft and msgRight. inDataLeft/Right is no longer allocated directly
      if(existsLeft)
	  inDataLeft = new double[2*numExpectedX*numPoints];
      else
	  inDataLeft=NULL;
      if(existsRight)
	  inDataRight = new double[2*numExpectedY*numPoints];
      else
	  inDataRight=NULL;
      orthoCookies= new CkSectionInfo[numOrtho*2];
      orthoCB= new CkCallback[numOrtho*2];
      if(PCstreamBWout && !collectAllTiles)
	{
	  columnCount= new int[numOrthoCol];
	  columnCountOther= new int[numOrthoCol];

	}

  }
  int i;
  for (i=0; i<grainSize; i++) p|resultCookies[i];
  if(notOnDiagonal)
    for (i=0; i<grainSize; i++) p|otherResultCookies[i];
  for (int i=0; i<grainSize; i++)
    CmiAssert(resultCookies[i].get_redNo() > 0);
  if(existsOut)
    p(outData, grainSize*grainSize);
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
  if(PCstreamBWout && !collectAllTiles)
    {
      PUParray(p,columnCount,numOrthoCol);
      PUParray(p,columnCountOther,numOrthoCol);
    }
#ifdef _PAIRCALC_DEBUG_
  if (p.isUnpacking())
    {
      CkPrintf("[%d,%d,%d,%d,%d] pup unpacking on %d resumed=%d memory %d\n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,symmetric,CkMyPe(),resumed, CmiMemoryUsage());
      CkPrintf("[%d,%d,%d,%d,%d] pupped : %d,%d,%d,%d,%d %d %d %d %d  %d %d cb cb_aid %d %d %d cb_lb inDataLeft inDataRight outData  %d \n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numRecd, numExpected, grainSize, numStates, numChunks, numPoints, symmetric, conserveMemory, lbpaircalc, cpreduce, cb_ep, existsLeft, existsRight,  resumed);

    }
  else
    CkPrintf("[%d,%d,%d,%d,%d] pup called on %d\n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,symmetric,CkMyPe());
#endif


}

PairCalculator::~PairCalculator()
{

#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] destructs on [%d]\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric, CkMyPe());
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

  int maxorthostateindex=(numStates/orthoGrainSize-1)*orthoGrainSize;
  //  int orthoIndexX=(msg->orthoX*orthoGrainSize-thisIndex.x)/orthoGrainSize;
  //  int orthoIndexY=(msg->orthoY*orthoGrainSize-thisIndex.y)/orthoGrainSize;
  int orthoIndexX=msg->orthoX*orthoGrainSize;
  orthoIndexX= (orthoIndexX>maxorthostateindex) ? maxorthostateindex : orthoIndexX;
  int orthoIndexY=msg->orthoY*orthoGrainSize;
  orthoIndexY= (orthoIndexY>maxorthostateindex) ? maxorthostateindex : orthoIndexY;
  orthoIndexX=(orthoIndexX-thisIndex.x)/orthoGrainSize;
  orthoIndexY=(orthoIndexY-thisIndex.y)/orthoGrainSize;

  int orthoIndex=orthoIndexX*numOrthoCol+orthoIndexY;

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] initGRed ox %d oy %d oindex %d oxindex %d oyindex %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric,msg->orthoX, msg->orthoY,orthoIndex, orthoIndexX, orthoIndexY);
#endif
 // numOrtho here is numOrtho per sGrain
  CkAssert(orthoIndex<numOrtho*2);
  CkGetSectionInfo(orthoCookies[orthoIndex],msg);
  orthoCB[orthoIndex]=msg->cb;
  mCastGrpIdOrtho=msg->mCastGrpId;
  /*  cpreduce=section;
  if(msg->lbsync)
  {
      int foo=1;
      contribute(sizeof(int), &foo , CkReduction::sum_int, msg->synccb);
  }

  */

  if(!symmetric && ++numRecd==numOrtho)
  {
    struct s_array { //(sendArray, bcastArray)
      int ep; //Entry point to call
      CkGroupID id; //Array ID to call it on
      CkArrayIndexStruct idx; //Index to send to (if any)
    } array, *ap;


      contribute(sizeof(int), &numRecd , CkReduction::sum_int, msg->synccb, instance);
      numRecd=0;
  }

  //  do not delete nokeep msg
}



void PairCalculator::initResultSection(initResultMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================

  CkAssert(msg->offset<grainSize);
  if(msg->dest == thisIndex.x && thisIndex.x != thisIndex.y)
  {
      CkGetSectionInfo(otherResultCookies[msg->offset],msg);

#ifdef _PAIRCALC_DEBUG_SPROXY_
      CkPrintf("[%d,%d,%d,%d,%d] other initResultSection for dest %d offset %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric,msg->dest, msg->offset);
#endif
  }
  else
  {
      CkGetSectionInfo(resultCookies[msg->offset],msg);

#ifdef _PAIRCALC_DEBUG_SPROXY_
      CkPrintf("[%d,%d,%d,%d,%d] initResultSection for dest %d offset %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric,msg->dest, msg->offset);
#endif
  }
  rck++;

  mCastGrpId=msg->mCastGrpId;
  //to force synchronize in lb resumption
  if(msg->lbsync && rck==grainSize)
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
      CkPrintf("[%d,%d,%d,%d,%d] resumes from sync\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric);
#endif
  }
}



void PairCalculator::acceptLeftData(paircalcInputMsg *msg) 
{
    const double *data = msg->data();
    const int numRows  = msg->numRows();
    const int numCols  = msg->numCols();
	/// Assert that data is a valid pointer
	CkAssert(data != NULL);
	/// Assert that numRows is as expected
	CkAssert(numRows == numExpectedX);
	/// Check data validity
	#ifdef _NAN_CHECK_
		for(int i=0; i < numRows*numCols; i++)
			CkAssert(finite(data[i]));
	#endif
    /// Once the basic checks have passed, and if we're debugging print status info
	#ifdef _PAIRCALC_DEBUG_
		CkPrintf("[%d,%d,%d,%d,%d] Received left matrix block of size %d x %d at %p\n",
                                thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numRows,numCols,data);
	#endif

	/// Set member data pertinent to the left block
    msgLeft      = msg;
	inDataLeft   = const_cast<double*> (data);
	existsLeft   = true;
	numRecd     += numRows;
	numPoints    = numCols/2;  ///< @note: GSpace sends a complex for each point, but PC treats them as 2 doubles.
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
                dbgStr<<"["<<thisIndex.w<<","<<thisIndex.x<<","<<thisIndex.y<<","<<thisIndex.z<<","<<symmetric<<"]"
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
    const double *data = msg->data();
    const int numRows  = msg->numRows();
    const int numCols  = msg->numCols();
	/// Assert that data is a valid pointer
	CkAssert(data != NULL);
	/// Assert that numRows is as expected
	CkAssert(numRows == numExpectedY);
	/// Check data validity
	#ifdef _NAN_CHECK_
		for(int i=0; i < numRows*numCols; i++)
			CkAssert(finite(data[i]));
	#endif
    /// Once the basic checks have passed, and if we're debugging print status info
	#ifdef _PAIRCALC_DEBUG_
		CkPrintf("[%d,%d,%d,%d,%d] Received right matrix block of size %d x %d at %p\n",
                                thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numRows,numCols,data);
	#endif

	/// Set member data pertinent to the right block
    msgRight     = msg;
	inDataRight  = const_cast<double*> (data);
	existsRight  = true;
	numRecd     += numRows;
	numPoints    = numCols/2;  ///< @note: GSpace sends a complex for each point, but PC treats them as 2 doubles.
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
                    dbgStr<<"["<<thisIndex.w<<","<<thisIndex.x<<","<<thisIndex.y<<","<<thisIndex.z<<","<<symmetric<<"]"
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
                                thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,
                                numRecd, numExpected, numExpectedX, numExpectedY);
    #endif
    /// Ensure that we're really ready to launch computations
    CkAssert(numRecd == numExpected);
    blkSize = aMsg->blkSize;
    
    // If this is not a PsiV loop, trigger the forward path for just the non-phantom chares
    if(!aMsg->doPsiV)
    {
        // This iteration is a normal loop. Hence normal behavior
        actionType = NORMALPC;
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
                CkPrintf("\nGamma beat OrthoT. Waiting for T to arrive before proceeding with forward path");
        }
        else
        {
            /// Do nothing for the phantom chare, non-psiv loops. Computation will be triggered only in the backward path
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
    numRecd = 0;
    if (!amPhantom)
        isLeftReady = false;
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
  int tilesq=orthoGrainSize*orthoGrainSize;
  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpIdOrtho).ckLocalBranch();
  //  for(int orthoX=0; orthoX<numOrtho; orthoX++)
  //    for(int orthoY=0; orthoY<numOrtho; orthoY++)

  int progcounter=0;
  for(int orthoIndex=0;orthoIndex<numOrtho;orthoIndex++)
      {
	// copy into submatrix, contribute
	// we need to stride by grainSize and copy by orthoGrainSize
	//	int orthoIndex=orthoX*numOrthoCol+orthoY;

	if(touchedTiles[orthoIndex]==tilesq)
	  {
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	    if(symmetric)
	      CkPrintf("[%d,%d,%d,%d,%d]: contributes %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, orthoIndex);
#endif
	    //	    if (flag_dp) {
	    //	      for (int i = 0; i < orthoGrainSize*orthoGrainSize; i++)
	    //		outTiles[orthoIndex][i] *= 2.0;
	    //	    }
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	    int orthoX=orthoIndex/numOrthoCol;
	    int orthoY=orthoIndex%numOrthoCol;
	    char filename[80];
	    snprintf(filename,80,"fwoutTile_%d_%d:",orthoX,orthoY);
	    dumpMatrixDouble(filename, outTiles[orthoIndex], orthoGrainSize, orthoGrainSize,thisIndex.x+orthoX*orthoGrainSize, thisIndex.y+orthoY*orthoGrainSize);
#endif

#ifdef _NAN_CHECK_
	    for(int i=0; i<orthoGrainSize*orthoGrainSize; i++)
	      CkAssert(finite(outTiles[orthoIndex][i]));
#endif

	    mcastGrp->contribute(orthoGrainSize*orthoGrainSize*sizeof(double), outTiles[orthoIndex], sumMatrixDoubleType, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
	    //mcastGrp->contribute(orthoGrainSize*orthoGrainSize*sizeof(double), outTiles[orthoIndex], CkReduction::sum_double, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
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
	      CkPrintf("[%d,%d,%d,%d,%d]: %i not ready with %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, orthoIndex, touchedTiles[orthoIndex]);
#endif
	  }
      }
}


/**
 * Forward path multiply.
 *
 *   * (numExpectedX X numPoints) X (numPoints X numExpectedY) = (numExpectedX X numExpectedY)
 * To make this work, we transpose the first matrix (A).
 In C++ it appears to be:
 * (ydima X ydimb) = (ydima X xdima) X (xdimb X ydimb)
 * Which would be wrong, this works because we're using fortran
 * BLAS, which has a transposed perspective (column major), so the
 * actual multiplication is:
 *
 * (xdima X xdimb) = (xdima X ydima) X (ydimb X xdimb)
 *
 * Since xdima==xdimb==numExpected==grainSize this gives us the
 * solution matrix we want in one step.
 *
 * In the border case it gives numExpectedX x numExpectedY which is
 * what we want on the borders.  In non borders numExpectedX==numExpectedY.
 *
 */
void
PairCalculator::multiplyForward(bool flag_dp)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] PairCalculator::multiplyForward() Starting forward path computations.\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
#endif
  if(!existsOut){
    CkAssert(outData==NULL);
    existsOut=true;
    outData = new double[grainSizeX * grainSizeY];
    bzero(outData, sizeof(double)* grainSizeX * grainSizeY);
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Allocated outData %d * %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,grainSizeX, grainSizeY);
#endif
  }
  char transform='N';
  int doubleN=2*numPoints;
  char transformT='T';
  // A is right B is left C= alpha*right  x left;
  int m_in=numExpectedY;  //rows of op(A) rows of C
  int n_in=numExpectedX;  // columns of op(B) columns of C
  int k_in=doubleN;       // columns of op(A) rows of op (B)
  //int ldc=numExpectedX;   //leading dimension C
  int ldc=numExpectedY;   //leading dimension C
  double alpha=double(1.0);//multiplicative identity
  if (flag_dp)  // scaling factor for psi
    alpha=2.0;
  double beta=double(0.0); // C is unset
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif


  double *matrixA;
  if(!symmetricOnDiagonal)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] matrixA is inDataRight\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
#endif
      matrixA=inDataRight;
    }
  else  //  symmetric diagonal
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] matrixA is inDataLeft\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
#endif
      matrixA=inDataLeft;
      // these are redundant, numExpectedX==numExpectedY
      m_in=numExpectedX;  // columns of op(B) columns of C
      ldc=numExpectedX;   //leading dimension C
    }

#if PC_FWD_DGEMM_SPLIT > 0
  double betap = 1.0;
  int Ksplit_m = gemmSplitFWk;
  int Ksplit   = ( (k_in > Ksplit_m) ? Ksplit_m : k_in);
  int Krem     = (k_in % Ksplit);
  int Kloop    = k_in/Ksplit-1;
#ifdef TEST_ALIGN
  CkAssert((unsigned int)matrixA%16==0);
  CkAssert((unsigned int)inDataLeft %16==0);
  CkAssert((unsigned int)outData%16==0);
#endif
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, Ksplit, alpha, beta, k_in, k_in, ldc);
#endif
  DGEMM(&transformT, &transform, &m_in, &n_in, &Ksplit, &alpha, matrixA , &k_in, inDataLeft, &k_in, &beta, outData, &ldc);
  CmiNetworkProgress();

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif

  for(int i=1;i<=Kloop;i++){
    int off = i*Ksplit;
    int KsplitU = (i==Kloop ? Ksplit+Krem : Ksplit);
    // if(i==Kloop){Ksplit+=Krem;}

#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif
#ifdef TEST_ALIGN
    CkAssert((unsigned int)&(matrixA[off])%16==0);
    CkAssert((unsigned int)&(inDataLeft[off]) %16==0);
    CkAssert((unsigned int)outData%16==0);
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, KsplitU, alpha, beta, k_in, k_in, ldc);
#endif
    DGEMM(&transformT, &transform, &m_in, &n_in, &KsplitU, &alpha, &matrixA[off], &k_in, &inDataLeft[off], &k_in, &betap, outData, &ldc);
    CmiNetworkProgress();

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif

  }//endfor

#else  // not split

  int lda=doubleN;   //leading dimension A
  int ldb=doubleN;   //leading dimension B

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  dumpMatrixDouble("fwlmdata", inDataLeft, numExpectedX, numPoints*2, thisIndex.x, 0);
  if(inDataRight!=NULL)
    dumpMatrixDouble("fwrmdata", inDataRight, numExpectedY, numPoints*2, thisIndex.y, 0);
#endif
  CkAssert(matrixA!=NULL);
  CkAssert(inDataLeft!=NULL);
  CkAssert(outData!=NULL);

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, k_in, k_in, ldc);
#endif

    DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha,
	  matrixA, &lda, inDataLeft, &ldb, &beta, outData, &ldc);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif  // SPLIT

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  dumpMatrixDouble("fwgmodata",outData,grainSizeX, grainSizeY,thisIndex.x, thisIndex.y);
#endif


#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  // do the slicing and dicing to send bits to Ortho
  contributeSubTiles(outData);
#ifdef _CP_SUBSTEP_TIMING_
  if((UpairCalcID1[instance].forwardTimerID>0)||(UpairCalcID2[instance].forwardTimerID>0))
    {
      double pstart=CmiWallTimer();
      if(symmetric)
	contribute(sizeof(double),&pstart,CkReduction::max_double, UpairCalcID1[instance].endTimerCB , UpairCalcID1[instance].forwardTimerID);

      else
	contribute(sizeof(double),&pstart,CkReduction::max_double, UpairCalcID2[instance].endTimerCB , UpairCalcID2[instance].forwardTimerID);

    }
#endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
  if(phantomSym && symmetric && notOnDiagonal) //mirror our data to the phantom
    {
      CkAssert(existsRight);
      //      CkPrintf("[%d,%d,%d,%d,%d] fw sending phantom\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
      paircalcInputMsg *msg2phantom = new (numExpectedY*numPoints, 8*sizeof(int)) paircalcInputMsg(numPoints,0,false,flag_dp,(complex*)inDataRight,false,blkSize,numExpectedY);
      bool prioPhan=false;
      if(prioPhan)
	{
	  CkSetQueueing(msg2phantom, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg2phantom) = 1; // just make it slower than non prioritized
	}
      thisProxy(thisIndex.w,thisIndex.y, thisIndex.x,thisIndex.z).acceptRightData(msg2phantom);
    }

	  /** If this is an asymmetric loop, dynamics case AND Ortho has already sent T, 
	   * call byMultiplyDynOrthoT() as we must also multiply orthoT by Fpsi
	   * 
	   * @note: This if condition originally lived in acceptPairData(). Has been shoveled here 
	   * to reduce branching over there.
	   */
	  if(expectOrthoT && numRecdBWOT==numOrtho)
	  {
	      thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).bwMultiplyDynOrthoT();
	      numRecdBWOT=0;
	  }
}




void
PairCalculator::contributeSubTiles(double *fullOutput)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: contributeSubTiles \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric);
#endif
  // Done:
  // necessary changes:
  // 1: border tiles are not uniform size
  // 2: offset calculation must handle non-uniformity
  // solutions, use sGrainSize for all row skips
  // use orthoGrainSizeX or orthoGrainSizeY as needed for row and column
  // iteration

  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpIdOrtho).ckLocalBranch();

#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  dumpMatrixDouble("fullOutput", fullOutput, grainSizeX, grainSizeY, thisIndex.x, thisIndex.y);
#endif
  double *outTile;
  bool reuseTile = false;
  bool borderX = orthoGrainSizeRemX != 0;
  bool borderY = orthoGrainSizeRemY != 0;
  bool borderXY = borderX && borderY;

  // remainder logic happens if you are a border sGrain
  if(! borderX && ! borderY)
    { // only do once to cut down on new/delete
      reuseTile = true;
      outTile  = new double[orthoGrainSize*orthoGrainSize];
      bzero(outTile,sizeof(double)*orthoGrainSize*orthoGrainSize);
    }
  // forward multiply ldc
  int bigGindex=grainSizeY;

  for(int orthoX = 0; orthoX < numOrthoCol ; orthoX++)
    {
      // advance tilestart to new column
      // only the size of tiles in last column is affected
      int orthoGrainSizeX=(orthoX == numOrthoCol-1) ? orthoGrainSize + orthoGrainSizeRemX : orthoGrainSize;
      int orthoXoff=orthoGrainSize*bigGindex*orthoX;
      for(int orthoY=0; orthoY<numOrthoRow; orthoY++)

	{

	  int orthoYoff=orthoY*orthoGrainSize;
	  int tileStart=orthoYoff+orthoXoff;
	  // only the last row is affected

	  int orthoGrainSizeY=(orthoY==numOrthoRow-1) ? orthoGrainSize+orthoGrainSizeRemY : orthoGrainSize;
	  int tileSize=orthoGrainSizeX*orthoGrainSizeY;
	  int bigOindex=orthoGrainSizeY;
	  int ocopySize=bigOindex*sizeof(double);
	  int orthoIndex=orthoX*numOrthoCol+orthoY;


	  //	  CkPrintf("[%d,%d,%d,%d,%d]: orthoX %d orthoY %d tileStart %d ogx %d ogy %d grainSizeX %d outx %d outy %d ogxrem %d ogyrem %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,orthoX, orthoY,tileStart, orthoGrainSizeX, orthoGrainSizeY,grainSizeX,thisIndex.x+orthoX*orthoGrainSize, thisIndex.y+orthoY*orthoGrainSize, orthoGrainSizeRemX, orthoGrainSizeRemY);
	  if(! reuseTile)
	    {
	      outTile=new double[tileSize];
	      bzero(outTile,sizeof(double)*tileSize);
	    }
	  // copy into submatrix, contribute
	  // we need to stride by grainSize and copy by orthoGrainSize

	  int itileStart=tileStart;
	  CkAssert(orthoIndex<numOrtho);
	  //	  CkAssert(tileSize*8==bigOindex*ocopySize);
	  //	  CkPrintf("itileStart %d + bigGindex %d * (tileSize %d/ bigOindex %d)== %d  grainSizeX %d * grainSizeY %d == %d\n",itileStart,bigGindex,tileSize,bigOindex, itileStart+bigGindex*(tileSize/bigOindex),grainSizeX,grainSizeY,grainSizeX*grainSizeY);
	  //	  int count=0;
	  //	  int copiedsofar=0;
	  for(int ystart=0; ystart<tileSize; ystart+=bigOindex, itileStart+=bigGindex){

	    //	    CkPrintf("count %d ystart %d, itileStart %d tileSize %d itileEnd %d ocopySize %d copiedsofar %d\n",count++, ystart,itileStart, tileSize, itileStart+ocopySize/8,ocopySize,copiedsofar);
	    CmiMemcpy(&(outTile[ystart]),&(fullOutput[itileStart]),ocopySize);
	    //	    copiedsofar+=ocopySize/8;
	    //	    CkAssert(itileStart+ocopySize/8<=grainSizeX*grainSizeY);

	  }
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	  char filename[80];
	  snprintf(filename,80,"fwoutTile_%d_%d:",orthoX,orthoY);
	  dumpMatrixDouble(filename, outTile, orthoGrainSizeX, orthoGrainSizeY,thisIndex.x+orthoX*orthoGrainSize, thisIndex.y+orthoY*orthoGrainSize);
#endif
	  mcastGrp->contribute(tileSize*sizeof(double), outTile, sumMatrixDoubleType, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
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
  //  CkPrintf("[%d,%d,%d,%d,%d] acceptOrthoT, numRecdBWOT now %d \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numRecdBWOT);
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
  int maxorthostateindex=(numStates/orthoGrainSize-1)*orthoGrainSize;
  // find our tile indices within this sGrain
  int orthoX=msg->orthoX*orthoGrainSize;

  int orthoY=msg->orthoY*orthoGrainSize;
  ///? @todo Document this after talking with EB. Shouldnt it be an error if orthoX/Y > maxorthostateindex?
  orthoX= (orthoX>maxorthostateindex) ? maxorthostateindex : orthoX;
  orthoY= (orthoY>maxorthostateindex) ? maxorthostateindex : orthoY;
  orthoX=(orthoX-thisIndex.x)/orthoGrainSize;
  orthoY=(orthoY-thisIndex.y)/orthoGrainSize;

  int orthoGrainSizeY=(orthoY==numOrthoRow-1) ? orthoGrainSize+orthoGrainSizeRemY : orthoGrainSize;

  int orthoGrainSizeX=(orthoX == numOrthoCol-1) ? orthoGrainSize + orthoGrainSizeRemX : orthoGrainSize;
  int matrixSize=grainSizeX*grainSizeY;

  collectTile(false, true, false,orthoX, orthoY, orthoGrainSizeX, orthoGrainSizeY, numRecdBWOT, matrixSize, matrix2, matrix1);
  if ((numRecdBWOT==numOrtho) && (numRecd == numExpected)) ///< @todo: resetting numRecd in launchComputations can have effects here. Fix this or change launchComp()
    { // forward path beat orthoT
      CkPrintf("GAMMA beat orthoT, multiplying\n");
      actionType=0;
      bool myfalse=false;
      thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).multiplyForward(myfalse);
      thisProxy(thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z).bwMultiplyDynOrthoT();
      numRecdBWOT=0;
    }

}



//PairCalculator::multiplyResult(int size, double *matrix1, double *matrix2)
void
PairCalculator::multiplyResultI(multiplyResultMsg *msg)
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
void
PairCalculator::multiplyPsiV()
{
	#ifdef DEBUG_CP_PAIRCALC_PSIV
		CkPrintf("[%d,%d,%d,%d,%d] In multiplyPsiV\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
	#endif

  // If I am a non-phantom chare in the symmetric instance, send a copy of my data to my mirror phantom chare
  if(!amPhantom && phantomSym && symmetric && notOnDiagonal) 
    {
      CkAssert(existsRight);
      //      CkPrintf("[%d,%d,%d,%d,%d] fw sending phantom\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
      paircalcInputMsg *msg2phantom = new (numExpectedY*numPoints, 8*sizeof(int)) paircalcInputMsg(numPoints,0,false,true,(complex*)inDataRight,true,blkSize,numExpectedY);
      bool prioPhan=false;
      if(prioPhan)
	{
	  CkSetQueueing(msg2phantom, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg2phantom) = 1; // just make it slower than non prioritized
	}
      thisProxy(thisIndex.w,thisIndex.y, thisIndex.x,thisIndex.z).acceptRightData(msg2phantom);
    }
  // we do not need to go through the all of multiplyresult for psiv
  // all we really need is the setup for the multiplyHelper

     // call helper function to do the math
  int  size=grainSizeX*grainSizeY;
  bool unitcoef=true;
  // TODO: figure out relationship between n_in k_in and grainSizeX grainSizeY
  int m_in=numPoints*2;   // rows of op(A)==rows C

  int n_in=grainSizeY;     // columns of op(B)==columns C
  int k_in=grainSizeX;     // columns op(A) == rows op(B)
  /*   if(amPhantom)
    {
      n_in=grainSizeX;
      k_in=grainSizeY;
      }*/
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

  sendBWsignalMsg *sigmsg=new (8*sizeof(int)) sendBWsignalMsg;
  if(PCdelayBWSend)
    {
      CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(sigmsg) = 1; // just make it slower
    }
  //collapse this into 1 flag
  if(amPhantom)
    sigmsg->otherdata= true;
  else if(((!phantomSym && symmetric) || !unitcoef) && notOnDiagonal)
    sigmsg->otherdata=true;
  else
    sigmsg->otherdata= false;
  if(gSpaceSum)
    thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResultDirect(sigmsg);
  else
    thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);
}




/**
 * Backward path multiplication
 */
void
PairCalculator::multiplyResult(multiplyResultMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: MultiplyResult from orthoX %d orthoY %d size %d numRecd %d actionType %d amPhantom %d notOnDiagonal %d phantomSym %d symmetricOnDiagonal %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->orthoX, msg->orthoY, msg->size, numRecdBW, msg->actionType, amPhantom, notOnDiagonal, phantomSym, symmetricOnDiagonal);
#endif
#ifdef _CP_SUBSTEP_TIMING_
  if((UpairCalcID1[instance].backwardTimerID>0)||(UpairCalcID1[instance].backwardTimerID>0))
  if(numRecdBW==0)
    {
      double pstart=CmiWallTimer();
      if(symmetric)
	contribute(sizeof(double),&pstart,CkReduction::min_double, UpairCalcID1[instance].beginTimerCB , UpairCalcID1[instance].backwardTimerID);

      else
	contribute(sizeof(double),&pstart,CkReduction::min_double, UpairCalcID2[instance].beginTimerCB , UpairCalcID2[instance].backwardTimerID);
    }
#endif
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


  int maxorthostateindex=(numStates/orthoGrainSize-1)*orthoGrainSize;
  // find our tile indices within this sGrain
  int orthoX=msg->orthoX*orthoGrainSize;

  int orthoY=msg->orthoY*orthoGrainSize;
  orthoX= (orthoX>maxorthostateindex) ? maxorthostateindex : orthoX;
  orthoY= (orthoY>maxorthostateindex) ? maxorthostateindex : orthoY;
  if(amPhantom)
    {
      orthoX=(orthoX-thisIndex.y)/orthoGrainSize;
      orthoY=(orthoY-thisIndex.x)/orthoGrainSize;
      int swap=orthoY;
      orthoY=orthoX;
      orthoX=swap;
      //      CkPrintf("[%d,%d,%d,%d,%d]: phantom MultiplyResult with size %d numRecd %d actionType %d numPoints %d orthoX %d orthoY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size, numRecdBW, msg->actionType, numPoints, orthoX, orthoY);
    }
  else
  {
      orthoX=(orthoX-thisIndex.x)/orthoGrainSize;
      orthoY=(orthoY-thisIndex.y)/orthoGrainSize;
  }

  int orthoGrainSizeY=(orthoY==numOrthoRow-1) ? orthoGrainSize+orthoGrainSizeRemY : orthoGrainSize;

  int orthoGrainSizeX=(orthoX == numOrthoCol-1) ? orthoGrainSize + orthoGrainSizeRemX : orthoGrainSize;
  //  CkPrintf("orthoGrainSizeX %d orthoGrainSizeY %d orthoX %d orthoY %d e1 %.10g\n",orthoGrainSizeX, orthoGrainSizeY, orthoX, orthoY, msg->matrix1[0]);
  //  CkPrintf("orthoGrainSizeX*orthoGrainSizeY = %d msg->size %d\n",orthoGrainSizeY*orthoGrainSizeX, msg->size);
  if(matrix2==NULL||size2<1)
    {
      unitcoef = true;
    }
  if(amPhantom && !existsRight)
    { //our mirror data is delayed
      collectAllTiles=true;
      PCstreamBWout=false;
      CkPrintf("[%d,%d,%d,%d,%d] Warning! phantom got bw before fw, forcing tile collection\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
    }
  int matrixSize=grainSizeX*grainSizeY;

  //ASSUMING TMATRIX IS REAL (LOSS OF GENERALITY)

  double *amatrix=NULL;
  double *amatrix2=matrix2;  // may be overridden later

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
  CkPrintf("orthoGrainSizeX %d orthoGrainSizeY %d orthoX %d orthoY %d e1 %.10g\n",orthoGrainSizeX, orthoGrainSizeY, orthoX, orthoY, msg->matrix1[0]);
  if(grainSize==orthoGrainSize)
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
  if(expectOrthoT && !notOnDiagonal)
    beta=1.0;  // need to subtract fpsi*orthoT
  // default these to 0, will be set for streaming comp if !collectAllTiles
  //BTransform=T offsets for C and A matrices
  int BTCoffset=0;
  int BTAoffset=0;
  //BTransform=N offsets for C and A matrices
  int BNCoffset=0;
  int BNAoffset=0;

  if(numStates==grainSize)// all at once no malloc
    {
      amatrix=matrix1;  // index is 0 in this case, so this is silly
    }

  if (orthoGrainSize==grainSize)
    { // you were sent the correct section only
      amatrix=matrix1;
      // the other tiles were already collected for PSIV
      numRecdBW=numOrtho;
    }
  else if(actionType==PSIV)
    {
      amatrix=matrix1;
      // the other tiles were already collected for PSIV
      numRecdBW=numOrtho;
    }
  else if (collectAllTiles)
    {
      collectTile(true, !unitcoef, false,orthoX, orthoY, orthoGrainSizeX, orthoGrainSizeY, numRecdBW, matrixSize, matrix1, matrix2);
      amatrix = inResult1;
      amatrix2 = inResult2;
    } //else !collectAllTiles
  else
    {  // settings for streaming computation on each tile

      /* For Reference to collect tiles we do this
	 if(symmetric && notOnDiagonal) //swap the non diagonals
	 {  // this is trickier to do than one might expect
	 // because the tiles can have funny shapes
	 bigGindex=grainSizeX;
	 bigOindex=orthoGrainSizeX;
	 orthoXoff=orthoGrainSize*bigGindex*orthoY;
	 orthoYoff=orthoX*orthoGrainSize;

	 }
      // which means we set up our multiply for oY lines of oX size
      // embedded in gY lines of gX size
      */

      amatrix=matrix1;
      if(!unitcoef)
	amatrix2=matrix2;
      // fix n_in and k_in
      n_in=orthoGrainSizeY;
      k_in=orthoGrainSizeX;
      if(symmetric && notOnDiagonal) //swap the non diagonals
	{
	  if(!amPhantom)
	    {
	      int swap=orthoX;
	      orthoX=orthoY;
	      orthoY=swap;
	    }

	  //	  int swap=orthoX;
	  //	  orthoX=orthoY;
	  //	  orthoY=swap;
	  //	  if(amPhantom)
	  // {
	  //    n_in=orthoGrainSizeX;
	  //     k_in=orthoGrainSizeY;
	  //  }
	  //
	  //	  k_in=(orthoX==numOrthoRow-1) ? orthoGrainSize+orthoGrainSizeRemX : orthoGrainSize;

	  //	  n_in=(orthoY == numOrthoCol-1) ? orthoGrainSize + orthoGrainSizeRemY : orthoGrainSize;

	}

      // skip to the rows which apply to this ortho
      BTCoffset=orthoY * m_in * orthoGrainSize;
      BTAoffset=orthoX * m_in * orthoGrainSize;
      BNCoffset=orthoX * m_in * orthoGrainSize;
      BNAoffset=orthoY * m_in * orthoGrainSize;

      /*      if(symmetricOnDiagonal)
	{
	  BNCoffset=orthoY * m_in * orthoGrainSize;
	  BNAoffset=orthoX * m_in * orthoGrainSize;
	}
      */
      beta=1.0;  // need to sum over tiles within orthoY columns
    }

  if(orthoGrainSize==grainSize || numRecdBW==numOrtho || !collectAllTiles || actionType==PSIV)
   { // have all the input we need
     // call helper function to do the math

     if(actionType!=PSIV && !collectAllTiles && n_in*k_in>size)
       {
	 CkPrintf("[%d,%d,%d,%d,%d] Warning! your n_in %d and k_in %d is larger than matrix1->size %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric, n_in,k_in,size);
	 n_in=orthoGrainSize;
	 k_in=orthoGrainSize;
       }
     bwMultiplyHelper(size, matrix1, matrix2, amatrix, amatrix2,  unitcoef, m_in, n_in, k_in, BNAoffset, BNCoffset, BTAoffset, BTCoffset, orthoX, orthoY, beta, orthoGrainSizeX, orthoGrainSizeY);
   }
  //#define SKIP_PARTIAL_SEND
#ifndef SKIP_PARTIAL_SEND
  if(PCstreamBWout && !collectAllTiles && !useBWBarrier && actionType!=PSIV)
    { // send results which are complete and not yet sent
      //      if(symmetric && notOnDiagonal && !amPhantom) //swap the
      //      non diagonals

      /*      if(symmetric && notOnDiagonal &&!amPhantom) //swap the non diagonals
	{
	  k_in=(orthoX==numOrthoRow-1) ? orthoGrainSize+orthoGrainSizeRemX : orthoGrainSize;

	  n_in=(orthoY == numOrthoCol-1) ? orthoGrainSize + orthoGrainSizeRemY : orthoGrainSize;
	  }
      */
      bwSendHelper( orthoX, orthoY, k_in, n_in, orthoGrainSizeX, orthoGrainSizeY);
      //      bwSendHelper( orthoX, orthoY, k_in, n_in, k_in, n_in);
    }

#else // dump them all after multiply complete

  if((PCstreamBWout && !collectAllTiles && !useBWBarrier && actionType!=PSIV) && numRecdBW==numOrtho)
    {
      sendBWsignalMsg *sigmsg=new (8*sizeof(int)) sendBWsignalMsg;
      if(PCdelayBWSend)
	{
	  CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(sigmsg) = 1; // just make it slower
	}
      //collapse this into 1 flag
      if(amPhantom)
	sigmsg->otherdata= true;
      else if(((!phantomSym && symmetric) || !unitcoef) && notOnDiagonal)
	sigmsg->otherdata=true;
      else
	sigmsg->otherdata= false;

      if(gSpaceSum)
	thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResultDirect(sigmsg);
      else
	thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);
      numRecdBW=0;
    }
  else
    {
      CkPrintf("[%d,%d,%d,%d,%d] not sending PCstreamBWout %d collectAllTiles %d useBWBarrier %d actionType %d  numRecdBW %d numOrtho%d \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,PCstreamBWout,collectAllTiles,useBWBarrier,actionType, numRecdBW,numOrtho);
    }
#endif
  // check to see if we're all done
  if(((!PCstreamBWout && collectAllTiles) && (orthoGrainSize==grainSize || numRecdBW==numOrtho)) || (useBWBarrier && (orthoGrainSize==grainSize || numRecdBW==numOrtho))|| actionType==PSIV)
    { // clean up
      if(collectAllTiles || !unitcoef)
	{  // only safe to do this if we allocated them
	  if(conserveMemory>=0)
	    {
	      if(inResult2!=NULL)
		delete [] inResult2;
	      if(inResult1!=NULL)
		delete [] inResult1;
	      inResult1=NULL;
	      inResult2=NULL;
	    }
	}

      if(useBWBarrier){

	int wehaveours=1;
	contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
		   CkCallback(CkIndex_PairCalculator::bwbarrier(NULL),thisProxy));
      }
      else
	{
	  sendBWsignalMsg *sigmsg=new (8*sizeof(int)) sendBWsignalMsg;
	  if(PCdelayBWSend)
	    {
	      CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
	      *(int*)CkPriorityPtr(sigmsg) = 1; // just make it slower
	    }
	  //collapse this into 1 flag
	  if(amPhantom)
	    sigmsg->otherdata= true;
	  else if(((!phantomSym && symmetric) || !unitcoef) && notOnDiagonal)
	    sigmsg->otherdata=true;
	  else
	    sigmsg->otherdata= false;

	  if(gSpaceSum)
	    thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResultDirect(sigmsg);
	  else
	    thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);
	}
      numRecdBW=0;
      if(PCstreamBWout)
	{
	  bzero(columnCount, sizeof(int) * numOrthoCol);
	  bzero(columnCountOther, sizeof(int) * numOrthoCol);
	}
      if(conserveMemory>0)
	{
	  // clear the right and left they'll get reallocated on the next pass
#ifndef PC_USE_RDMA
	  // we really don't want to reregister this every phase
	  delete msgLeft;
      msgLeft = 0;
	  inDataLeft=NULL;
	  if(!symmetric || (symmetric && notOnDiagonal)) {
	    delete msgRight;
        msgRight = 0;
	    inDataRight = NULL;
	  }
	  existsLeft=false;
	  existsRight=false;  ///< @todo: Shouldnt this be inside the previous if block. Check if it makes a difference while running with conserveMem>0
#endif
	  if(outData!=NULL && actionType!=KEEPORTHO)
	    {
	      delete [] outData;
	      outData = NULL;
	      existsOut=false;
	    }
	}
    }

}



void PairCalculator::collectTile(bool doMatrix1, bool doMatrix2, bool doOrthoT, int orthoX, int orthoY, int orthoGrainSizeX, int orthoGrainSizeY, int numRecdBW, int matrixSize, double *matrix1, double* matrix2)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: collectTile aggregating numRecdBW %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numRecdBW);
#endif


  //do strided copy
  // stride amatrix by grainSize,
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

  int orthoXoff=orthoGrainSize*bigGindex*orthoX;
  int orthoYoff=orthoY*orthoGrainSize;

  if(symmetric && notOnDiagonal) //swap the non diagonals
    {
      if(amPhantom)
	{
	  bigGindex=grainSizeY;
	  bigOindex=orthoGrainSizeY;
	  orthoXoff=orthoGrainSize*bigGindex*orthoX;
	  orthoYoff=orthoY*orthoGrainSize;
	}
      else{
	bigGindex=grainSizeX;
	bigOindex=orthoGrainSizeX;
	orthoXoff=orthoGrainSize*bigGindex*orthoY;
	orthoYoff=orthoX*orthoGrainSize;
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
    dumpMatrixDouble(filename, matrix2, orthoGrainSizeX, orthoGrainSizeY,orthoX*orthoGrainSize, orthoY*orthoGrainSize);
#endif

  }
  double *dest= (doOrthoT) ? outData : inResult1;

  tileStart=savetileStart;

  if(doMatrix1)
    {
#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
  CkPrintf("[%d,%d,%d,%d,%d]: collectTile copying tile bigOindex %d bigGindex %d tileSize %d tileStart %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, bigOindex, bigGindex, tileSize,tileStart);
#endif

      for(int i=0; i< tileSize; i+=bigOindex,tileStart+=bigGindex)
	CmiMemcpy(&(dest[tileStart]),&(matrix1[i]),ocopySize);
#ifdef _PAIRCALC_DEBUG_PARANOID_BW_

      char filename[80];
      if(doOrthoT)
	snprintf(filename,80,"bworthoT_%d_%d:",orthoX,orthoY);
      else
	snprintf(filename,80,"bwinResult1_%d_%d:",orthoX,orthoY);
      dumpMatrixDouble(filename, matrix1, orthoGrainSizeX, orthoGrainSizeY,orthoX*orthoGrainSize, orthoY*orthoGrainSize);
#endif
    }


}



void PairCalculator::bwbarrier(CkReductionMsg *msg)
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



// PRE: amatrix contains the entire matrix for the multiplication.
//      matrix contains just what was sent in the trigger message.
void PairCalculator::bwMultiplyHelper(int size, double *matrix1, double *matrix2, double *amatrix, double *amatrix2, bool unitcoef, int m_in, int n_in, int k_in, int BNAoffset, int BNCoffset, int BTAoffset, int BTCoffset, int orthoX, int orthoY, double beta, int orthoGrainSizeX, int orthoGrainSizeY)
{


#ifdef _PAIRCALC_DEBUG_

  if(numRecdBW==numOrtho|| orthoGrainSize==grainSize || collectAllTiles)
    {
      streamCaughtR++;
    }
  CkPrintf("[%d,%d,%d,%d,%d]: bwMultiplyHelper with size %d numRecd %d actionType %d orthoX %d orthoY %d orthoGrainSizeX %d orthoGrainSizeY %d BTCoffset %d BNCoffset %d m_in %d n_in %d k_in %d iter %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, size, numRecdBW, actionType, orthoX, orthoY,orthoGrainSizeX, orthoGrainSizeY, BTCoffset, BNCoffset, m_in, n_in, k_in, streamCaughtR);
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
  if(orthoGrainSize==grainSize || collectAllTiles)
    {
      dumpMatrixDouble("bwm1cidata",amatrix,grainSizeX,grainSizeY,0,0,0,streamCaughtR);
      if(!unitcoef)
	{ // CG non minimization case
	  dumpMatrixDouble("bwm2cidata",amatrix2,grainSizeX, grainSizeY,0,0,0,streamCaughtR);
	}
    }
#endif
  int  matrixSize=grainSizeX*grainSize;
  if(symmetric && actionType==KEEPORTHO) // there will be a psiV step following
    {

      if(outData==NULL)
	{ //reuse outData to hold onto ortho
	  CkAssert(!existsOut);
	  outData=new double[matrixSize];
	  bzero(outData,sizeof(double)*matrixSize);
	  existsOut=true;
	}
      //keep the orthoT we just received in matrix1
      if(!collectAllTiles && orthoGrainSize!=grainSize)
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
      CkPrintf("[%d,%d,%d,%d,%d] Allocated mynewData %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpectedY,numPoints);
#endif

      mynewData = new complex[numPoints*numExpectedY];
      bzero(mynewData,numPoints*numExpectedY* sizeof(complex));
      existsNew=true;
      if(!amPhantom && ((symmetric || !unitcoef) && notOnDiagonal)){
	if(othernewData==NULL)
	  {
#ifdef _PAIRCALC_DEBUG_
	    CkPrintf("[%d,%d,%d,%d,%d] Allocated other %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpectedX,numPoints);
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
	    CkPrintf("[%d,%d,%d,%d,%d] Allocated other %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpectedY,numPoints);
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

  int chunksize=blkSize/numChunks;
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
      if(symmetricOnDiagonal && orthoX!=orthoY && !collectAllTiles && actionType!=PSIV )
	{
	  // here is where we need some more juggling
	  lk_in=n_in;
	  ln_in=k_in;
	}


#if PC_BWD_DGEMM_SPLIT > 0
      if(symmetric)
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

      if(symmetric) {
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
      if(symmetric)
	dumpMatrixComplex(snark,mynewData,numExpectedY,numPoints,0,ystart,streamCaughtR);
      else
	dumpMatrixComplex(snark,mynewData,numExpectedY,numPoints,0,ystart,streamCaughtR);
#endif
    }// end of !amPhantom

#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  // Multiply to compensate for the missing triangle in symmetric case
  if((amPhantom || (!phantomSym && symmetric && notOnDiagonal)) && existsRight)
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

    CkAssert(collectAllTiles || orthoGrainSize==grainSize);

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
      if(!collectAllTiles && orthoGrainSize!=grainSize)
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
  CkPrintf("[PAIRCALC] [%d,%d,%d,%d,%d] backward gemm out %.10g %.10g %.10g %.10g \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric, mynewDatad[0],mynewDatad[1],mynewData[numPoints*grainSizeX-1].re,mynewData[numPoints*grainSizeX-1].im);
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

    // output modified by subtracting an application of orthoT

    // C = alpha*A*B + beta*C
    // C= -1 * inRight * orthoT + C

  if(othernewData==NULL)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] Allocated other %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpectedX,numPoints);
#endif
      othernewData = new complex[numPoints*numExpectedX];
      bzero(othernewData,numPoints*numExpectedX * sizeof(complex));
      existsNew=true;
    }
  if(mynewData==NULL)
    {
      CkAssert(numPoints>0);
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d,%d,%d,%d,%d] Allocated mynewData %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpectedY,numPoints);
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
    int chunksize=blkSize/numChunks;
    int ystart=chunksize*thisIndex.z;
    if(notOnDiagonal)
      dumpMatrixComplex("bwg2modata",othernewData,numExpectedX,numPoints,0,ystart,streamCaughtR);
#endif

}




void PairCalculator::bwSendHelper(int orthoX, int orthoY, int sizeX, int sizeY, int orthoGrainSizeX, int orthoGrainSizeY)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper with  numRecd %d actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,  numRecdBW, actionType);
#endif


  // Check to see if this column is  complete, if so send it.

  // symmetric uses BNCoffset on newData for diagonal.
  // symmetric also BTCoffset on otherNewData for off diagonal.
  // asymmetric uses BTC offset on newData

  // determine the orthograinsizeXY for the orthoX,orthoY tile
  // do the calculation based on numExpectedX numExpectedY
  // to account for last row col remainder expansion in PC


  int maxorthostateindex=(numStates/orthoGrainSize-1)*orthoGrainSize;

  int orthoXgrain=orthoX*orthoGrainSize;
  int orthoYgrain=orthoY*orthoGrainSize;

  orthoXgrain= (orthoXgrain>maxorthostateindex) ? maxorthostateindex : orthoXgrain;
  orthoYgrain= (orthoYgrain>maxorthostateindex) ? maxorthostateindex : orthoYgrain;

  //  CkAssert(orthoGrainSizeY==sizeY);
  //  CkAssert(orthoGrainSizeX==sizeX);
  if(symmetric)
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
	      CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper !amPhantom orthoXgrain %d sizeX %d orthoYgrain %d sizeY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,  orthoXgrain, sizeX, orthoYgrain, sizeY);
#endif

	      int   startGrain=orthoYgrain;
	      int   endGrain=startGrain+sizeY;
	      if(notOnDiagonal)
		{
		  startGrain=orthoXgrain;
		  endGrain=startGrain+sizeY;
		}
	      // send orthoX in newData
	      if(gSpaceSum)
		sendBWResultColumnDirect(false, startGrain, endGrain);
	      else
		sendBWResultColumn(false, startGrain, endGrain);
	    }
	  CkAssert(columnCount[index]<=numOrthoCol);
	}

    }
  else //asymm
    {
      if(++columnCount[orthoY]==numOrthoCol) // BTC
	{
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper asymm orthoYgrain %d sizeY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,  orthoYgrain, sizeY);
#endif

	  //send orthoY in newdata
	  if(gSpaceSum)
	    sendBWResultColumnDirect(false, orthoYgrain, orthoYgrain+sizeY);
	  else
	    sendBWResultColumn(false, orthoYgrain, orthoYgrain+sizeY);
	}
      CkAssert(columnCount[orthoY]<=numOrthoCol);
    }
  // send the othernewData for symm off diag (including phantom) and
  // asymm off diag dynamics
  if((amPhantom || (!phantomSym && (othernewData!=NULL)&& notOnDiagonal))&& existsRight)
    {
      int index =orthoX;
      if(symmetric)
      	index=orthoY;
      if(++columnCountOther[index]==numOrthoCol) // BTC
	{
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper amPhantom | other orthoXgrain %d sizeX %d orthoYgrain %d sizeY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,  orthoXgrain, sizeX, orthoYgrain, sizeY);
#endif


	    int startGrain=orthoXgrain;
	    int endGrain=sizeX+startGrain;

	    //int startGrain=orthoYgrain;
	    //int endGrain=sizeY+startGrain;
	    if(symmetric)
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
	  CkPrintf("[%d,%d,%d,%d,%d]: bwSendHelper amPhantom | other startGrain %d endGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,  startGrain, endGrain );
#endif
	  //send orthoY in otherNewData
	  if(gSpaceSum)
	    sendBWResultColumnDirect(true, startGrain, endGrain);
	  else
	    sendBWResultColumn(true, startGrain, endGrain);
	}
      CkAssert(columnCountOther[index]<=numOrthoCol);
    }

  // this could be refined to track an array of completed columns
  // and send them in some grouping scheme
  CkAssert(numRecdBW<=numOrtho);
  if(numRecdBW==numOrtho)
    { //all done clean up after ourselves
      //      CkPrintf("[PAIRCALC] [%d,%d,%d,%d,%d] cleaning up \n",
  // thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric);

      if(!amPhantom)
	{
	  if(conserveMemory>=0)
	    {
	      delete [] mynewData;
	      mynewData=NULL;
	      existsNew=false;
	    }
	  else
	    {
	      bzero(mynewData,numPoints*numExpectedY* sizeof(complex));
	    }
	}
      if(othernewData!=NULL)
	{
	  if(conserveMemory>=0)
	    {
#ifdef _PAIRCALC_DEBUG_
	    CkPrintf("[%d,%d,%d,%d,%d] deleting othernewData\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
#endif
	      delete [] othernewData;
	      othernewData=NULL;
	    }
	  else
	    {
	      if(amPhantom)
		bzero(othernewData,numPoints*numExpectedY * sizeof(complex));
	      else
		bzero(othernewData,numPoints*numExpectedX * sizeof(complex));
	    }
	}
	#ifdef PC_USE_RDMA
		// Unless a PsiV step is next, let the collators know that they should now expect the next batch of data via RDMA. If a PsiV step is next, then PsiV data will come in via traditional messages. expectNext() will be called later at the end of multiplyPsiV().
		if ( !(symmetric && actionType == KEEPORTHO) )
		{
			leftCollator->expectNext();
			rightCollator->expectNext();
		}
		else
		{
			#ifdef DEBUG_CP_PAIRCALC_PSIV
				CkPrintf("[%d,%d,%d,%d,%d]: Am NOT notifying the message handlers to expectNext() as a PsiV step is next (actionType=%d). Data should be arriving in messages. \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,  numRecdBW, actionType);
			#endif
		}
	#endif

      bzero(columnCount, sizeof(int) * numOrthoCol);
      bzero(columnCountOther, sizeof(int) * numOrthoCol);
      numRecdBW=0;
#ifdef _CP_SUBSTEP_TIMING_
      if((UpairCalcID1[instance].backwardTimerID>0)||(UpairCalcID2[instance].backwardTimerID>0))
	{
	  double pstart=CmiWallTimer();
	  if(symmetric)
	    contribute(sizeof(double),&pstart,CkReduction::max_double, UpairCalcID1[instance].endTimerCB , UpairCalcID1[instance].backwardTimerID);
	  else
	    contribute(sizeof(double),&pstart,CkReduction::max_double, UpairCalcID2[instance].endTimerCB , UpairCalcID2[instance].backwardTimerID);
	}
#endif

    }
}




void
PairCalculator::sendBWResultColumnDirect(bool otherdata, int startGrain, int endGrain  )
{
#ifdef _PAIRCALC_DEBUG_CONTRIB_
  CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumnDirect with otherdata %d actionType %d startGrain %d endGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, otherdata, actionType, startGrain, endGrain);
#endif


  int cp_entry=cb_ep;
  if(actionType==PSIV)
    {
      cp_entry= cb_ep_tol;
    }
  CkAssert(endGrain<=numStates);
  int numToSend=endGrain-startGrain;
  int permuter=(int) ((float) thisIndex.z/ (float) numChunks) * (float) numToSend;
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
	if(gpriority)
	  {
	    *((int*)CkPriorityPtr(msg)) = gpriority;
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  }
	msg->init(thisIndex,numPoints, thisIndex.z, computed);

#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d]: sending partial other of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,jPermuted,index+jPermuted,thisIndex.w);
#endif

#ifdef _NAN_CHECK_
	for(int i=0;i<msg->N ;i++)
	  {
	    if(!finite(msg->result[i].re)|| !finite(msg->result[i].im))
	      {
		CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
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
		fprintf(stderr,"[%d,%d,%d,%d,%d]: sendBWResultColumnDirect nan in computed at %d %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, jPermuted,i);
		CkAssert(finite(computed[i].re));
		CkAssert(finite(computed[i].im));
	      }
	  }
#endif
	int index=thisIndex.y;
	CkCallback mycb(cp_entry, CkArrayIndex2D(jPermuted+index ,thisIndex.w), cb_aid);
	partialResultMsg *msg=new (numPoints, 8*sizeof(int) ) partialResultMsg;
	if(gpriority)
	  {
	    *((int*)CkPriorityPtr(msg)) = gpriority;
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  }
	msg->init(thisIndex,numPoints, thisIndex.z, computed);
	/*
	  msg->N=numPoints;
	  msg->myoffset = thisIndex.z; // chunkth
	  CmiMemcpy(msg->result,mynewData+jPermuted*numPoints,msg->N*sizeof(complex));
	*/
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d]:sending partial of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,jPermuted,index+jPermuted,thisIndex.w);
#endif

#ifdef _NAN_CHECK_
	for(int i=0;i<msg->N ;i++)
	  {
	    if(!finite(msg->result[i].re)|| !finite(msg->result[i].im))
	      {
		CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
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
  CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultColumn with actionType %d startGrain %d sendGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType, startGrain, endGrain);
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
	CkPrintf("[%d,%d,%d,%d,%d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
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
	  CkPrintf("[%d,%d,%d,%d,%d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
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
		CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultDirect with actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType);
	#endif
	// Now we have results in mynewData and if(symmetric) othernewData
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
			if(gpriority)
			{
				*((int*)CkPriorityPtr(omsg)) = gpriority;
				CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
			}
			omsg->init(thisIndex,numPoints, thisIndex.z, computed);
			#ifdef _PAIRCALC_DEBUG_
				CkPrintf("[%d,%d,%d,%d,%d]:sending other partial of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,j,index+j,thisIndex.w);
			#endif
			#ifdef _NAN_CHECK_
				for(int i=0;i<omsg->N ;i++)
				{
					if(!finite(omsg->result[i].re) || !finite(omsg->result[i].im))
						CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
					CkAssert(finite(omsg->result[i].re));
					CkAssert(finite(omsg->result[i].im));
				}
			#endif
			mycb.send(omsg);
		}
		if(otherdata)delete [] othernewData;
		othernewData=NULL;
		existsNew=false;
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
			if(gpriority)
			{
				*((int*)CkPriorityPtr(omsg)) = gpriority;
				CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
			}
			omsg->init(thisIndex,numPoints, thisIndex.z, computed);
			#ifdef _PAIRCALC_DEBUG_
				CkPrintf("[%d,%d,%d,%d,%d]:sending partial of size %d offset %d to [%d %d]\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,numPoints,j,index+j,thisIndex.w);
			#endif
			#ifdef _NAN_CHECK_
				for(int i=0;i<omsg->N ;i++)
				{
					if(!finite(omsg->result[i].re) || !finite(omsg->result[i].im))
						CkPrintf("[%d,%d,%d,%d,%d]: sendBWResultDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
					CkAssert(finite(omsg->result[i].re));
					CkAssert(finite(omsg->result[i].im));
				}
			#endif
			mycb.send(omsg);
		}
		if(mynewData!=NULL)
			delete [] mynewData;
		mynewData=NULL;
		existsNew=false;
	}
}




// entry method to allow us to delay this outbound communication
// to minimize brain dead BG/L interference we have a signal to prioritize this
void
PairCalculator::sendBWResult(sendBWsignalMsg *msg)
{

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d]: sendBWResult with actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType);
#endif
  bool otherdata=msg->otherdata;
  delete msg;
  // Now we have results in mynewData and if(symmetric) othernewData
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
	CkPrintf("[%d,%d,%d,%d,%d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
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
	  CkPrintf("[%d,%d,%d,%d,%d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
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
  if(!amPhantom)
    {
      if(mynewData!=NULL)
	delete [] mynewData;
      mynewData=NULL;
      existsNew=false;
    }
  if(otherdata)
    delete [] othernewData;
  othernewData=NULL;
}



void PairCalculator::dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim, int xstart, int ystart, int xtra1, int xtra2)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d_%d_%d.out",999);
  sprintf(filename, fmt, thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,
	  xtra1, xtra2, symmetric);
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
  sprintf(filename, fmt, thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, iter ,symmetric);
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
	//	if(symmetric)
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


#include "paircalcMessages.def.h"
#include "InputDataHandler.h"
#include "inputDataHandler.def.h"
#include "ckPairCalculator.def.h"

