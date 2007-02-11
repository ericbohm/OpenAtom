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
 * backward path of the asymmetric pairculator we will receive 2 input
 * matrices, gamma and orthoT.  Where orthoT came from the
 * orthonormalization following the symmetric paircalculator.  And
 * gamma was produced by a multiplication with orthoT in the Ortho
 * object.
 *
 * If ortho falls out of tolerance then Ortho will signal the GSP that
 * a tolerance update is needed.  We then proceed with the psi
 * calculation as normal.  On receipt of newpsi, Gspace will then
 * react by sending the PC the Psi velocities (PsiV) in the same way
 * (acceptPairData) that it sends Psi, but with the psiv flag set
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
 * paircalculator decomposition.  
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
 * to the paircalculator.  
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
 * sGrainSize. sGrainSize must be an even mod of the number of states.
 * The former restriction avoids overlap and boundary issues for tiles
 * between ortho and the calculator.  The former avoids some rather
 * gnarly issues in the setup of multicasts and reductions. The latter
 * avoids the necessity of remainder logic and irregularity in the
 * state decomposition.  Strictly speaking the latter restriction
 * rules out problems with an prime number of states and diminishes
 * decomposition flexibility for problems with poorly factorized
 * numbers of states.  In practice this isn't a problem because their
 * is enough flexibility in creating the simulation system such that
 * users can be taught to avoid creating such systems.
 *
 *
 * NOTE: The magic number 2 appears in 2 contexts.  Either we have
 * twice as many inputs in the non diagonal elements of the matrix.
 * Or in the case of transforming our arrays of complex into arrays of
 * doubles.  This transformation is done so we can use DGEMM instead
 * of ZGEMM.  Movated by the empirical discovery that BLAS
 * implementors do all their work on DGEMM and leave ZGEMM out in the
 * unoptimized cold.  The latter issue crops up everywhere we have to
 * do allocation or manipulations of the input.  Could arguably be
 * abstracted away by absorbing it into a size field.
*/
//**************************************************************************


#include "ckPairCalculator.h"
#include "pairCalculator.h"
#ifdef CMK_VERSION_BLUEGENE
#include "builtins.h"
// void __alignx(int alignment,  const void *address);
// void _alignx(int alignment,  const void *address);
// void alignx(int alignment,  const void *address);
#endif
//#define PRINT_DGEMM_PARAMS

ComlibInstanceHandle mcastInstanceCP;
ComlibInstanceHandle mcastInstanceACP;

extern PairCalcID pairCalcID1;
extern PairCalcID pairCalcID2;

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
#ifdef CMK_VERSION_BLUEGENE
      __alignx(16,ret);
#endif
  int size0=msgs[0]->getSize();
  int size=size0/sizeof(double);

  double *inmatrix;
  //  int progcount=0;
  if(nMsg>3) // switch loops and unroll
    {  
      int i=1;
      // idea here is to have only 1 store for 4 loads
#ifdef CMK_VERSION_BLUEGENE
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
#ifdef CMK_VERSION_BLUEGENE
      __alignx(16,inmatrix);
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

PairCalculator::PairCalculator(bool sym, int grainSize, int s, int numChunks, CkCallback cb, CkArrayID cb_aid, int _cb_ep, int _cb_ep_tol, bool conserveMemory, bool lbpaircalc,  redtypes _cpreduce, int _orthoGrainSize, bool _collectTiles, bool _PCstreamBWout, bool _PCdelayBWSend, int _streamFW, bool _gSpaceSum, int _gpriority, bool _phantomSym, bool _useBWBarrier, int _gemmSplitFWk, int _gemmSplitFWm, int _gemmSplitBW)
{
#ifdef _PAIRCALC_DEBUG_PLACE_
  CkPrintf("[PAIRCALC] [%d %d %d %d %d] inited on pe %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,sym, CkMyPe());
#endif
  this->conserveMemory=conserveMemory;
  this->symmetric = sym;
  this->grainSize = grainSize;
  this->numStates = s;
  this->numChunks = numChunks;
  this->numPoints = -1;
  this->cb_aid = cb_aid;
  this->cb_ep = _cb_ep;
  this->cb_ep_tol = _cb_ep_tol;
  useBWBarrier=_useBWBarrier;
  PCstreamBWout=_PCstreamBWout;
  PCdelayBWSend=_PCdelayBWSend;
  orthoGrainSize=_orthoGrainSize;
  collectAllTiles=_collectTiles;
  streamFW=_streamFW;
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
  numRecRight = 0;
  numRecLeft = 0;
  streamCaughtR=0;
  streamCaughtL=0;
  numExpected = grainSize;
  cpreduce=_cpreduce;
  resumed=true;

  touchedTiles=NULL;
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
  if(phantomSym && thisIndex.y<thisIndex.x && symmetric )
    amPhantom=true;
  else
    amPhantom=false;
  if(lbpaircalc)
      setMigratable(true);
  else
      setMigratable(false);
  resultCookies=NULL;
  otherResultCookies=NULL;
  resultCookies=new CkSectionInfo[grainSize];
  int numOrthoCol=grainSize/orthoGrainSize;
  int numOrtho=numOrthoCol*numOrthoCol;
  orthoCookies=new CkSectionInfo[numOrtho];
  orthoCB=new CkCallback[numOrtho];
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
  if(streamFW>0)
    {
      LeftOffsets= new int[numExpected];
      RightOffsets= new int[numExpected];
      LeftRev = new int[numExpected];
      RightRev = new int[numExpected];
      touchedTiles= new int[numOrtho];
      bzero(touchedTiles,numOrtho*sizeof(int));
      outTiles = new double *[numOrtho];
      for(int i=0;i<numOrtho;i++)
	outTiles[i]= new double[orthoGrainSize*orthoGrainSize];
    }
  if(thisIndex.x != thisIndex.y)  
    // we don't actually use these in the asymmetric minimization case
    otherResultCookies=new CkSectionInfo[grainSize];

}

void
PairCalculator::pup(PUP::er &p)
{
  ArrayElement4D::pup(p);
  p|numRecd;
  p|numRecdBW;
  p|numExpected;
  p|grainSize;
  p|numStates;
  p|numChunks;
  p|numPoints;
  p|symmetric;
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
  p|PCstreamBWout;
  p|PCdelayBWSend;
  p|gSpaceSum;
  p|gpriority;
  p|phantomSym;
  p|amPhantom;
  p|useBWBarrier;
  int numOrthoCol=grainSize/orthoGrainSize;
  int numOrtho=numOrthoCol*numOrthoCol;
  if (p.isUnpacking()) 
  {
      mynewData=NULL;
      othernewData=NULL;
      resultCookies=new CkSectionInfo[grainSize];
      if(thisIndex.x != thisIndex.y)
	otherResultCookies= new CkSectionInfo[grainSize];
      else
	otherResultCookies=NULL;
      if(existsOut)
	  outData= new double[grainSize*grainSize];
      else
	  outData=NULL;
      if(existsLeft)
	  inDataLeft = new double[2*numExpected*numPoints];
      else
	  inDataLeft=NULL;
      if(existsRight)
	  inDataRight = new double[2*numExpected*numPoints];
      else
	  inDataRight=NULL;
      orthoCookies= new CkSectionInfo[numOrtho];
      orthoCB= new CkCallback[numOrtho];
      if(PCstreamBWout && !collectAllTiles) 
	{
	  columnCount= new int[numOrthoCol];
	  columnCountOther= new int[numOrthoCol];

	}

  }
  int i;
  for (i=0; i<grainSize; i++) p|resultCookies[i];
  if(thisIndex.x != thisIndex.y)
    for (i=0; i<grainSize; i++) p|otherResultCookies[i];
  for (int i=0; i<grainSize; i++)
    CmiAssert(resultCookies[i].get_redNo() > 0);
  if(existsOut)
    p(outData, grainSize*grainSize);
  if(existsLeft)
    p(inDataLeft, numExpected * numPoints * 2);
  if(existsRight)
    p(inDataRight, numExpected* numPoints * 2);
  PUParray(p,orthoCookies,numOrtho);
  PUParray(p,orthoCB,numOrtho);
  if(PCstreamBWout && !collectAllTiles) 
    {
      PUParray(p,columnCount,numOrthoCol);
      PUParray(p,columnCountOther,numOrthoCol);
    }
#ifdef _PAIRCALC_DEBUG_
  if (p.isUnpacking())
    {
      CkPrintf("[%d,%d,%d,%d,%d] pup unpacking on %d resumed=%d memory %d\n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,symmetric,CkMyPe(),resumed, CmiMemoryUsage());
      CkPrintf("[%d,%d,%d,%d,%d] pupped : %d %d %d %d %d %d %d %d %d  %d %d cb cb_aid %d %d %d cb_lb inDataLeft inDataRight outData  %d \n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numRecd, numExpected, grainSize, numStates, numChunks, numPoints, symmetric, conserveMemory, lbpaircalc, cpreduce, cb_ep, existsLeft, existsRight,  resumed);

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
  if(inDataLeft!=NULL)
    delete [] inDataLeft;
  if(inDataRight!=NULL)
    delete [] inDataRight;
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
  if(thisIndex.x != thisIndex.y && otherResultCookies !=NULL)
    delete [] otherResultCookies;
}


// initialize the section cookie and the reduction client
void PairCalculator::initGRed(initGRedMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
 
  int numOrtho=grainSize/orthoGrainSize;
  int orthoIndexX=(msg->orthoX*orthoGrainSize-thisIndex.x)/orthoGrainSize;
  int orthoIndexY=(msg->orthoY*orthoGrainSize-thisIndex.y)/orthoGrainSize;
  int orthoIndex=orthoIndexX*numOrtho+orthoIndexY;
  //int orthoIndex=orthoIndexX*numOrtho+orthoIndexY;
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] initGRed ox %d oy %d oindex %d oxindex %d oyindex %d\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric,msg->orthoX, msg->orthoY,orthoIndex, orthoIndexX, orthoIndexY);
#endif

  CkAssert(orthoIndex<numOrtho*numOrtho);
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
  if(!symmetric && ++numRecd==numOrtho*numOrtho)
    {
      //      CkPrintf("[%d,%d,%d,%d,%d] contributes to doneInit with %d numRecd \n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric,numRecd);
      contribute(sizeof(int), &numRecd , CkReduction::sum_int, msg->synccb);
      numRecd=0;
    }
  
  //  do not delete nokeep msg
}

//!initialize the section cookie for each slice of the result
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


/**
 * Forward path multiply.  Accumulates rows and columns when all
 * arrive it calls the multiply
 */
void
PairCalculator::acceptPairData(calculatePairsMsg *msg)
{
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
 
#ifdef _PAIRCALC_DEBUG_

  CkPrintf(" symm=%d    pairCalc[%d %d %d %d] got from [%d %d] with size {%d}, from=%d, count=%d, resumed=%d\n", symmetric, thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,  thisIndex.w, msg->sender, msg->size, msg->fromRow, numRecd,resumed);
#endif
#ifdef _CP_SUBSTEP_TIMING_
  if(numRecd==0)
    {
      double pstart=CmiWallTimer();
      if(symmetric)
	contribute(sizeof(double),&pstart,CkReduction::min_double, pairCalcID1.beginTimerCB , pairCalcID1.forwardTimerID);
      else
	contribute(sizeof(double),&pstart,CkReduction::min_double, pairCalcID2.beginTimerCB , pairCalcID2.forwardTimerID);


    }
#endif
  if(!resumed)
    {
      ResumeFromSync();
    }
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size;i++)
    {
      CkAssert(finite(msg->points[i].re));
      CkAssert(finite(msg->points[i].im));
    }
#endif
  blkSize=msg->blkSize;
  numRecd++;   // increment the number of received counts
  int offset = -1;
  if (msg->fromRow) {   // This could be the symmetric diagonal case
    offset = msg->sender - thisIndex.x;
    if (!existsLeft)
      { // now that we know N we can allocate contiguous space
	CkAssert(inDataLeft==NULL);
	existsLeft=true;
	numPoints = msg->size; // numPoints is init here with the size of the data chunk.
	
	inDataLeft = new double[numExpected*numPoints*2];
	bzero(inDataLeft,numExpected*numPoints*2*sizeof(double));
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d] Allocated Left %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric, numExpected,numPoints);
#endif

      }
    if(streamFW>0 && allCaughtLeft==NULL)
      {
	numRecLeft=0;
	streamCaughtL=0;
	allCaughtLeft=inDataLeft;
      }
    CkAssert(numPoints==msg->size);
    CkAssert(offset<numExpected);
    if(streamFW>0)
      { //record offset, copy data
	memcpy(&(allCaughtLeft[numRecLeft *numPoints*2]), msg->points, numPoints * 2 *sizeof(double));
	LeftOffsets[numRecLeft]=offset;
	LeftRev[offset]=numRecLeft++;
	streamCaughtL++;
      }
    else
      {
	memcpy(&(inDataLeft[offset*numPoints*2]), msg->points, numPoints * 2 *sizeof(double));
      }
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Copying into offset*numPoints %d * %d numPoints *2 %d points start %.12g end %.12g\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z, symmetric, offset, numPoints, numPoints*2,msg->points[0].re, msg->points[numPoints-1].im);
#endif
  }
  else { //not from row
    offset = msg->sender - thisIndex.y;
    if (!existsRight)
      { // now that we know numPoints we can allocate contiguous space
	CkAssert(inDataRight==NULL);
	existsRight=true;
	numPoints = msg->size; // numPoints is init here with the size of the data chunk.
	inDataRight = new double[numExpected*numPoints*2];
	bzero(inDataRight,numExpected*numPoints*2*sizeof(double));
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d] Allocated right %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpected,numPoints);
#endif
      }
    CkAssert(numPoints==msg->size);
    CkAssert(offset<numExpected);
    if(streamFW>0 && allCaughtRight==NULL)
      {
	numRecRight=0;
	streamCaughtR=0;
	allCaughtRight=inDataRight;
      }


    if(streamFW>0)
      { //record offset, copy data
	memcpy(&(allCaughtRight[numRecRight *numPoints*2]), msg->points, numPoints * 2 *sizeof(double));
	RightOffsets[numRecRight]=offset;
	RightRev[offset]=numRecRight++;
	streamCaughtR++;
      }
    else
      {
	memcpy(&(inDataRight[offset*numPoints*2]), msg->points, numPoints * 2 *sizeof(double));
      }

#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Copying into offset*numPoints %d * %d numPoints *2 %d points start %.12g end %.12g\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,offset,numPoints, numPoints*2,msg->points[0].re, msg->points[numPoints-1].im);
#endif

  }
  if(streamFW)
    CkAssert(numRecRight+numRecLeft==numRecd);

  /*
   *  NOTE: For this to work the data chunks of the same plane across
   *  all states must be of the same size
   */
  // copy the input into our matrix
  /*
   * Once we have accumulated all rows  we gemm it.
   */
  bool streamready;

  if(symmetric && thisIndex.x==thisIndex.y) //only left
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

  if(streamready || ((streamFW>0) && (numRecd == numExpected * 2 || (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected))) && (!msg->doPsiV) )
    {
	multiplyForwardStream(msg->flag_dp);	
	// not yet supported for dynamic psiV
    }
  else if (numRecd == numExpected * 2 || (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected)) 
    {
    if(!msg->doPsiV)
      {  //normal behavior
	actionType=0;
	multiplyForward(msg->flag_dp);	
      }
    else
      {
	// tolerance correction psiV
	multiplyPsiV();
      }
    }
  
  //nokeep message not deleted
}

/**
 * Forward path multiply.  Streaming the computation by computing on
 * chunks as they arrive.
 *
 * One could naively capture the idea as:

 * (streamCaughtR X numPoints) * (numPoints X streamCaughtL) = (streamCaughtR X streamCaughtL)
 * 
 *
 * But we actually need to multiply previous R by new L, previous L by
 * new R new R by new L to generate all values in the solution matrix.
 * This can be accomplished in exactly 2 multiplies.
 *
 * New is appended to old to create AllCaught.
 *
 * In the 2 multiplies the new * new must only happen once.  So we use
 * allcaught in first multiply and oldcaught in second.
 *
 * Asymmetric is: 
 *
 * (allCaughtR X numPoints) * (numPoints X streamCaughtL) =
 * (allCaughtR X streamCaughtL)
 * AND
 * (streamCaughtR X numPoints) * (numPoints X oldCaughtL) =
 * (streamCaughtR X oldCaughtL) 
 *
 * vs the symmetric case which is:
 * 
 * (streamCaughtL X numPoints) * (numPoints X allCaughtL) =
 * (streamCaughtL X allCaughtL) 
 * AND
 * (oldCaughtL X numPoints) * (numPoints X streamCaughtL) =
 * (oldCaughtL X streamCaughtL) 
 *
 * So the sizes of the dgemms rise with each hop in the stream.
 * 
 * These results must then be copied to the correct rows and columns
 * in the tile output.  As each tile is filled we can contribute it
 * and thereby get some streaming output.  This is slightly different
 * from the backward path because here our stream begins as a trickle
 * and ends as a river.  Which is still an improvement over the old
 * mechanism which was just a flood.
 *
*/
void
PairCalculator::multiplyForwardStream(bool flag_dp)
{
  //  CkPrintf("[%d,%d,%d,%d,%d] multi fw stream\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);


  int actualPoints= numPoints*2;
  char transform='N';
  char transformT='T';
  int m_in;              // rows of op(A)==rows C	 
  int n_in;              // columns of op(B)==columns C 
  int k_in=actualPoints; // columns op(A) == rows op(B) 
  int lda=actualPoints;  //leading dimension A
  int ldb=actualPoints;  //leading dimension B
  int ldc;               //leading dimension C  <--- m_in

  double alpha=double(1.0);
  if(symmetric && flag_dp)
    alpha=2.0;
  double beta=double(0.0); // C is unset
  double *outData1=NULL;
  double *outData2=NULL;
  // if asym deal with right by left
  int oldCaughtLeft=numRecLeft-streamCaughtL;
  int oldCaughtRight=numRecRight-streamCaughtR;
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

  if(!symmetric || (thisIndex.x !=thisIndex.y))
  {
    CkAssert(streamCaughtL>0 || (streamCaughtR>0&&oldCaughtLeft>0));

    // As these multiplies create two different shapes of disjoint results
    // we'll be saner if we allocate the results in two different arrays

    m_in= numRecRight;
    n_in= streamCaughtL;
    ldc = m_in; 

    double *leftNewTemp = &(allCaughtLeft[oldCaughtLeft*actualPoints]);

    //multiply all right with new left
    // if we have new left data
    if(streamCaughtL>0)
      {
	outData1= new double[m_in*n_in];
	bzero(outData1,sizeof(double)*m_in*n_in);

#if PC_FWD_DGEMM_SPLIT > 0
        dgemmSplitFwdStreamMK(m_in, n_in, k_in, &transform, &transformT, &alpha, allCaughtRight, &lda, leftNewTemp, &ldb, outData1, &ldc);
#else  // not split

#ifndef CMK_OPTIMIZE
        StartTime=CmiWallTimer();
#endif

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, lda, ldb, ldc);
#endif
	DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, allCaughtRight, &lda, leftNewTemp, &ldb, &beta, outData1, &ldc);
#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
	//kick off progress before next dgemm
	CmiNetworkProgress();
#endif // end of split 

	copyIntoTiles(outData1, outTiles, n_in, m_in,  RightOffsets, &(LeftOffsets[oldCaughtLeft]),  touchedTiles, orthoGrainSize, grainSize / orthoGrainSize);
      }

    //multiply new right with old left if they exist
    if(oldCaughtLeft && streamCaughtR){
      m_in= streamCaughtR;
      n_in= oldCaughtLeft;
      ldc = m_in; 

      // right newTemp is last bit of old
      double *rightNewTemp = &(allCaughtRight[(oldCaughtRight)*actualPoints]);
      outData2= new double[m_in*n_in];
      bzero(outData2,sizeof(double)*m_in*n_in);

#if PC_FWD_DGEMM_SPLIT > 0
      dgemmSplitFwdStreamNK(m_in, n_in, k_in, &transform, &transformT, &alpha, rightNewTemp, &lda, allCaughtLeft, &ldb, outData2, &ldc);
#else  // not split

#ifndef CMK_OPTIMIZE
      StartTime=CmiWallTimer();
#endif
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, lda, ldb, ldc);
#endif
      DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, rightNewTemp, &lda, allCaughtLeft, &ldb, &beta, outData2, &ldc);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif // end of split

      copyIntoTiles(outData2, outTiles, n_in, m_in, &(RightOffsets[oldCaughtRight]), LeftOffsets, touchedTiles, orthoGrainSize, grainSize / orthoGrainSize);
    }
  }
  else // in sym there is only left by left
    {
      // multiply all left by new left
      // left newTemp is last bit of old
      m_in= numRecLeft;
      n_in= streamCaughtL;
      ldc = m_in; 
      outData1= new double[m_in*n_in];
      bzero(outData1,sizeof(double)*m_in*n_in);
      double *leftNewTemp = &(allCaughtLeft[oldCaughtLeft*actualPoints]);

#if PC_FWD_DGEMM_SPLIT > 0
      dgemmSplitFwdStreamMK(m_in, n_in, k_in, &transform, &transformT, &alpha, allCaughtLeft, &lda, leftNewTemp, &ldb, outData1, &ldc);
#else  // not split

#ifndef CMK_OPTIMIZE
      StartTime=CmiWallTimer();
#endif
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, lda, ldb, ldc);
#endif
      DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, allCaughtLeft, &lda, leftNewTemp, &ldb, &beta, outData1, &ldc);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
      //kick off progress before next dgemm
      CmiNetworkProgress();
#endif  // end of split

      copyIntoTiles(outData1, outTiles, n_in, m_in, LeftOffsets, &(LeftOffsets[oldCaughtLeft]), touchedTiles, orthoGrainSize, grainSize / orthoGrainSize);
      
      // multiply new left by old left
      if(oldCaughtLeft)
	{
	  m_in= streamCaughtL;
	  n_in= oldCaughtLeft;
	  ldc = m_in; 
	  outData2= new double[m_in*n_in];
	  bzero(outData2,sizeof(double)*m_in*n_in);
	  // oldCaught is the same pointer as Allcaught, we just decrease n.

#if PC_FWD_DGEMM_SPLIT > 0
	  dgemmSplitFwdStreamNK(m_in, n_in, k_in, &transform, &transformT, &alpha, leftNewTemp, &lda, allCaughtLeft, &ldb, outData2, &ldc);
#else 

#ifndef CMK_OPTIMIZE
	  StartTime=CmiWallTimer();
#endif
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, lda, ldb, ldc);
#endif
	  DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, leftNewTemp, &lda, allCaughtLeft, &ldb, &beta, outData2, &ldc);

#ifndef CMK_OPTIMIZE
            traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif

#endif // split 

	  copyIntoTiles(outData2, outTiles, n_in, m_in, &(LeftOffsets[oldCaughtLeft]), LeftOffsets, touchedTiles, orthoGrainSize, grainSize / orthoGrainSize);
	}
    }

  // check all touched tiles for completeness, send completed tiles
  sendTiles(flag_dp);
  // cleanup
  if(outData1!=NULL)
    delete [] outData1;
  if(outData2!=NULL)
    delete [] outData2;
  streamCaughtL=0;
  streamCaughtR=0;
  if((numRecd == numExpected * 2 || (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected)))
    { // we are done

      int numOrtho=grainSize/orthoGrainSize;
      numOrtho*=numOrtho;

      if(streamFW)
	{
	  double *scratch= new double[actualPoints];

	  if(existsLeft){
	    reorder(LeftOffsets, LeftRev, allCaughtLeft, scratch);
#ifdef _NAN_CHECK_
	  for(int i=0;i<numExpected*numPoints*2;i++)
	    {
	      CkAssert(finite(allCaughtLeft[i]));
	    }
#endif
	  }
	  if(existsRight){
	    reorder(RightOffsets, RightRev, allCaughtRight, scratch);
#ifdef _NAN_CHECK_
	    for(int i=0;i<numExpected*numPoints*2;i++)
	      {
		CkAssert(finite(allCaughtRight[i]));
	      }
#endif
	  }
	  delete [] scratch;
	}
      numRecLeft= numRecRight=numRecd=0;
#ifdef _PAIRCALC_DEBUG_CONTRIB_
      for(int i=0 ; i<numOrtho ; i++)
	{
	  CkAssert(touchedTiles[i]==0);
	}
#endif
      if(phantomSym && symmetric && thisIndex.x!=thisIndex.y) //mirror our data to the phantom
	{
	  CkAssert(existsRight);
	  //	  CkPrintf("[%d,%d,%d,%d,%d] multi fw stream sending phantom\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
	  phantomMsg *phantom=new ( numExpected*numPoints*2, 8*sizeof(int)) phantomMsg;
	  bool prioPhan=false;
	  if(prioPhan)
	    {
	      CkSetQueueing(phantom, CK_QUEUEING_IFIFO);
	      *(int*)CkPriorityPtr(phantom) = 1; // just make it slower
					    // than non prioritized
	    }
	  phantom->init(numExpected*numPoints*2, numPoints, flag_dp, inDataRight, blkSize);
	  thisProxy(thisIndex.w,thisIndex.y, thisIndex.x,thisIndex.z).acceptPhantomData(phantom);

	}
    }
}

void PairCalculator::reorder(int * offsetMap, int *revOffsetMap, double *data, double *scratch)
{

  int actualPoints= numPoints*2;
  // rather ugly, but each iteration correctly places 2 rows
  int datasize=actualPoints*sizeof(double);
  for(int off=0;off<numExpected;off++)
    {
      if(off!=offsetMap[off])
	{  // gotta move
	  int want = revOffsetMap[off];
	  int found = offsetMap[off];
	  int currentOffset = off * actualPoints;
	  int wantOffset = want * actualPoints;
		    
	  memcpy(scratch, &(data[currentOffset]), datasize);
	  memcpy(&(data[currentOffset]), &(data[wantOffset]), datasize);
	  if(want==found) //simple exchange
	    {
	      memcpy(&(data[wantOffset]),scratch, datasize);
	    }
	  else  // three card shuffle
	    {
	      int foundOffset = found * actualPoints;
	      memcpy(&(data[wantOffset]), &(data[foundOffset]), datasize);
	      memcpy(&(data[foundOffset]),scratch, datasize);
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

/**
 * Forward path multiply.  
 *
 *   * (numExpected X numPoints) X (numPoints X numExpected) = (numExpected X numExpected)
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
 */
void
PairCalculator::multiplyForward(bool flag_dp)
{

  CkAssert(numRecd == numExpected * 2 || (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected));
  if(!existsOut){
    CkAssert(outData==NULL);
    existsOut=true;
    outData = new double[grainSize * grainSize];
    bzero(outData, sizeof(double)* grainSize * grainSize);
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Allocated outData %d * %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,grainSize, grainSize);
#endif
  }
  char transform='N';
  int doubleN=2*numPoints;
  char transformT='T';
  int m_in=numExpected;
  int n_in=numExpected;
  int k_in=doubleN;
  int ldc=numExpected;   //leading dimension C
  double alpha=double(1.0);//multiplicative identity
  if (flag_dp)  // scaling factor for psi
    alpha=2.0;
  double beta=double(0.0); // C is unset
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] gemming %c %c %d %d %d %f A %d B %d %f C %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,transformT,transform, m_in, n_in, k_in, alpha, lda, ldb, beta, ldc);
#endif

  double *matrixA;
  if( numRecd == numExpected * 2)
    {
      matrixA=inDataRight;
    }
  else
    {
      matrixA=inDataLeft;
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
  dumpMatrixDouble("fwlmdata", inDataLeft, numExpected, numPoints*2);
  if(inDataRight!=NULL)
    dumpMatrixDouble("fwrmdata", inDataRight, numExpected, numPoints*2);
#endif
  CkAssert(matrixA!=NULL);
  CkAssert(inDataLeft!=NULL);
  CkAssert(outData!=NULL);

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transformT, transform, m_in, n_in, k_in, alpha, beta, k_in, k_in, ldc);
#endif
  DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, matrixA, &lda, inDataLeft, &ldb, &beta, outData, &ldc);
  
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
#endif  // SPLIT
  if( (numRecd == numExpected * 2 )|| (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected))
    {
      numRecd = 0; 
    }

      
#ifdef _PAIRCALC_DEBUG_PARANOID_
  dumpMatrixDouble("fwgmodata",outData,grainSize, grainSize);
#endif


#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  // do the slicing and dicing to send bits to Ortho
  contributeSubTiles(outData);
#ifdef _CP_SUBSTEP_TIMING_
  double pstart=CmiWallTimer();
  if(symmetric)
    contribute(sizeof(double),&pstart,CkReduction::max_double, pairCalcID1.endTimerCB , pairCalcID1.forwardTimerID);

  else
    contribute(sizeof(double),&pstart,CkReduction::max_double, pairCalcID2.endTimerCB , pairCalcID2.forwardTimerID);


#endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
  if(phantomSym && symmetric && thisIndex.x!=thisIndex.y) //mirror our data to the phantom
    {
      CkAssert(existsRight);
      //      CkPrintf("[%d,%d,%d,%d,%d] fw sending phantom\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
      phantomMsg *phantom=new ( numExpected*numPoints*2, 8*sizeof(int)) phantomMsg;
      bool prioPhan=false;
      if(prioPhan)
	{
	  CkSetQueueing(phantom, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(phantom) = 1; // just make it slower
	  // than non prioritized
	}
      phantom->init(numExpected*numPoints*2, numPoints, flag_dp, inDataRight, blkSize);
      thisProxy(thisIndex.w,thisIndex.y, thisIndex.x,thisIndex.z).acceptPhantomData(phantom);
    }

  }


void
PairCalculator::sendTiles(bool flag_dp)
{
  int tilesq=orthoGrainSize*orthoGrainSize;
  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpIdOrtho).ckLocalBranch();
  int numOrtho=grainSize/orthoGrainSize;
  //  for(int orthoX=0; orthoX<numOrtho; orthoX++)
  //    for(int orthoY=0; orthoY<numOrtho; orthoY++)

  int progcounter=0;
  for(int orthoIndex=0;orthoIndex<numOrtho*numOrtho;orthoIndex++)
      {
	// copy into submatrix, contribute
	// we need to stride by grainSize and copy by orthoGrainSize
	//	int orthoIndex=orthoX*numOrtho+orthoY;

	if(touchedTiles[orthoIndex]==tilesq)
	  {
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	    if(symmetric)
	    CkPrintf("[%d %d %d %d %d]: contributes %i \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
#endif
	    //	    if (flag_dp) {
	    //	      for (int i = 0; i < orthoGrainSize*orthoGrainSize; i++)
	    //		outTiles[orthoIndex][i] *= 2.0;
	    //	    }
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	    int orthoX=orthoIndex/numOrtho;	   
	    int orthoY=orthoIndex%numOrtho;
	    char filename[80];
	    snprintf(filename,80,"fwoutTile_%d_%d:",orthoX,orthoY);
	    dumpMatrixDouble(filename, outTiles[orthoIndex], orthoGrainSize, orthoGrainSize,orthoX*orthoGrainSize, orthoY*orthoGrainSize);
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
	      CkPrintf("[%d %d %d %d %d]: %i not ready with %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, orthoIndex, touchedTiles[orthoIndex]);
#endif
	  }
      }
}


void
PairCalculator::contributeSubTiles(double *fullOutput)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d %d %d %d %d]: contributeSubTiles \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric);
#endif

  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpIdOrtho).ckLocalBranch();
  double *outTile=new double[orthoGrainSize*orthoGrainSize];
  bzero(outTile,sizeof(double)*orthoGrainSize*orthoGrainSize);
  //reuse the same tile each time as contribute makes its own copy
  int numOrtho=grainSize/orthoGrainSize;
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  dumpMatrixDouble("fullOutput", fullOutput, grainSize, grainSize);
#endif
  
  for(int orthoX=0; orthoX<numOrtho; orthoX++)
    for(int orthoY=0; orthoY<numOrtho; orthoY++)
      {
	// copy into submatrix, contribute
	// we need to stride by grainSize and copy by orthoGrainSize
	int orthoIndex=orthoX*numOrtho+orthoY;
	int tileStart=orthoX*orthoGrainSize*grainSize+orthoY*orthoGrainSize;
	CkAssert(orthoIndex<numOrtho*numOrtho);
	for(int ystart=0; ystart<orthoGrainSize*orthoGrainSize; ystart+=orthoGrainSize, tileStart+=grainSize)
	  for(int x=0;x<orthoGrainSize;x++)
	    {
	      outTile[ystart+x]=fullOutput[tileStart+x];
	    }
	
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
	char filename[80];
	snprintf(filename,80,"fwoutTile_%d_%d:",orthoX,orthoY);
	dumpMatrixDouble(filename, outTile, orthoGrainSize, orthoGrainSize,orthoX*orthoGrainSize, orthoY*orthoGrainSize);
#endif

	mcastGrp->contribute(orthoGrainSize*orthoGrainSize*sizeof(double), outTile, sumMatrixDoubleType, orthoCookies[orthoIndex], orthoCB[orthoIndex]);
	//mcastGrp->contribute(orthoGrainSize*orthoGrainSize*sizeof(double), outTile, CkReduction::sum_double, orthoCookies[orthoIndex], orthoCB[orthoIndex]);

      }
  delete [] outTile;
}


/**
 * Forward path multiply.  Accumulates rows and columns when all
 * arrive it calls the multiply
 */
void
PairCalculator::acceptPhantomData(phantomMsg *msg)
{
#ifdef _PAIRCALC_DEBUG_
  /  CkPrintf("[%d,%d,%d,%d,%d] accept Phantom\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
#endif


#ifdef _CP_SUBSTEP_TIMING_
  //satisfy fw path reduction for this triangle
  double pstart=CmiWallTimer();
  contribute(sizeof(double),&pstart,CkReduction::min_double, pairCalcID1.beginTimerCB , pairCalcID1.forwardTimerID);
  pstart=CmiWallTimer();
  contribute(sizeof(double),&pstart,CkReduction::max_double, pairCalcID1.endTimerCB , pairCalcID1.forwardTimerID);
#endif

  blkSize=msg->blkSize;
  CkAssert(amPhantom);
  if (!existsRight)
      { // now that we know numPoints we can allocate contiguous space
	CkAssert(inDataRight==NULL);
	existsRight=true;
	numPoints = msg->numPoints; // numPoints is init here with the size of the data chunk.
	inDataRight = new double[numExpected*numPoints*2];
	bzero(inDataRight,numExpected*numPoints*2*sizeof(double));
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d] Allocated right %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpected,numPoints);
#endif
      }
  CkAssert(numPoints== msg->numPoints); 
  memcpy(inDataRight, msg->points, numExpected*numPoints*2* sizeof(double));
  delete msg;
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
 */
void 
PairCalculator::multiplyPsiV()
{
  //reset for next time
  numRecd=0;

  // this is the same as the regular one matrix backward path
  // with the following exceptions

  // inDataLeft and inDataRight contain PsiV
  
  // outData contains the orthoT from the previous (standard) backward
  // path invocation

  // make a message
  multiplyResultMsg *psiV= new ( grainSize*grainSize,0,0 ) multiplyResultMsg;
  // copy in outData
  CkAssert(outData!=NULL);
  // this is slightly silly but we'll live with some memory inefficiency
  // in the infrequently called tolerance check
  psiV->init1(grainSize*grainSize, outData,0,0,PSIV);
  // convert this to an inline call once this all works
  thisProxy(thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z).multiplyResult(psiV);

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
  CkPrintf("[%d %d %d %d %d]: MultiplyResult with size %d numRecd %d actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size, numRecdBW, msg->actionType);
#endif
#ifdef _CP_SUBSTEP_TIMING_
  if(numRecdBW==0)
    {
      double pstart=CmiWallTimer();
      if(symmetric)
	contribute(sizeof(double),&pstart,CkReduction::min_double, pairCalcID1.beginTimerCB , pairCalcID1.backwardTimerID);

      else
	contribute(sizeof(double),&pstart,CkReduction::min_double, pairCalcID2.beginTimerCB , pairCalcID2.backwardTimerID);
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
  int numOrthoCol=grainSize/orthoGrainSize;


  // find our tile indices within this sGrain
  int orthoX=(msg->orthoX*orthoGrainSize-thisIndex.x)/orthoGrainSize;
  int orthoY=(msg->orthoY*orthoGrainSize-thisIndex.y)/orthoGrainSize;
  if(amPhantom)
    {
      orthoX=(msg->orthoX*orthoGrainSize-thisIndex.y)/orthoGrainSize;
      orthoY=(msg->orthoY*orthoGrainSize-thisIndex.x)/orthoGrainSize;

      //      CkPrintf("[%d %d %d %d %d]: phantom MultiplyResult with size %d numRecd %d actionType %d numPoints %d orthoX %d orthoY %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size, numRecdBW, msg->actionType, numPoints, orthoX, orthoY);
    }

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
  int matrixSize=grainSize*grainSize;

  //ASSUMING TMATRIX IS REAL (LOSS OF GENERALITY)

  double *amatrix=NULL;
  double *amatrix2=matrix2;  // may be overridden later

  // default DGEMM for non streaming comp case
  int m_in=numPoints*2;   // rows of op(A)==rows C	 
  int n_in=grainSize;     // columns of op(B)==columns C 
  int k_in=grainSize;     // columns op(A) == rows op(B) 
  double beta(0.0);

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

  if (orthoGrainSize==grainSize||actionType==PSIV)
    { // you were sent the correct section only
      amatrix=matrix1;
    }
  else if (collectAllTiles)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d %d %d %d %d]: MultiplyResult aggregating numRecdBW %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size, numRecdBW);
#endif      
      
      //do strided copy 
      // stride amatrix by grainSize,
      // stride matrix1 by orthoGrainSize
      // copy out the orthograinSize section
      // goffset=orthoX*orthoGrainSize+orthoY
      // into amatrix[goffset]
      // goffset+=grainSize
      if(numRecdBW==1) //alloc on first receipt
	{
	  CkAssert(inResult1==NULL);
	  inResult1 = new double[matrixSize];
	  bzero(inResult1,sizeof(double)*matrixSize);
	}
      if(!unitcoef){ // CG non minimization case have GAMMA      
	if(numRecdBW==1) //alloc on first receipt
	  {
	    inResult2 = new double[matrixSize];
	    bzero(inResult2,sizeof(double)*matrixSize);
	  }
	amatrix2 = inResult2;
	int tileStart=orthoX*orthoGrainSize*grainSize+orthoY*orthoGrainSize;
	if(symmetric && (thisIndex.x!=thisIndex.y)) //swap the non diagonals
	  tileStart=orthoY*orthoGrainSize*grainSize+orthoX*orthoGrainSize;
	for(int i=0; i<orthoGrainSize*orthoGrainSize; i+=orthoGrainSize,tileStart+=grainSize)
	  for(int j=0; j<orthoGrainSize; j++)
	    inResult2[tileStart+j] = matrix2[i+j];
      }
      amatrix = inResult1;
      int tileStart=orthoX*orthoGrainSize*grainSize+orthoY*orthoGrainSize;
      if(symmetric && (thisIndex.x!=thisIndex.y)) //swap the non diagonals
	tileStart=orthoY*orthoGrainSize*grainSize+orthoX*orthoGrainSize;
      for(int i=0; i<orthoGrainSize*orthoGrainSize; i+=orthoGrainSize,tileStart+=grainSize)
	for(int j=0; j<orthoGrainSize; j++)
	  inResult1[tileStart+j] = matrix1[i+j];
#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
      char filename[80];
      snprintf(filename,80,"bwinResult1_%d_%d:",orthoX,orthoY);
      dumpMatrixDouble(filename, matrix1, orthoGrainSize, orthoGrainSize,orthoX*orthoGrainSize, orthoY*orthoGrainSize);
#endif
    } //else !collectAllTiles
  else
    {  // settings for streaming computation on each tile
      amatrix=matrix1;
      if(!unitcoef)
	amatrix2=matrix2;
      // fix n_in and k_in 
      n_in=orthoGrainSize;
      k_in=orthoGrainSize;
      if(symmetric && (thisIndex.x!=thisIndex.y)) //swap the non diagonals
	{
	  int swap=orthoX;
	  orthoX=orthoY;
	  orthoY=swap;
	}
      BTCoffset=orthoY * m_in * orthoGrainSize;
      BTAoffset=orthoX * m_in * orthoGrainSize;

      BNCoffset=orthoX * m_in * orthoGrainSize;
      BNAoffset=orthoY * m_in * orthoGrainSize;

      beta=1.0;  // need to sum over tiles within orthoY columns
    }

  int numOrtho=numOrthoCol*numOrthoCol;
  if(orthoGrainSize==grainSize || numRecdBW==numOrtho || !collectAllTiles) // have all the input we need  
   {
     // call helper function to do the math
     // purely a readability improver

     bwMultiplyHelper(size, matrix1, matrix2, amatrix, amatrix2,  unitcoef, m_in, n_in, k_in, BNAoffset, BNCoffset, BTAoffset, BTCoffset, orthoX, orthoY, beta);
   }

  if(PCstreamBWout && !collectAllTiles && !useBWBarrier)  // send results which are complete and not yet sent
    {
      // just a readability improver

      bwSendHelper( orthoX, orthoY, numOrthoCol, numOrtho);
    }

  // check to see if we're all done
  if(((!PCstreamBWout && collectAllTiles) && (orthoGrainSize==grainSize || numRecdBW==numOrtho)) || (useBWBarrier && (orthoGrainSize==grainSize || numRecdBW==numOrtho))) 
    { // clean up
      if(collectAllTiles)
	{  // only safe to do this if we allocated them
	  if(inResult2!=NULL)
	    delete [] inResult2;
	  if(inResult1!=NULL)
	    delete [] inResult1;
	  inResult1=NULL;
	  inResult2=NULL;
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
	  else if(((!phantomSym && symmetric) || !unitcoef) && (thisIndex.x != thisIndex.y))
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
      if(conserveMemory)
	{
	  // clear the right and left they'll get reallocated on the next pass
	  delete [] inDataLeft;
	  inDataLeft=NULL;
	  if(!symmetric || (symmetric && thisIndex.x!=thisIndex.y )) {
	    delete [] inDataRight;
	    inDataRight = NULL;
	  }
	  existsLeft=false;
	  existsRight=false;
	  if(outData!=NULL && actionType!=KEEPORTHO)
	    {
	      delete [] outData;
	      outData = NULL;
	      existsOut=false;
	    }
	}
    }

}

void PairCalculator::bwMultiplyHelper(int size, double *matrix1, double *matrix2, double *amatrix, double *amatrix2, bool unitcoef, int m_in, int n_in, int k_in, int BNAoffset, int BNCoffset, int BTAoffset, int BTCoffset, int orthoX, int orthoY, double beta)
{

  if(symmetric && actionType==KEEPORTHO) // there will be a psiV step following
    {
      if(outData==NULL)
	{
	  CkAssert(!existsOut);
	  outData=new double[grainSize*grainSize];
	  bzero(outData,sizeof(double)*grainSize*grainSize);
	  existsOut=true;
	}
      CkAssert(size==grainSize*grainSize);
      //keep the orthoT we just received in matrix1
      if(!collectAllTiles && orthoGrainSize!=grainSize)
	{ // copy this tile
	  int tileStart=orthoX*orthoGrainSize*grainSize+orthoY*orthoGrainSize;
	  if(symmetric && (thisIndex.x!=thisIndex.y)) //swap the non diagonals
	    tileStart=orthoY*orthoGrainSize*grainSize+orthoX*orthoGrainSize;
	  for(int i=0; i<orthoGrainSize*orthoGrainSize; i+=orthoGrainSize,tileStart+=grainSize)
	    for(int j=0; j<orthoGrainSize; j++)
	      outData[tileStart+j] = matrix1[i+j];
	}
      else
	{
	  memcpy(outData, amatrix, size*sizeof(double));
	}
      // it is safe to reuse this memory 
      // normal backward path has no use for outData
      // forward path won't be called again until after we're done with outData.
    }
  if(!amPhantom && mynewData==NULL)
    {
      CkAssert(numPoints>0);
      mynewData = new complex[numPoints*grainSize];
      existsNew=true;
      bzero(mynewData,numPoints*grainSize* sizeof(complex));
    }
  if(amPhantom || ((symmetric || !unitcoef) && (thisIndex.x != thisIndex.y))){
    CkAssert(numPoints>0);
    if(othernewData==NULL)
      {
	othernewData = new complex[numPoints*grainSize];
	bzero(othernewData,numPoints*grainSize* sizeof(complex));
      }
  }
  else
    othernewData=NULL;

  double *mynewDatad= reinterpret_cast <double *> (mynewData);
  double alpha(1.0);
  char transform='N';
  char transformT='T';

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_

  int chunksize=blkSize/numChunks;
  int ystart=chunksize*thisIndex.z;
  if(!amPhantom)
    dumpMatrixDouble("bwmlodata",inDataLeft,numExpected,numPoints*2,0,ystart);
  if(!unitcoef){ // CG non minimization case
    dumpMatrixDouble("bwmrodata",inDataRight,numExpected,numPoints*2,0,ystart);
  }
  if(grainSize==orthoGrainSize)
    {
      dumpMatrixDouble("bwm1idata",matrix1,grainSize,grainSize);
      if(!unitcoef)
	{ // CG non minimization case
	  dumpMatrixDouble("bwm2idata",matrix2,grainSize,grainSize);
	}
    }
  else
    {
      dumpMatrixDouble("bwm1idata",inResult1,grainSize,grainSize);
      if(!unitcoef)
	{ // CG non minimization case
	  dumpMatrixDouble("bwm2idata",inResult2,grainSize,grainSize);
	}
    }
#endif

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

  //first multiply to apply the T or L matrix
  if(!amPhantom)
    {
#if PC_BWD_DGEMM_SPLIT > 0
      if(symmetric)
	dgemmSplitBwdM(m_in, n_in, k_in, &transform, &transform, &alpha, &(inDataLeft[BNAoffset]),  amatrix, &beta, &(mynewDatad[BNCoffset]));
      else
	dgemmSplitBwdM(m_in, n_in, k_in, &transform, &transformT, &alpha, &(inDataLeft[BTAoffset]),  amatrix, &beta, &(mynewDatad[BTCoffset]));
#else // no split

#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(inDataLeft[BNAoffset] )%16==0);
      CkAssert((unsigned int) amatrix %16==0);
      CkAssert((unsigned int)&(mynewDatad[BNCoffset] )%16==0);
#endif

      if(symmetric) {
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transform, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif
	DGEMM(&transform, &transform, &m_in, &n_in, &k_in, &alpha, &(inDataLeft[BNAoffset]), &m_in, amatrix, &k_in, &beta, &(mynewDatad[BNCoffset]), &m_in);
      }
      else {
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transformT, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif
	DGEMM(&transform, &transformT, &m_in, &n_in, &k_in, &alpha, &(inDataLeft[BTAoffset]), &m_in,  amatrix, &k_in, &beta, &(mynewDatad[BTCoffset]), &m_in);
      }
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(230, StartTime, CmiWallTimer());
#endif
#endif // end of split
    }

#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  if((amPhantom) || (!phantomSym && symmetric && (thisIndex.x !=thisIndex.y)) && existsRight)
    {
      double *othernewDatad= reinterpret_cast <double *> (othernewData);

#if PC_BWD_DGEMM_SPLIT > 0
      dgemmSplitBwdM(m_in, n_in, k_in, &transform, &transformT, &alpha, &(inDataRight[BTAoffset]),  amatrix, &beta, &(othernewDatad[BTCoffset]));
#else  // no split

      CmiNetworkProgress();
#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(inDataRight[BTAoffset] )%16==0);
      CkAssert((unsigned int) amatrix %16==0);
      CkAssert((unsigned int)&(othernewDatad[BTCoffset] )%16==0);
#endif
 		 
#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transformT, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif
      DGEMM(&transform, &transformT, &m_in, &n_in, &k_in, &alpha, &(inDataRight[BTAoffset]), &m_in,  amatrix, &k_in, &beta, &(othernewDatad[BTCoffset]), &m_in);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(250, StartTime, CmiWallTimer());
#endif
#endif  // end of split
    }

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
  if(!amPhantom)
    dumpMatrixComplex("bwgmodata",mynewData,grainSize,numPoints,0,ystart);
#endif


  if(!unitcoef){ // CG non minimization case  GAMMA

    // output modified by subtracting an application of orthoT

    // TODO: figure out if this can be streamed, the beta=0 on off
    // diagonal looks like a sticky correctness killer

    // for now, force no streaming "all matrixes arrived" condition
    //	CkAssert(collectAllTiles);
    // C = alpha*A*B + beta*C
    // C= -1 * inRight * orthoT + C
    double *othernewDatad;
#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif
    alpha=-1.0;  //comes in with a minus sign
    if(thisIndex.x!=thisIndex.y){
      beta=0.0; // new contribution off-diagonal
      othernewDatad= reinterpret_cast <double *> (othernewData);    
    }else{
      beta=1.0; //subtract contribution from existing on diagonal
      othernewDatad=mynewDatad;
    }//endif

    CmiNetworkProgress();

#ifdef PRINT_DGEMM_PARAMS
  CkPrintf("HEY-DGEMM %c %c %d %d %d %f %f %d %d %d\n", transform, transform, m_in, n_in, k_in, alpha, beta, m_in, k_in, m_in);
#endif
    DGEMM(&transform, &transform, &m_in, &n_in, &k_in, &alpha, &(inDataRight[BNAoffset]), 
	  &m_in, amatrix2, &k_in, &beta, &(othernewDatad[BNCoffset]), &m_in);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(240, StartTime, CmiWallTimer());
#endif

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
    dumpMatrixComplex("bwg2modata",othernewData,grainSize,numPoints);
#endif

  } // end  CG case

#ifdef _PAIRCALC_VALID_OUT_
  CkPrintf("[PAIRCALC] [%d %d %d %d %d] backward gemm out %.10g %.10g %.10g %.10g \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric, mynewDatad[0],mynewDatad[1],mynewData[numPoints*grainSize-1].re,mynewData[numPoints*grainSize-1].im);
#endif

}

void PairCalculator::bwSendHelper(int orthoX, int orthoY, int numOrthoCol, int numOrtho)
{

  // not supported in dynamics until we figure out how to stream that
  // computation.
      
  // Check to see if this column is  complete, if so send it.
      
      
  // symmetric uses BNCoffset on newData for diagonal.      
  // symmetric also BTCoffset on otherNewData for off diagonal.
    // asymmetric uses BTC offset on newData
  int orthoXgrain=orthoX*orthoGrainSize;
  int orthoYgrain=orthoY*orthoGrainSize;
  if(symmetric)
    {
      if(!amPhantom)
	{
	  if(++columnCount[orthoX]==numOrthoCol ) // BNC
	    {
	      // send orthoX in newData
	      if(gSpaceSum)
		sendBWResultColumnDirect(false, orthoXgrain, orthoXgrain+orthoGrainSize);
	      else
		sendBWResultColumn(false, orthoXgrain, orthoXgrain+orthoGrainSize);
	    }
	}
      CkAssert(columnCount[orthoY]<=numOrthoCol);
    }
  else //asymm
    {
      if(++columnCount[orthoY]==numOrthoCol) // BTC 
	{
	  //send orthoY in newdata
	  if(gSpaceSum)
	    sendBWResultColumnDirect(false, orthoYgrain, orthoYgrain+orthoGrainSize);
	  else
	    sendBWResultColumn(false, orthoYgrain, orthoYgrain+orthoGrainSize);
	}
      CkAssert(columnCount[orthoY]<=numOrthoCol);
    }
  if((amPhantom) || (!phantomSym && symmetric && (thisIndex.x !=thisIndex.y)) && existsRight)
    {
      if(++columnCountOther[orthoY]==numOrthoCol) // BTC
	{
	  //send orthoY in otherNewData
	  if(gSpaceSum)
	    sendBWResultColumnDirect(true, orthoYgrain, orthoYgrain+orthoGrainSize);
	  else
	    sendBWResultColumn(true, orthoYgrain, orthoYgrain+orthoGrainSize);
	}
      CkAssert(columnCountOther[orthoY]<=numOrthoCol);
    }

  // this could be refined to track an array of completed columns
  // and send them in some grouping scheme
  CkAssert(numRecdBW<=numOrtho);
  if(numRecdBW==numOrtho)
    { //all done clean up after ourselves
      //      CkPrintf("[PAIRCALC] [%d %d %d %d %d] cleaning up \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric);
      if(!amPhantom)
	{
	  delete [] mynewData;
	  mynewData=NULL;
	  existsNew=false;
	}
      if(othernewData!=NULL)
	delete [] othernewData;
      othernewData=NULL;
      bzero(columnCount, sizeof(int) * numOrthoCol);
      bzero(columnCountOther, sizeof(int) * numOrthoCol);
      numRecdBW=0;
#ifdef _CP_SUBSTEP_TIMING_
      double pstart=CmiWallTimer();
      if(symmetric)
	contribute(sizeof(double),&pstart,CkReduction::max_double, pairCalcID1.endTimerCB , pairCalcID1.backwardTimerID);
      else
	contribute(sizeof(double),&pstart,CkReduction::max_double, pairCalcID2.endTimerCB , pairCalcID2.backwardTimerID);


#endif

    }
}

void
PairCalculator::sendBWResultColumnDirect(bool otherdata, int startGrain, int endGrain  )
{
#ifdef _PAIRCALC_DEBUG_CONTRIB_
  CkPrintf("[%d %d %d %d %d]: sendBWResultColumnDirect with actionType %d startGrain %d sendGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType, startGrain, endGrain);
#endif


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
	complex *computed=&(othernewData[j*numPoints]);
	CkCallback mycb(cp_entry, CkArrayIndex2D(j+ index ,thisIndex.w), cb_aid);
	partialResultMsg *msg=new (numPoints, 8*sizeof(int) ) partialResultMsg;
	if(gpriority)
	  {
	    *((int*)CkPriorityPtr(msg)) = gpriority;
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  }
	msg->init(numPoints, thisIndex.z, computed);

#ifdef _PAIRCALC_DEBUG_
	CkPrintf("sending partial of size %d offset %d to [%d %d]\n",numPoints,j,thisIndex.y+j,thisIndex.w);
#endif

#ifdef _NAN_CHECK_
	for(int i=0;i<msg->N ;i++)
	  {
	    if(!finite(msg->result[i].re)|| !finite(msg->result[i].im))
	      {
		CkPrintf("[%d %d %d %d %d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
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
	complex *computed=&(mynewData[j*numPoints]);
#ifdef _NAN_CHECK_
	for(int i=0;i<numPoints ;i++)
	  {
	    if(!finite(computed[i].re) || !finite(computed[i].im))
	      {
		fprintf(stderr,"[%d %d %d %d %d]: sendBWResultColumnDirect nan in computed at %d %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, j,i);
		CkAssert(finite(computed[i].re));
		CkAssert(finite(computed[i].im));
	      }
	  }
#endif

	CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);
	partialResultMsg *msg=new (numPoints, 8*sizeof(int) ) partialResultMsg;
	if(gpriority)
	  {
	    *((int*)CkPriorityPtr(msg)) = gpriority;
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  }
	msg->init(numPoints, thisIndex.z, computed);
	/*
	  msg->N=numPoints;
	  msg->myoffset = thisIndex.z; // chunkth
	  memcpy(msg->result,mynewData+j*numPoints,msg->N*sizeof(complex));
	*/
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("sending partial of size %d offset %d to [%d %d]\n",numPoints,j,thisIndex.y+j,thisIndex.w);
#endif

#ifdef _NAN_CHECK_
	for(int i=0;i<msg->N ;i++)
	  {
	    if(!finite(msg->result[i].re)|| !finite(msg->result[i].im))
	      {
		CkPrintf("[%d %d %d %d %d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
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

void
PairCalculator::sendBWResultColumn(bool otherdata, int startGrain, int endGrain  )
{


#ifdef _PAIRCALC_DEBUG_CONTRIB_
  CkPrintf("[%d %d %d %d %d]: sendBWResultColumn with actionType %d startGrain %d sendGrain %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType, startGrain, endGrain);
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
	CkPrintf("[%d %d %d %d %d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
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
      for(int j=startGrain;j<endGrain;j++) //mynewdata
	{
	complex *computed=&(mynewData[j*numPoints]);
#ifndef CMK_OPTIMIZE
	double StartTime=CmiWallTimer();
#endif
	CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);

#ifdef _PAIRCALC_DEBUG_CONTRIB_
	  CkPrintf("[%d %d %d %d %d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
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

void
PairCalculator::sendBWResultDirect(sendBWsignalMsg *msg)
{

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d %d %d %d %d]: sendBWResultDirect with actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType);
#endif
  bool otherdata=msg->otherdata;
  delete msg;  
  // Now we have results in mynewData and if(symmetric) othernewData
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
    int index=thisIndex.x;
    if(amPhantom)
      index=thisIndex.y;
    for(int j=0;j<grainSize;j++)
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
	omsg->init(numPoints, thisIndex.z, computed);
	/*	omsg->N=numPoints;
		omsg->myoffset = thisIndex.z; // chunkth
		memcpy(omsg->result,othernewData+j*numPoints,omsg->N*sizeof(complex));
	*/
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("sending partial of size %d offset %d to [%d %d]\n",numPoints,j,index+j,thisIndex.w);
#endif
#ifdef _NAN_CHECK_
	for(int i=0;i<omsg->N ;i++)
	  {
	    if(!finite(omsg->result[i].re) || !finite(omsg->result[i].im))
	      {
		CkPrintf("[%d %d %d %d %d]: sendBWResultDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
	      }
	    CkAssert(finite(omsg->result[i].re));
	    CkAssert(finite(omsg->result[i].im));
	  }
#endif

	mycb.send(omsg);

      }
    if(otherdata)
      delete [] othernewData;
    othernewData=NULL;
    existsNew=false;
  }
  if(!amPhantom)
    {
      CkAssert(mynewData!=NULL);
      for(int j=0;j<grainSize;j++) //mynewdata
	{
	  complex *computed=&(mynewData[j*numPoints]);
	  CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);
	  partialResultMsg *omsg= new (numPoints, 8*sizeof(int) ) partialResultMsg;
	  if(gpriority)
	    {
	      *((int*)CkPriorityPtr(omsg)) = gpriority;
	      CkSetQueueing(omsg, CK_QUEUEING_IFIFO);
	    }
	  omsg->init(numPoints, thisIndex.z, computed);
	  /*
	    omsg->N=numPoints;
	    omsg->myoffset = thisIndex.z; // chunkth
	    memcpy(omsg->result, mynewData+j*numPoints, omsg->N*sizeof(complex) );
	  */
#ifdef _PAIRCALC_DEBUG_
	  CkPrintf("sending partial of size %d offset %d to [%d %d]\n",numPoints,j,thisIndex.y+j,thisIndex.w);
#endif
#ifdef _NAN_CHECK_
	for(int i=0;i<omsg->N ;i++)
	  {
	    if(!finite(omsg->result[i].re) || !finite(omsg->result[i].im))
	      {
		CkPrintf("[%d %d %d %d %d]: sendBWResultColumnDirect nan at %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, i);
	      }
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
/**
 * Send the results via multiple reductions as triggered by a prioritized message
 */
void
PairCalculator::sendBWResult(sendBWsignalMsg *msg)
{

#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d %d %d %d %d]: sendBWResult with actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, actionType);
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
    for(int j=0;j<grainSize;j++)
      {
	//this callback creation could be obviated by keeping an
	//array of callbacks, not clearly worth doing
	complex *computed=&(othernewData[j*numPoints]);
	CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.x ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	CkPrintf("[%d %d %d %d %d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
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
      CkAssert(mynewData!=NULL);
      for(int j=0;j<grainSize;j++) //mynewdata
	{
	complex *computed=&(mynewData[j*numPoints]);
	  CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	  CkPrintf("[%d %d %d %d %d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
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




void PairCalculator::dumpMatrixDouble(const char *infilename, double *matrix, int xdim, int ydim, int xstart, int ystart)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename,fmt, thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
  FILE *loutfile = fopen(filename, "w");
  for(int i=0;i<xdim;i++)
    for(int j=0;j<ydim;j++)
      fprintf(loutfile,"%d %d %.12g\n",i+xstart,j+ystart,matrix[i*ydim+j]);
  fclose(loutfile);
}

void PairCalculator::dumpMatrixComplex(const char *infilename, complex *matrix, int xdim, int ydim, int xstart, int ystart)
{
  char fmt[1000];
  char filename[1000];
  strncpy(fmt,infilename,999);
  strncat(fmt,"_%d_%d_%d_%d_%d.out",999);
  sprintf(filename, fmt, thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
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
	int numOrthoCol=grainSize/orthoGrainSize;
	int numOrtho=numOrthoCol*numOrthoCol;
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
        DGEMM(trans, transT, &Msplit, &n, &k, alpha, A, &m, B, &k, bt, C, &m);
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
        DGEMM(trans, transT, &Msplit, &n, &k, alpha, &A[off], &m, B, &k, bt, &C[off], &m);

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

#include "ckPairCalculator.def.h"


