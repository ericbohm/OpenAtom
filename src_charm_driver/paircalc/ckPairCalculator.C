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
 * The elan code is a specialized machine reduction/broadcast which runs   
 * much faster on lemieux, and theoretically on any other elan machine     
 * It operates by one dummy reduction which is used to indicate that all   
 * calculators on a PE machine have reported in.  Then the machine         
 * reduction is triggered.                                                 
 *                                                                         
 * The paircalculator is a 4 dimensional array.  Those dimensions are:     
 *            w: gspace state plane (the second index of the 2D gspace)    
 *            x: coordinate offset within plane (a factor of grainsize)    
 *            y: coordinate offset within plane (a factor of grainsize)    
 *            z: blocksize, unused always 0                                
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
 * true.  These will be recored in the left data array.  We will then
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
 * together in its reduction client.  A functionality that is probably
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
 * matrix of points.
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

#define PARTITION_SIZE 500

ComlibInstanceHandle mcastInstanceCP;

CkReduction::reducerType sumMatrixDoubleType;

void registersumMatrixDouble(void)
{ 
  sumMatrixDoubleType=CkReduction::addReducer(sumMatrixDouble);
}


// sum together matrices of doubles
// possibly faster than sum_double due to minimizing copies
inline CkReductionMsg *sumMatrixDouble(int nMsg, CkReductionMsg **msgs)
{
  double *ret=(double *)msgs[0]->getData();

  int size0=msgs[0]->getSize();
  int size=size0/sizeof(double);

  double *inmatrix;
  for(int i=1; i<nMsg;i++)
    {
      inmatrix=(double *) msgs[i]->getData();
      for(int d=0;d<size;d++)
	ret[d]+=inmatrix[d];
    }
  return CkReductionMsg::buildNew(size*sizeof(double),ret);
}


PairCalculator::PairCalculator(CkMigrateMessage *m) { }

PairCalculator::PairCalculator(bool sym, int grainSize, int s, int numChunks, CkCallback cb, CkArrayID cb_aid, int _cb_ep, int _cb_ep_tol, bool conserveMemory, bool lbpaircalc,  redtypes _cpreduce, int _orthoGrainSize)
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
  orthoGrainSize=_orthoGrainSize;
  existsLeft=false;
  existsRight=false;
  existsOut=false;
  existsNew=false;
  numRecd = 0;
  numExpected = grainSize;
  cpreduce=_cpreduce;
  resumed=true;

  inDataLeft = NULL;
  inDataRight = NULL;

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
  resultCookies=new CkSectionInfo[grainSize];
  int numOrtho=grainSize/orthoGrainSize;
  numOrtho=numOrtho*numOrtho;
  orthoCookies=new CkSectionInfo[numOrtho];
  orthoCB=new CkCallback[numOrtho];
  if(thisIndex.x != thisIndex.y)  
    // we don't actually use these in the asymmetric minimization case
    otherResultCookies=new CkSectionInfo[grainSize];

}

void
PairCalculator::pup(PUP::er &p)
{
  ArrayElement4D::pup(p);
  p|numRecd;
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
  p|resumed;
  p|rck;
  p|actionType;
  p|orthoGrainSize;
  int numOrtho=grainSize/orthoGrainSize;
  numOrtho=numOrtho*numOrtho;
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
  mCastGrpId=msg->mCastGrpId; //redundant
  /*  cpreduce=section;
  if(msg->lbsync)
  {
      int foo=1;
      contribute(sizeof(int), &foo , CkReduction::sum_int, msg->synccb);
  }
  */
  //    delete msg; do not delete nokeep method
}

//!initialize the section cookie for each slice of the result
void PairCalculator::initResultSection(initResultMsg *msg)
{

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

  mCastGrpId=msg->mCastGrpId;  //arguably redundant
  //to force synchronize in lb resumption
  if(msg->lbsync && rck==grainSize)
  {
      contribute(sizeof(int), &rck, CkReduction::sum_int, msg->synccb);
  }
  //    delete msg; do not delete nokeep method
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
#ifdef _PAIRCALC_DEBUG_
  CkPrintf(" symm=%d    pairCalc[%d %d %d %d] got from [%d %d] with size {%d}, from=%d, count=%d, resumed=%d\n", symmetric, thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,  thisIndex.w, msg->sender, msg->size, msg->fromRow, numRecd,resumed);
#endif
  if(!resumed)
    {
      ResumeFromSync();
    }
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
    CkAssert(numPoints==msg->size);
    CkAssert(offset<numExpected);
    memcpy(&(inDataLeft[offset*numPoints*2]), msg->points, numPoints * 2 *sizeof(double));
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Copying into offset*numPoints %d * %d numPoints *2 %d points start %.12g end %.12g\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z, symmetric, offset, numPoints, numPoints*2,msg->points[0].re, msg->points[numPoints-1].im);
#endif
#ifdef _PAIRCALC_DEBUG_STUPID_PARANOID_
    CkPrintf("copying in \n");
    double re;
    double im;
    for(int i=0;i<numPoints;i++)
      {
	re=msg->points[i].re;
	im=msg->points[i].im;
	CkPrintf("CL %d %g %g\n",i,re,im);
	if(fabs(re)>0.0)
	  CkAssert(fabs(re)>1.0e-300);
	if(fabs(im)>0.0)
	  CkAssert(fabs(im)>1.0e-300);
      }
#endif
  }
  else {
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
    memcpy(&(inDataRight[offset*numPoints*2]), msg->points, numPoints * 2 *sizeof(double));
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Copying into offset*numPoints %d * %d numPoints *2 %d points start %.12g end %.12g\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,offset,numPoints, numPoints*2,msg->points[0].re, msg->points[numPoints-1].im);
#endif
#ifdef _PAIRCALC_DEBUG_STUPID_PARANOID_
    CkPrintf("copying in \n");
    double re;
    double im;
    for(int i=0;i<numPoints;i++)
      {
	re=msg->points[i].re;
	im=msg->points[i].im;
	CkPrintf("CR %d %g %g\n",i,re,im);
	if(fabs(re)>0.0)
	  CkAssert(fabs(re)>1.0e-300);
	if(fabs(im)>0.0)
	  CkAssert(fabs(im)>1.0e-300);
      }
#endif
  }


  /*
   *  NOTE: For this to work the data chunks of the same plane across
   *  all states must be of the same size
   */
  // copy the input into our matrix
  /*
   * Once we have accumulated all rows  we gemm it.
   */

  if (numRecd == numExpected * 2 || (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected)) 
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
  int lda=doubleN;   //leading dimension A
  int ldb=doubleN;   //leading dimension B
  int ldc=numExpected;   //leading dimension C
  double alpha=double(1.0);//multiplicative identity
  double beta=double(0.0); // C is unset
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] gemming %c %c %d %d %d %f A %d B %d %f C %d\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,transformT,transform, m_in, n_in, k_in, alpha, lda, ldb, beta, ldc);
#endif
  if( numRecd == numExpected * 2)
    {
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
      dumpMatrixDouble("fwlmdata", inDataLeft, numExpected, numPoints*2);
      dumpMatrixDouble("fwrmdata", inDataRight, numExpected, numPoints*2);
#endif
      CkAssert(inDataRight!=NULL);
      CkAssert(inDataLeft!=NULL);
      CkAssert(outData!=NULL);
      DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, inDataRight, &lda, inDataLeft, &ldb, &beta, outData, &ldc);

    }
  else if (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected)
    {
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
      dumpMatrixDouble("fwlmdata",inDataLeft,numExpected, numPoints*2);
#endif
      DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, inDataLeft, &lda, inDataLeft, &ldb, &beta, outData, &ldc);
    }

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif

  numRecd = 0; 
  if (flag_dp) {
    for (int i = 0; i < grainSize*grainSize; i++)
      outData[i] *= 2.0;
  }

#ifdef _PAIRCALC_DEBUG_PARANOID_
  dumpMatrixDouble("fwgmodata",outData,grainSize, grainSize);
#endif


#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  // do the slicing and dicing to send bits to Ortho
  contributeSubTiles(outData);
  //  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  // mcastGrp->contribute(grainSize*grainSize*sizeof(double), outData, sumMatrixDoubleType, cookie, cb);
  //mcastGrp->contribute(grainSize*grainSize*sizeof(double), outData, CkReduction::sum_double, cookie, cb);


#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
}

void
PairCalculator::contributeSubTiles(double *fullOutput)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d %d %d %d %d]: contributeSubTiles \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric);
#endif

  CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  double *outTile=new double[orthoGrainSize*orthoGrainSize];
  //reuse the same tile each time as contribute makes its own copy
  int numOrtho=grainSize/orthoGrainSize;
#ifdef _PAIRCALC_DEBUG_PARANOID_FW_
  dumpMatrixDouble("fullOutput", fullOutput, grainSize, grainSize);
#endif
  //tranpose the output first
  
  for(int orthoX=0; orthoX<numOrtho; orthoX++)
    for(int orthoY=0; orthoY<numOrtho; orthoY++)
      {
	// copy into submatrix, contribute
	// we need to stride by sGrainSize and copy by orthoGrainSize
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
	/* cheap dirty transpose
	  double *dest= outTile;
	int chunksize=orthoGrainSize;
	double tmp;
	for(int i = 0; i < chunksize; i++)
	  for(int j = i + 1; j < chunksize; j++){
	    tmp = dest[i * chunksize + j];
	    dest[i * chunksize + j] = dest[j*chunksize + i];
	    dest[j * chunksize + i] = tmp;
	  }
	*/

	mcastGrp->contribute(orthoGrainSize*orthoGrainSize*sizeof(double), outTile, sumMatrixDoubleType, orthoCookies[orthoIndex], orthoCB[orthoIndex]);



      }
  delete [] outTile;
}

//PairCalculator::multiplyResult(int size, double *matrix1, double *matrix2)
void
PairCalculator::multiplyResultI(multiplyResultMsg *msg)
{
  
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
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d %d %d %d %d]: MultiplyResult with size %d numRecd %d actionType %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size, numRecd, msg->actionType);
#endif
  numRecd++; 
  int size=msg->size;
  int size2=msg->size2;
  double *matrix1=msg->matrix1;
  double *matrix2=msg->matrix2;
  actionType=msg->actionType;
  bool unitcoef = false;
  int numOrtho=grainSize/orthoGrainSize;
  int orthoX=(msg->orthoX*orthoGrainSize-thisIndex.x)/orthoGrainSize;
  int orthoY=(msg->orthoY*orthoGrainSize-thisIndex.y)/orthoGrainSize;
  int orthoIndex=orthoY*numOrtho+orthoX;
  if(matrix2==NULL||size2<1) 
    {
      unitcoef = true;
    }

  int offset = 0, index = thisIndex.y*numStates + thisIndex.x;

  if(!symmetric)
    index = thisIndex.x*numStates + thisIndex.y;
  int matrixSize=grainSize*grainSize;
  //ASSUMING TMATRIX IS REAL (LOSS OF GENERALITY)
  register double m=0;
  bool makeLocalCopy=false;
  double *amatrix=NULL;
  double *amatrix2=matrix2;  // may be overridden later
  
  if(numStates==grainSize)// all at once no malloc
    {
      amatrix=matrix1+index;  // index is 0 in this case, so this is silly
    }
  else if (orthoGrainSize==grainSize||actionType==PSIV)
    { // you were sent the correct section only
      amatrix=matrix1;
    }
  else
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("[%d %d %d %d %d]: MultiplyResult aggregating numRecd %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size, numRecd);
#endif      
      //do strided copy 
      // stride amatrix by grainSize,
      // stride matrix1 by orthoGrainSize
      // copy out the orthograinSize section
      // goffset=orthoX*orthoGrainSize+orthoY
      // into amatrix[goffset]
      // goffset+=grainSize
      if(numRecd==1) //alloc on first receipt
	{
	  CkAssert(inResult1==NULL);
	  inResult1 = new double[matrixSize];
	}
      if(!unitcoef){ // CG non minimization case have GAMMA      
	if(numRecd==1) //alloc on first receipt
	  inResult2 = new double[matrixSize];
	amatrix2 = inResult2;
	int tileStart=orthoX*orthoGrainSize*grainSize+orthoY*orthoGrainSize;
	for(int i=0; i<orthoGrainSize*orthoGrainSize; i+=orthoGrainSize,tileStart+=grainSize)
	  for(int j=0; j<orthoGrainSize; j++)
	    inResult2[tileStart+j] = matrix2[i+j];
	//slightly evil, but safe
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
    } //else

  numOrtho=numOrtho*numOrtho;
  if(orthoGrainSize==grainSize || numRecd==numOrtho)
    { // we have all the result matrices we need
      // if(cpreduce==section) //could update cookie from multicast
      // we rely on manual resets for this in migration case
      //      CkGetSectionInfo(cookie,msg);
      if(symmetric && actionType==KEEPORTHO) // there will be a psiV step following
	{
	  if(outData==NULL)
	    {
	      CkAssert(!existsOut);
	      outData=new double[grainSize*grainSize];
	      existsOut=true;
	    }
	  CkAssert(size==grainSize*grainSize);
	  //keep the orthoT we just received in matrix1
	  memcpy(outData, amatrix, size*sizeof(double));
	  // it is safe to reuse this memory 
	  // normal backward path has no use for outData
	  // forward path won't be called again until after we're done with outData.
	}
      CkAssert(mynewData==NULL);
      mynewData = new complex[numPoints*grainSize];
      existsNew=true;
      bzero(mynewData,numPoints*grainSize* sizeof(complex));
      if((symmetric || !unitcoef) && (thisIndex.x != thisIndex.y)){
	if(othernewData==NULL)
	  othernewData = new complex[numPoints*grainSize];
	bzero(othernewData,numPoints*grainSize* sizeof(complex));
      }
      else
	othernewData=NULL;

      double *localMatrix;
      double *outMatrix;

      int m_in=grainSize;
      int n_in=numPoints*2;
      int k_in=grainSize;
      double *mynewDatad= reinterpret_cast <double *> (mynewData);
      double alphad(1.0);
      double betad(0.0);
      char transform='N';
      char transformT='T';

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_

      int chunksize=blkSize/numChunks;
      int ystart=chunksize*thisIndex.z;

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
      if(symmetric)
	DGEMM(&transform, &transform, &n_in, &m_in, &k_in, &alphad, inDataLeft, &n_in,  amatrix, &k_in, &betad, mynewDatad, &n_in);
      else
	DGEMM(&transform, &transformT, &n_in, &m_in, &k_in, &alphad, inDataLeft, &n_in,  amatrix, &k_in, &betad, mynewDatad, &n_in);

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(230, StartTime, CmiWallTimer());
#endif
      if(symmetric && (thisIndex.x !=thisIndex.y) && existsRight)
	{
	  //#ifdef CMK_VERSION_BLUEGENE
	  CmiNetworkProgressAfter(1);
	  //#endif

#ifndef CMK_OPTIMIZE
	  StartTime=CmiWallTimer();
#endif
	  double *othernewDatad= reinterpret_cast <double *> (othernewData);
	  DGEMM(&transform, &transformT, &n_in, &m_in, &k_in, &alphad, inDataRight, &n_in,  amatrix, &k_in, &betad, othernewDatad, &n_in);
#ifndef CMK_OPTIMIZE
	  traceUserBracketEvent(250, StartTime, CmiWallTimer());
#endif
	}


#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
      dumpMatrixComplex("bwgmodata",mynewData,grainSize,numPoints,0,ystart);
#endif


      if(!unitcoef){ // CG non minimization case  GAMMA
	// C = alpha*A*B + beta*C
	// C= -1 * inRight * orthoT + C
	double *othernewDatad;
#ifndef CMK_OPTIMIZE
	StartTime=CmiWallTimer();
#endif
	alphad=-1.0;  //comes in with a minus sign
	if(thisIndex.x!=thisIndex.y){
	  betad=0.0; // new contribution off-diagonal
	  othernewDatad= reinterpret_cast <double *> (othernewData);    
	}else{
	  betad=1.0; //subtract contribution from existing on diagonal
	  othernewDatad=mynewDatad;
	}//endif

	DGEMM(&transform, &transform, &n_in, &m_in, &k_in, &alphad, inDataRight, 
	      &n_in,  amatrix2, &k_in, &betad, othernewDatad, &n_in);

#ifdef _PAIRCALC_DEBUG_PARANOID_BW_
	dumpMatrixComplex("bwg2modata",othernewData,grainSize,numPoints);
#endif

#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(240, StartTime, CmiWallTimer());
#endif
      } // end  CG case

#ifdef _PAIRCALC_VALID_OUT_
      CkPrintf("[PAIRCALC] [%d %d %d %d %d] backward gemm out %.10g %.10g %.10g %.10g \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric, mynewDatad[0],mynewDatad[1],mynewData[numPoints*grainSize-1].re,mynewData[numPoints*grainSize-1].im);
#endif

      sendBWsignalMsg *sigmsg=new (8*sizeof(int)) sendBWsignalMsg;
      //collapse this into 1 flag
      sigmsg->otherdata= 
	((!unitcoef || symmetric) &&(thisIndex.x !=thisIndex.y)) ? true : false;

      CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
      *(int*)CkPriorityPtr(sigmsg) = 1; // just make it slower than non prioritized
      thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);
      //sendBWResult(sigmsg);
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
      if(inResult2!=NULL)
	delete [] inResult2;
      if(inResult1!=NULL)
	delete [] inResult1;
      inResult1=NULL;
      inResult2=NULL;
      numRecd=0;
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
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
  if(otherdata){  // we have this othernewdata issue for the symmetric case
    // and the asymmetric dynamic case
    // for the off diagonal elements
    CkAssert(othernewData!=NULL);
    for(int j=0;j<grainSize;j++)
      {
	//this callback creation could be obviated by keeping an
	//array of callbacks, not clearly worth doing
	CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.x ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	CkPrintf("[%d %d %d %d %d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numPoints,j,thisIndex.x+j,thisIndex.w);
#endif
#ifndef CMK_OPTIMIZE
	StartTime=CmiWallTimer();
#endif
	int outOffset=thisIndex.z;
	mcastGrp->contribute(numPoints*sizeof(complex),othernewData+j*numPoints, sumMatrixDoubleType, otherResultCookies[j], mycb, outOffset);
#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif

	//mcastGrp->contribute(numPoints*sizeof(complex),othernewData+j*numPoints, CkReduction::sum_double, otherResultCookies[j], mycb);

      }
  }
  CkAssert(mynewData!=NULL);
  for(int j=0;j<grainSize;j++) //mynewdata
    {
      CkCallback mycb(cp_entry, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
      CkPrintf("[%d %d %d %d %d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,numPoints,j,thisIndex.y+j,thisIndex.w);
#endif
#ifndef CMK_OPTIMIZE
      StartTime=CmiWallTimer();
#endif
      int outOffset=thisIndex.z;
      mcastGrp->contribute(numPoints*sizeof(complex), mynewData+j*numPoints, sumMatrixDoubleType, resultCookies[j], mycb, outOffset);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
      //mcastGrp->contribute(numPoints*sizeof(complex), mynewData+j*numPoints, CkReduction::sum_double, resultCookies[j], mycb);
    }
  delete [] mynewData;
  mynewData=NULL;
  existsNew=false;
  if(otherdata)
    delete [] othernewData;
  othernewData=NULL;
}


#if CONVERSE_VERSION_ELAN

typedef void (* ELAN_REDUCER)(void *in, void *inout, int *count, void *handle);

extern "C" void elan_machine_reduce(int nelem, int size, void * data,
                                    void *dest, ELAN_REDUCER fn, int root);
extern "C" void elan_machine_allreduce(int nelem, int size, void * data,
                                       void *dest, ELAN_REDUCER fn);
#endif

void PairCalcReducer::acceptContribute(int size, double* matrix, CkCallback cb,
                                       bool isAllReduce, bool symmtype, int offx, int offy, int grainSize)
{
    this->isAllReduce = isAllReduce;
    this->size = size;
    this->symmtype = symmtype;
    this->cb = cb;

#if CONVERSE_VERSION_ELAN
    reduction_elementCount ++;
#ifdef _ELAN_PAIRCALC_DEBUG_
    CkPrintf("Accept contribute called on %d count is %d out of %d\n",CkMyPe(),reduction_elementCount, localElements[symmtype].length());
#endif
    int red_size = size *sizeof(double);
    if(tmp_matrix == NULL) {
        tmp_matrix = new double[size];
	bzero(tmp_matrix,red_size);
    }
    int off=offx*grainSize+offy;
    for(int i = 0; i < grainSize*grainSize ; i ++){
        tmp_matrix[i+off] += matrix[i];
    }

    if(reduction_elementCount == localElements[symmtype].length()) {

#ifdef _ELAN_PAIRCALC_DEBUG_
    CkPrintf("[%d] Contributing %d \n",CkMyPe(),reduction_elementCount);

#endif
        reduction_elementCount = 0;
        contribute(sizeof(int),&reduction_elementCount,CkReduction::sum_int);

    }
#else
    CkAbort("Converse Version Is not ELAN, h/w reduction is not supported");
#endif
}


void PairCalcReducer::startMachineReduction() {
#if CONVERSE_VERSION_ELAN
    double * dst_matrix =  NULL;
#ifdef _ELAN_PAIRCALC_DEBUG_
    CkPrintf("Starting machine reduction %d\n",CkMyPe());
#endif
    if(isAllReduce) {
        dst_matrix = new double[size];
        bzero(dst_matrix, size * sizeof(double));
        elan_machine_allreduce(size, sizeof(double), tmp_matrix, dst_matrix, add_double);
    }
    else {
        int pe = CkNumPes()/2; //HACK FOO BAR, GET IT FROM CALLBACK cb

        if(pe == CkMyPe()) {
            dst_matrix = new double[size];
            bzero(dst_matrix, size * sizeof(double));
        }

        elan_machine_reduce(size, sizeof(double), tmp_matrix, dst_matrix, add_double, pe);
    }

    if(isAllReduce) {
         broadcastEntireResult(size, dst_matrix,  symmtype);
        delete [] dst_matrix;
    }
    else {
        if(dst_matrix != NULL) {
            cb.send(size *sizeof(double), dst_matrix);
            delete [] dst_matrix;
        }
    }
    delete [] tmp_matrix;
    tmp_matrix = NULL;
#ifdef _ELAN_PAIRCALC_DEBUG_
    CkPrintf("Leaving machine reduction %d\n",CkMyPe());
#endif
#else
    CkAbort("Converse Version Is not ELAN, h/w reduction is not supported");
#endif


}

void
PairCalcReducer::broadcastEntireResult(entireResultMsg *imsg)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("bcasteres:On Pe %d -- %d objects\n", CkMyPe(), localElements[imsg->symmetric].length());
#endif
  CkAbort("you should not be in this code!\n");
  /*  multiplyResultMsg *msg= new (size, 0, 0) multiplyResultMsg;
  msg->init1(imsg->size,imsg->matrix);
  for (int i = 0; i < localElements[imsg->symmetric].length(); i++)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("call accept on :");
      (localElements[imsg->symmetric])[i].dump();
#endif
      CkSendMsgArrayInline(CkIndex_PairCalculator::__idx_multiplyResultI_multiplyResultMsg, msg, (localElements[imsg->symmetric])[i].id, (localElements[imsg->symmetric])[i].idx, CK_MSG_KEEP);

    }
  delete imsg;
  */
}


void
PairCalcReducer::broadcastEntireResult(int size, double* matrix, bool symmtype)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("bcasteres:On Pe %d -- %d objects\n", CkMyPe(), localElements[symmtype].length());
#endif
  CkAbort("you should not be in this code\n");
  /*
  multiplyResultMsg *msg= new (size,0,0) multiplyResultMsg;
  msg->init1(size,matrix);
  for (int i = 0; i < localElements[symmtype].length(); i++)
    {
#ifdef _PAIRCALC_DEBUG_
      CkPrintf("call accept on :");
      (localElements[symmtype])[i].dump();
#endif
      CkSendMsgArrayInline(CkIndex_PairCalculator::__idx_multiplyResultI_multiplyResultMsg, msg, (localElements[symmtype])[i].id, (localElements[symmtype])[i].idx, CK_MSG_KEEP);

    }
  
  //  delete msg;
  */

}

void
PairCalcReducer::broadcastEntireResult(int size, double* matrix1, double* matrix2, bool symmtype){
    CkPrintf("On Pe %d -- %d objects\n", CkMyPe(), localElements[symmtype].length());
    CkAbort("you should not be i nthis code \n");
    /*
    multiplyResultMsg *msg= new (size,size,0) multiplyResultMsg;
    msg->init(size,size,matrix1,matrix2);
    for (int i = 0; i < localElements[symmtype].length(); i++)
      {
	CkSendMsgArrayInline(CkIndex_PairCalculator::__idx_multiplyResultI_multiplyResultMsg, msg, (localElements[symmtype])[0].id, (localElements[symmtype])[i].idx,CK_MSG_KEEP);
      }
    //    delete msg;
    */
}


void
PairCalcReducer::broadcastEntireResult(entireResultMsg2 *imsg)
{
    CkPrintf("On Pe %d -- %d objects\n", CkMyPe(), localElements[imsg->symmetric].length());
    CkAbort(" you should not be in this code \n");
    /*
    multiplyResultMsg *msg= new (imsg->size,imsg->size,0) multiplyResultMsg;
    msg->size=imsg->size;
    msg->size2=imsg->size;
    memcpy(msg->matrix1,imsg->matrix1,imsg->size* sizeof(double));
    memcpy(msg->matrix2,imsg->matrix2,msg->size2* sizeof(double));
    for (int i = 0; i < localElements[imsg->symmetric].length(); i++)
      {
	CkSendMsgArrayInline(CkIndex_PairCalculator::__idx_multiplyResultI_multiplyResultMsg, msg, (localElements[imsg->symmetric])[0].id, (localElements[imsg->symmetric])[i].idx,CK_MSG_KEEP);
      }
    delete imsg;
    */
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

#include "ckPairCalculator.def.h"


