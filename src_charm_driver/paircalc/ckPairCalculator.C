//**************************************************************************
 /** \file ckPairCalculator.C						   *
 * This is a matrix multiply library with extra frills to communicate the  *
 * results back to gspace or the calling ortho char as directed by the     *
 * callback.                                                               *
 *                                                                         *
 * The extra complications are for parallelization and the multiplication  *
 * of the forces and energies.                                             *
 *                                                                         *
 * In normal use the calculator is created.  Then forces are sent to it    *
 * and multiplied in a big dgemm.  Then this result is reduced to the      *
 * answer matrix and shipped back.  The received left and/or right data is *
 * retained for the backward pass which is triggered by the finishPairCalc *
 * call.  This carries in another set of matrices for multiplication.      *
 * The results are again reduced and cast back.  Thus terminating the life *
 * cycle of the data in the pair calculator.  As the calculator will be    *
 * reused again throughout each iteration the calculators themselves are   *
 * only created once.                                                      *
 *                                                                         *
 * The elan code is a specialized machine reduction/broadcast which runs   *
 * much faster on lemieux, and theoretically on any other elan machine     *
 * It operates by one dummy reduction which is used to indicate that all   *
 * calculators on a PE machine have reported in.  Then the machine         *
 * reduction is triggered.                                                 *
 *                                                                         *
 * The paircalculator is a 4 dimensional array.  Those dimensions are:     *
 *            w: gspace state plane (the second index of the 2D gspace)    *
 *            x: coordinate offset within plane (a factor of grainsize)    *
 *            y: coordinate offset within plane (a factor of grainsize)    *
 *            z: blocksize, unused always 0                                *
 *       So, for an example grainsize of 64 for a 128x128 problem:         *
 *        S/grainsize gives us a 2x2 decomposition.                        *
 *        1st quadrant ranges from [0,0]   to [63,63]    index [w,0,0,0]   *
 *        2nd quadrant ranges from [0,64]  to [63,127]   index [w,0,64,0]  *
 *        3rd quadrant ranges from [64,0]  to [127,63]   index [w,64,0,0]  *
 *        4th quadrant ranges from [64,64] to [127,127]  index [w,64,64,0] *
 *                                                                         *
 *       0   64   127                                                      *
 *     0 _________                                                         *
 *       |   |   |                                                         *
 *       | 1 | 2 |                                                         *
 *    64 ---------                                                         *
 *       |   |   |                                                         *
 *       | 3 | 4 |                                                         *
 *   127 ---------                                                         *
 *                                                                         *
 *                                                                         *
 *                                                                         *
 * Further complication arises from the fact that each plane is a          *
 * cross-section of a sphere.  So the actual data is sparse and is         *
 * represented by a contiguous set of the nonzero elements.  This is       * 
 * commonly referred to as N or size within the calculator.                *
 *                                                                         */
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

PairCalculator::PairCalculator(bool sym, int grainSize, int s, int blkSize, CkCallback cb, CkArrayID cb_aid, int cb_ep, bool conserveMemory, bool lbpaircalc,  redtypes _cpreduce)
{
#ifdef _PAIRCALC_DEBUG_PLACE_
  CkPrintf("[PAIRCALC] [%d %d %d %d %d] inited on pe %d \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,sym, CkMyPe());
#endif
  this->conserveMemory=conserveMemory;
  this->symmetric = sym;
  this->grainSize = grainSize;
  this->S = s;
  this->blkSize = blkSize;
  this->cb = cb;
  this->N = -1;
  this->cb_aid = cb_aid;
  this->cb_ep = cb_ep;
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
  usesAtSync=true;
  if(lbpaircalc)
      setMigratable(true);
  else
      setMigratable(false);
  resultCookies=NULL;
  otherResultCookies=NULL;
  resultCookies=new CkSectionInfo[grainSize];
  if(symmetric && (thisIndex.x != thisIndex.y))
    otherResultCookies=new CkSectionInfo[grainSize];

}

void
PairCalculator::pup(PUP::er &p)
{
  ArrayElement4D::pup(p);
  p|numRecd;
  p|numExpected;
  p|grainSize;
  p|S;
  p|blkSize;
  p|N;
  p|symmetric;
  p|conserveMemory;
  p|lbpaircalc;
  p|cpreduce;
  p|cb;
  p|cb_aid;
  p|cb_ep;
  p|existsLeft;
  p|existsRight;
  p|existsOut;
  p|existsNew;
  p|mCastGrpId;
  p|cookie;
  p|resumed;
  p|rck;
  if (p.isUnpacking()) 
  {
      mynewData=NULL;
      othernewData=NULL;
      resultCookies=new CkSectionInfo[grainSize];
      if(symmetric && (thisIndex.x != thisIndex.y))
	  otherResultCookies= new CkSectionInfo[grainSize];
      else
	  otherResultCookies=NULL;
      if(existsOut)
	  outData= new double[grainSize*grainSize];
      else
	  outData=NULL;
      if(existsLeft)
	  inDataLeft = new double[2*numExpected*N];
      else
	  inDataLeft=NULL;
      if(existsRight)
	  inDataRight = new double[2*numExpected*N];
      else
	  inDataRight=NULL;

  }
  int i;
  for (i=0; i<grainSize; i++) p|resultCookies[i];
  if(symmetric && (thisIndex.x != thisIndex.y))
    for (i=0; i<grainSize; i++) p|otherResultCookies[i];
  for (int i=0; i<grainSize; i++)
    CmiAssert(resultCookies[i].get_redNo() > 0);
  if(existsOut)
    p(outData, grainSize*grainSize);
  if(existsLeft)
    p(inDataLeft, numExpected * N * 2);
  if(existsRight)
    p(inDataRight, numExpected* N * 2);

#ifdef _PAIRCALC_DEBUG_
  if (p.isUnpacking())
    {
      CkPrintf("[%d,%d,%d,%d,%d] pup unpacking on %d resumed=%d memory %d\n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z,symmetric,CkMyPe(),resumed, CmiMemoryUsage());
      CkPrintf("[%d,%d,%d,%d,%d] pupped : %d %d %d %d %d %d %d %d %d  %d %d cb cb_aid %d %d %d cb_lb inDataLeft inDataRight outData  %d \n",thisIndex.w,thisIndex.x, thisIndex.y, thisIndex.z, symmetric, numRecd, numExpected, grainSize, S, blkSize, N, symmetric, conserveMemory, lbpaircalc, cpreduce, cb_ep, existsLeft, existsRight,  resumed);

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
  if(resultCookies!=NULL)
    delete [] resultCookies;
  if(symmetric && (thisIndex.x != thisIndex.y) && otherResultCookies !=NULL)
    delete [] otherResultCookies;
}






// initialize the section cookie and the reduction client
void PairCalculator::initGRed(initGRedMsg *msg)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d,%d,%d,%d,%d] initGRed\n",thisIndex.w,thisIndex.x,thisIndex.y, thisIndex.z, symmetric);
#endif
  CkGetSectionInfo(cookie,msg);
  cb=msg->cb;
  mCastGrpId=msg->mCastGrpId; //redundant
  cpreduce=section;
  if(msg->lbsync)
  {
      int foo=1;
      contribute(sizeof(int), &foo , CkReduction::sum_int, msg->synccb);
  }
  //    delete msg; do not delete nokeep method
}

//!initialize the section cookie for each slice of the result
void PairCalculator::initResultSection(initResultMsg *msg)
{

  CkAssert(msg->offset<grainSize);
  if(symmetric && msg->dest == thisIndex.x && thisIndex.x != thisIndex.y)
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
 * arrive it multiplies and contributes.
 */
void
PairCalculator::calculatePairs_gemm(calculatePairsMsg *msg)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf(" symm=%d    pairCalc[%d %d %d %d] got from [%d %d] with size {%d}, from=%d, count=%d, resumed=%d\n", symmetric, thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,  thisIndex.w, msg->sender, msg->size, msg->fromRow, numRecd,resumed);
#endif
  if(!resumed)
    {
      ResumeFromSync();
    }
  numRecd++;   // increment the number of received counts
  int offset = -1;
  if (msg->fromRow) {   // This could be the symmetric diagonal case
    offset = msg->sender - thisIndex.x;
    if (!existsLeft)
      { // now that we know N we can allocate contiguous space
	CkAssert(inDataLeft==NULL);
	existsLeft=true;
	N = msg->size; // N is init here with the size of the data chunk.
	inDataLeft = new double[numExpected*N*2];
	bzero(inDataLeft,numExpected*N*2*sizeof(double));
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d] Allocated Left %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric, numExpected,N);

#endif
      }
    CkAssert(N==msg->size);
    CkAssert(offset<numExpected);
    memcpy(&(inDataLeft[offset*N*2]), msg->points, N * 2 *sizeof(double));
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Copying into offset*N %d * %d N *2 %d points start %g end %g\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z, symmetric, offset, N, N*2,msg->points[0].re, msg->points[N-1].im);
#endif
#ifdef _PAIRCALC_DEBUG_STUPID_PARANOID_
    CkPrintf("copying in \n");
    double re;
    double im;
    for(int i=0;i<N;i++)
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
      { // now that we know N we can allocate contiguous space
	CkAssert(inDataRight==NULL);
	existsRight=true;
	N = msg->size; // N is init here with the size of the data chunk.
	inDataRight = new double[numExpected*N*2];
	bzero(inDataRight,numExpected*N*2*sizeof(double));
#ifdef _PAIRCALC_DEBUG_
	CkPrintf("[%d,%d,%d,%d,%d] Allocated right %d * %d *2 \n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,numExpected,N);
#endif
      }
    CkAssert(N==msg->size);
    CkAssert(offset<numExpected);
    memcpy(&(inDataRight[offset*N*2]), msg->points, N * 2 *sizeof(double));
#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d,%d,%d,%d,%d] Copying into offset*N %d * %d N *2 %d points start %g end %g\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric,offset,N, N*2,msg->points[0].re, msg->points[N-1].im);
#endif
#ifdef _PAIRCALC_DEBUG_STUPID_PARANOID_
    CkPrintf("copying in \n");
    double re;
    double im;
    for(int i=0;i<N;i++)
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
   * (numExpected X N) X (N X numExpected) = (numExpected X numExpected)
   */
  /* To make this work, we transpose the first matrix (A).
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

  if (numRecd == numExpected * 2 || (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected)) {
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
    int doubleN=2*N;
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
#ifdef _PAIRCALC_DEBUG_PARANOID_
	char filename[80];
	sprintf(filename, "fwlmdata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
	FILE *loutfile = fopen(filename, "w");
	fprintf(loutfile,"[%d,%d,%d,%d,%d] inDataLeft\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
	for(int i=0;i<N*2;i++)
	  for(int j=0;j<numExpected;j++)
	    fprintf(loutfile,"%d %d %g \n",i,j,inDataLeft[i*numExpected+j]);
	fclose(loutfile);
	sprintf(filename, "fwrmdata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
	FILE *routfile = fopen(filename, "w");
	fprintf(routfile,"[%d,%d,%d,%d,%d] inDataRight\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
	for(int i=0;i<N*2;i++)
	  for(int j=0;j<numExpected;j++)
	    fprintf(routfile,"%d %d %g\n",i,j,inDataRight[i*numExpected+j]);
	fclose(routfile);
#endif

	DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, inDataRight, &lda, inDataLeft, &ldb, &beta, outData, &ldc);
	// switching solves a transposition problem
	//	  DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, inDataLeft, &lda, inDataRight, &ldb, &beta, outData, &ldc);
      }
    else if (symmetric && thisIndex.x==thisIndex.y && numRecd==numExpected)
      {
#ifdef _PAIRCALC_DEBUG_PARANOID_
	char filename[80];
	sprintf(filename, "fwlmdata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
	FILE *loutfile = fopen(filename, "w");
	fprintf(loutfile,"[%d,%d,%d,%d,%d] inDataLeft\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
	for(int i=0;i<N*2;i++)
	  for(int j=0;j<numExpected;j++)
	    fprintf(loutfile,"%d %d %g \n",i,j,inDataLeft[i*numExpected+j]);
	fclose(loutfile);
#endif
	DGEMM(&transformT, &transform, &m_in, &n_in, &k_in, &alpha, inDataLeft, &lda, inDataLeft, &ldb, &beta, outData, &ldc);
      } 

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(210, StartTime, CmiWallTimer());
#endif
    numRecd = 0; 
    if (msg->flag_dp) {
	for (int i = 0; i < grainSize*grainSize; i++)
	  outData[i] *= 2.0;
    }
#ifdef _PAIRCALC_DEBUG_PARANOID_
    char filename[80];
    sprintf(filename, "fwgmodata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
    FILE *outfile = fopen(filename, "w");


    fprintf(outfile,"[%d,%d,%d,%d,%d] outData=C\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
    for(int i=0;i<grainSize;i++)
      {
	for(int j=0;j<grainSize;j++)
	  fprintf(outfile," %g",outData[i*grainSize+j]);
	  fprintf(outfile,"\n");
      }
    fclose(outfile);
#endif

#ifdef _PAIRCALC_VALID_OUT_
    CkPrintf("[PAIRCALC] [%d %d %d %d %d] forward gemm out %.10g %.10g %.10g %.10g \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric, outData[0],outData[1],outData[grainSize*grainSize-2],outData[grainSize*grainSize-1]);
#endif
#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif

    CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    mcastGrp->contribute(grainSize*grainSize*sizeof(double), outData, sumMatrixDoubleType, cookie, cb);
    //mcastGrp->contribute(grainSize*grainSize*sizeof(double), outData, CkReduction::sum_double, cookie, cb);


#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(220, StartTime, CmiWallTimer());
#endif
  }

}

//PairCalculator::multiplyResult(int size, double *matrix1, double *matrix2)
void
PairCalculator::multiplyResultI(multiplyResultMsg *msg)
{
  
    multiplyResult(msg);
  }


// entry method to allow us to delay this outbound communication
// to minimize brain dead BG/L interference we have a signal to prioritize this
/* Send the results via multiple reductions */
void
PairCalculator::sendBWResult(sendBWsignalMsg *msg)
{

#ifdef _PAIRCALC_DEBUG_
    CkPrintf("[%d %d %d %d %d]: sendBWResult\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric);
#endif
    delete msg; // contains nothing
    // Now we have results in mynewData and if(symmetric) othernewData
    CkMulticastMgr *mcastGrp=CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();

    if(symmetric){  // we have this othernewdata issue for the symmetric case
	if (thisIndex.y != thisIndex.x) 

	{ //othernewdata exists for we are symmetric and not on the diagonal
	    CkAssert(othernewData!=NULL);
	    for(int j=0;j<grainSize;j++)
	    {

		CkCallback mycb(cb_ep, CkArrayIndex2D(j+thisIndex.x ,thisIndex.w), cb_aid);

#ifdef _PAIRCALC_DEBUG_CONTRIB_
		CkPrintf("[%d %d %d %d %d] contributing other %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, N,j,thisIndex.x+j,thisIndex.w);
#endif
		mcastGrp->contribute(N*sizeof(complex),othernewData+j*N, sumMatrixDoubleType, otherResultCookies[j], mycb);
		//mcastGrp->contribute(N*sizeof(complex),othernewData+j*N, CkReduction::sum_double, otherResultCookies[j], mycb);

	    }
	}
    }
    CkAssert(mynewData!=NULL);
    for(int j=0;j<grainSize;j++) //mynewdata
    {
	CkCallback mycb(cb_ep, CkArrayIndex2D(j+thisIndex.y ,thisIndex.w), cb_aid);
#ifdef _PAIRCALC_DEBUG_CONTRIB_
	CkPrintf("[%d %d %d %d %d] contributing %d offset %d to [%d %d]\n",thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric,N,j,thisIndex.y+j,thisIndex.w);
#endif
	mcastGrp->contribute(N*sizeof(complex), mynewData+j*N, sumMatrixDoubleType, resultCookies[j], mycb);
	//mcastGrp->contribute(N*sizeof(complex), mynewData+j*N, CkReduction::sum_double, resultCookies[j], mycb);
    }
    delete [] mynewData;
    mynewData=NULL;
    existsNew=false;
    if(symmetric && thisIndex.x != thisIndex.y)
	delete [] othernewData;
    othernewData=NULL;

}

void
PairCalculator::multiplyResult(multiplyResultMsg *msg)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("[%d %d %d %d %d]: Accept Result with size %d\n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z, symmetric, msg->size);
#endif
  int size=msg->size;
  int size2=msg->size2;
  double *matrix1=msg->matrix1;
  double *matrix2=msg->matrix2;
  bool unitcoef = false;
  if(cpreduce==section) //update cookie from multicast
      CkGetSectionInfo(cookie,msg);
  if(matrix2==NULL||size2<1) 
    {
      unitcoef = true;
    }
  else
      CkPrintf("Hey! matrix2 isn't null!\n");
  CkAssert(mynewData==NULL);
  mynewData = new complex[N*grainSize];
  existsNew=true;
  if(symmetric && (thisIndex.x != thisIndex.y)){
    othernewData = new complex[N*grainSize];
    bzero(othernewData,N*grainSize* sizeof(complex));
  }
  else
      othernewData=NULL;

  int offset = 0, index = thisIndex.y*S + thisIndex.x;

  if(!symmetric)
    index = thisIndex.x*S + thisIndex.y;
  int matrixSize=grainSize*grainSize;
  //ASSUMING TMATRIX IS REAL (LOSS OF GENERALITY)
  register double m=0;
  bool makeLocalCopy=false;
  double *amatrix=NULL;

  if(S!=grainSize && size!=matrixSize)
      CkAbort("we don't support PC striding anymore");
  
  if(S==grainSize)// all at once no malloc
    {
      amatrix=matrix1+index;
    }
  else
    { // you were sent the correct section only
      amatrix=matrix1;
    }


  double *localMatrix;
  double *outMatrix;

  int m_in=grainSize;
  int n_in=N*2;
  int k_in=grainSize;
  double *mynewDatad= reinterpret_cast <double *> (mynewData);
  double alphad(1.0);
  double betad(0.0);
  char transform='N';
  char transformT='T';

#ifdef _PAIRCALC_DEBUG_PARANOID_
  char filename[80];
  sprintf(filename, "bwlmodata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
  FILE *outfile = fopen(filename, "w");


  fprintf(outfile,"[%d,%d,%d,%d,%d] LM\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
  for(int i=0;i<N*grainSize;i++)
    {
      fprintf(outfile," %g %g",inDataLeft[i]);
      fprintf(outfile,"\n");
    }
  fclose(outfile);
  sprintf(filename, "bwmdata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
  FILE *moutfile = fopen(filename, "w");
  fprintf(moutfile,"[%d,%d,%d,%d,%d] amatrix\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
  for (int i = 0; i < grainSize; i++) {
    for (int j = 0; j < grainSize; j++){ 
      fprintf(moutfile,"%.10g\n",amatrix[i*grainSize+j]);
    }
  }
  fclose(moutfile);
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
      //this isn't right
      DGEMM(&transform, &transformT, &n_in, &m_in, &k_in, &alphad, inDataRight, &n_in,  amatrix, &k_in, &betad, othernewDatad, &n_in);
#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(250, StartTime, CmiWallTimer());
#endif
    }


#ifdef _PAIRCALC_DEBUG_PARANOID_
  sprintf(filename, "bwgmodata.%d_%d_%d_%d_%d", thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z, symmetric);
  FILE *ooutfile = fopen(filename, "w");


  fprintf(ooutfile,"[%d,%d,%d,%d,%d] outData=C\n",thisIndex.w,thisIndex.x,thisIndex.y,thisIndex.z,symmetric);
  for(int i=0;i<N*grainSize;i++)
    {
      fprintf(ooutfile," %g %g",mynewData[i].re,mynewData[i].im);
      fprintf(ooutfile,"\n");
    }
  fclose(ooutfile);
#endif


  if(!unitcoef){ // CG non minimization case
      amatrix=matrix2;
      double *rightd=reinterpret_cast <double *> (inDataRight);
    // C = alpha*A*B + beta*C
#ifndef CMK_OPTIMIZE
    StartTime=CmiWallTimer();
#endif
    betad=1.0;
    DGEMM(&transform, &transformT, &n_in, &m_in, &k_in, &alphad, inDataRight, &n_in,  amatrix, &k_in, &betad, mynewDatad, &n_in);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(240, StartTime, CmiWallTimer());
#endif
  } // end  CG case



#ifdef _PAIRCALC_VALID_OUT_
  CkPrintf("[PAIRCALC] [%d %d %d %d %d] backward gemm out %.10g %.10g %.10g %.10g \n", thisIndex.w, thisIndex.x, thisIndex.y, thisIndex.z,symmetric, mynewDatad[0],mynewDatad[1],mynewData[N*grainSize-1].re,mynewData[N*grainSize-1].im);
#endif

  sendBWsignalMsg *sigmsg=new (8*sizeof(int)) sendBWsignalMsg;
  CkSetQueueing(sigmsg, CK_QUEUEING_IFIFO);
  *(int*)CkPriorityPtr(sigmsg) = 1; // just make it slower than non prioritized
  thisProxy(thisIndex.w,thisIndex.x, thisIndex.y,thisIndex.z).sendBWResult(sigmsg);
  //sendBWResult(sigmsg);

  if(conserveMemory)
    {
      // clear the right and left they'll get reallocated on the next pass
      delete [] inDataLeft;
      inDataLeft=NULL;
      if(!symmetric || (symmetric&&thisIndex.x!=thisIndex.y)) {
	delete [] inDataRight;
	inDataRight = NULL;
      }
      existsLeft=false;
      existsRight=false;
      if(outData!=NULL)
	{
	  delete [] outData;
	  outData = NULL;
	  existsOut=false;
	}
    }
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
  multiplyResultMsg *msg= new (size,0,0) multiplyResultMsg;
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
}


void
PairCalcReducer::broadcastEntireResult(int size, double* matrix, bool symmtype)
{
#ifdef _PAIRCALC_DEBUG_
  CkPrintf("bcasteres:On Pe %d -- %d objects\n", CkMyPe(), localElements[symmtype].length());
#endif
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


}

void
PairCalcReducer::broadcastEntireResult(int size, double* matrix1, double* matrix2, bool symmtype){
    CkPrintf("On Pe %d -- %d objects\n", CkMyPe(), localElements[symmtype].length());
    multiplyResultMsg *msg= new (size,size,0) multiplyResultMsg;
    msg->init(size,size,matrix1,matrix2);
    for (int i = 0; i < localElements[symmtype].length(); i++)
      {
	CkSendMsgArrayInline(CkIndex_PairCalculator::__idx_multiplyResultI_multiplyResultMsg, msg, (localElements[symmtype])[0].id, (localElements[symmtype])[i].idx,CK_MSG_KEEP);
      }
    //    delete msg;

}

void
PairCalcReducer::broadcastEntireResult(entireResultMsg2 *imsg)
{
    CkPrintf("On Pe %d -- %d objects\n", CkMyPe(), localElements[imsg->symmetric].length());
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

}


#include "ckPairCalculator.def.h"


