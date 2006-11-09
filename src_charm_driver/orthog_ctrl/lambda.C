//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file lambda.C
 * The Lambda object is basically used to accept the rows of the "L"
 * matrix, so they can be summed across the PairCalculators.  Then
 * the summed (reduced) form is sent back to the calculators for
 * application to their input (gspace) for the asymmetric backward path.
 *
 * For the dynamics (non minization) case we have an additional
 * multiplication of gamma= lambda*lambdaT.  And we pass both gamma and
 * lambdaT to the PC instead.
 *
 * To accomplish that multiply we have to stitch together lambdas
 * corresponding to the parts of orthoT in the Ortho object.  Dynamics
 * will trigger a send of lambda into ortho to overlap with the S->T
 * process.  So that the gamma multiply can proceed once we have T.
 *
 */
//============================================================================

#include "charm++.h"
#include "../../include/debug_flags.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"
#include "lambda.h"
#include <unistd.h>
#include "../../src_mathlib/mathlib.h"


//============================================================================

extern Config config;
extern int nstates;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_CP_State_GSpacePlane gSpacePlaneProxy;
extern  PairCalcID PairCalcID1;
extern  PairCalcID pairCalcID2;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_CP_Rho_GHartExt rhoGHartExtProxy;



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Lambda::acceptAllLambda(CkReductionMsg *msg) {
    delete msg;
    CkAbort("do not call acceptAllLambda");
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Lambda::acceptSectionLambda(CkReductionMsg *msg) {
//============================================================================

  double *lambda = (double *)msg->getData();
  int lambdaCount = msg->getSize()/sizeof(double);

#ifdef _CP_DEBUG_LMAT_
  char lmstring[80];
  snprintf(lmstring,80,"lmatrix_t:%d_%d_%d.out",numGlobalIter,thisIndex.x,thisIndex.y);
  FILE *fp = fopen(lmstring,"w");
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fprintf(fp, "[%d %d] %.12g \n", i + n*thisIndex.x+1, j+n*thisIndex.y+1, lambda[i*n+j]);
    }
  }
  fclose(fp);
#endif
  // finish pair calc
  finishPairCalcSection(lambdaCount, lambda, &lPairCalcID2, thisIndex.x, thisIndex.y, 0, lPairCalcID2.priority+1);
#ifdef _CP_DEBUG_LAMBDA_VERBOSE_
  if(thisIndex.x==0 && thisIndex.y==0)
    CkPrintf("[%d,%d] finishing asymm\n",thisIndex.x, thisIndex.y);
#endif
  delete msg;

//----------------------------------------------------------------------------
}// end routine
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void Lambda::makeSections(int indexSize, int *indexZ){

  int s1=thisIndex.x*lambdaGrainSize;
  int s2=thisIndex.y*lambdaGrainSize;
  if(lambdaGrainSize!=sGrainSize)
    {
      // do something clever
      s1=s1/sGrainSize*sGrainSize;
      s2=s2/sGrainSize*sGrainSize;
    }

  // thisIndex.x and thisIndex.y range from 0 to nstates/lambdaGrainSize
    initOneRedSect(indexSize, indexZ, config.numChunksAsym, &lPairCalcID2, CkCallback(CkIndex_Lambda::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)) , s1, s2, thisIndex.x, thisIndex.y, lambdaGrainSize, false, false, config.useOrthoDirect);
     initOneRedSect(indexSize, indexZ, config.numChunksAsym, &lPairCalcID2, CkCallback(CkIndex_Lambda::acceptSectionLambda(NULL), thisProxy(thisIndex.x, thisIndex.y)) , s1, s2, thisIndex.x, thisIndex.y, lambdaGrainSize, false, true, config.useOrthoDirect);

//----------------------------------------------------------------------------
  }// end routine
//==============================================================================



