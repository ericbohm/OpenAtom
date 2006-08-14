// Things to do
// dofft should do the fft

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_State_RealSpacePlane.C
 * This is a description of the "life" of a CP_State_RealSpacePlane object.
 *
 * At the start of the program, the constructor CP_State_RealSpacePlane() is called.
 * The CP_State_GSpacePlane objects send data to CP_State_RealSpacePlane 
 * after doing FFT in the y direction. The data is sent through the doFFT() method. 
 * In this method FFTs in the z and x directions are performed. After this the 
 * squared magnitudes of the psi_r values are sent to RealSpaceDensity. The 
 * CP_State_RealSpacePlane object is idle until further messages are sent to it.
 *
 * The idle period of the CP_State_RealSpacePlane is ended by a message from 
 * RealSpaceDensity objects - doProduct(). In this method the Slab of data
 * in CP_State_RealSpacePlane, psi_r, is multiplied with the data sent from 
 * RealSpaceDensity. Then, inverse ffts in z and x directions are performed.
 * Then the slab of data is send to CP_State_GSpacePlanes so that the inverse fft
 * in the x direction can be performed.
 */
//============================================================================

#include "charm++.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "../../include/debug_flags.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"


//============================================================================

extern CProxy_CP_State_GSpacePlane    gSpacePlaneProxy;
extern CProxy_CP_Rho_RealSpacePlane   rhoRealProxy;
extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CPcharmParaInfoGrp      scProxy;
extern CProxy_main                    mainProxy;
extern CProxy_CP_State_ParticlePlane  particlePlaneProxy;
extern CProxy_FFTcache                fftCacheProxy;

extern CkGroupID            mCastGrpId;
extern ComlibInstanceHandle mssInstance;

extern int    sizeX;
extern Config config;

//============================================================================

//#define _CP_DEBUG_STATER_VERBOSE_

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
RTH_Routine_locals(CP_State_RealSpacePlane,run)
RTH_Routine_code(CP_State_RealSpacePlane,run) {
//============================================================================

  while(1) { 
    // constructor invokes run and then you suspend (no work yet)
    RTH_Suspend(); 
    c->doFFT();    // state(g,z) from gstate arrives in dofft(msg) which resumes
#ifndef _CP_DEBUG_RHO_OFF_
    RTH_Suspend(); // after doreduction sends data to rhoreal, suspend
#endif
    c->doProductThenFFT(); // vks(r) arrives in doproduct(msg) which resumes
    c->sendFPsiToGSP();
  } //end while not done

//--------------------------------------------------------------------------
   } RTH_Routine_end(CP_State_RealSpacePlane,run)
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::run () {
  run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_State_RealSpacePlane,run),this);
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_State_RealSpacePlane::CP_State_RealSpacePlane(size2d size, int gSpaceUnits, 
                  int realSpaceUnits, int _ngrida, int _ngridb, int _ngridc) {
//============================================================================
//  ckout << "State R Space Constructor : "
//	<< thisIndex.x << " " << thisIndex.y << " " <<CkMyPe() << endl;
//============================================================================

    count = 0;

    ngrida = _ngrida;
    ngridb = _ngridb;
    ngridc = _ngridc;
    csize = (ngrida/2 + 1)*ngridb; 
    rsize = (ngrida   + 2)*ngridb; ;

    initRealStateSlab(&rs, size, gSpaceUnits, realSpaceUnits, thisIndex.x, thisIndex.y);

    gproxy = gSpacePlaneProxy;
    if (config.useMssInsGP){
      ComlibAssociateProxy(&mssInstance,gproxy);
    }//endif

    setMigratable(false);

    run();

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);

  p|ngrida;
  p|ngridb;
  p|ngridc;
  p|count;
  p|csize;
  p|rsize;
  p|cookie;
  p|gproxy;

  rs.pup(p); // pup your plane, now, honey.

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::setNumPlanesToExpect(int num){
    rs.numPlanesToExpect = num;
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::doFFT(RSFFTMsg *msg) {
//============================================================================

    int size               = msg->size; 
    int Index              = msg->senderIndex;
    complex *partiallyFFTd = msg->data;
    int nchareG            = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
    int **tranUnpack       = scProxy.ckLocalBranch()->cpcharmParaInfo->index_tran_upack;
    int *nline_per_chareG  = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_per_chareG;

    int planeSize          = rs.size;

    count++;
    if (count > nchareG) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Mismatch in allowed gspace chare arrays : %d %d %d %d\n",
                count,nchareG,thisIndex.x,thisIndex.y);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//============================================================================
// Unpack the message
// You have received packed data (x,y) from processor sendIndex
// Every real space chare receives the same x,y indicies.
// For double pack, x=0,1,2,3,4 ...  y= {-K ... K}
// The x increase with processor number. The y are split.
// The rundescriptor contains all we need to unpack the data.
// For doublepack : nffty*run[i][j].x + run[i][j].y
// we store this stuff in the convenient package
// Pictorially a half cylinder is sent which is unpacked into
// a half cube for easy FFTing. Y is the inner index.


    // non-zero elements are set. Zero elements are zeroed here
    // planeSize also contains extra elements on the boundary for fftw
    if(config.conserveMemory && count==1){rs.allocate();}
    complex *planeArr = rs.planeArr;
    if(count==1){bzero(planeArr,planeSize*sizeof(complex));} 

    if(size!=nline_per_chareG[Index]){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Dude, %d != %d for chare %d %d\n",size,nline_per_chareG[Index],
                   thisIndex.y,Index);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    for(int i=0;i< size;i++){planeArr[tranUnpack[Index][i]] = partiallyFFTd[i];}

    delete msg;

//============================================================================
// If every chareG has reported then you can resume/go on and do the FFT

    if (count == nchareG) {
      count=0;
      RTH_Runtime_resume(run_thread);
    }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// After receiving from G-space, FFT to real space
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::doFFT(){
//============================================================================

#ifdef _CP_DEBUG_STATER_VERBOSE_
   ckout << "Real Space " << thisIndex.x << " " << thisIndex.y << " doing FFT" << endl;
#endif

//============================================================================
// Perform the FFT and get psi^2 which we can store in cache tmpData because
// we will blast it right off before losing control

    FFTcache *fftcache  = fftCacheProxy.ckLocalBranch();
    int nplane_x        = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
    complex *planeArr   = rs.planeArr;
    double  *planeArrR  = rs.planeArrR;

#ifndef CMK_OPTIMIZE    
    double StartTime=CmiWallTimer();
#endif

    fftcache->doStpFFTGtoR_Rchare(planeArr,planeArrR,nplane_x,ngrida,ngridb);

    double *data = fftcache->tmpDataR;
    for(int i=0,i2=0;i<ngridb;i++,i2+=2){
      for(int j=i*ngrida;j<(i+1)*ngrida;j++){
        data[j] = planeArrR[(j+i2)]*planeArrR[(j+i2)];
      }//endfor
    }//endfor
    CmiNetworkProgress();

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(doRealFwFFT_, StartTime, CmiWallTimer());
#endif    

//============================================================================
// If non-local itself is on and the ees method is to be used, launch

#ifndef _CP_DEBUG_SFNL_OFF_ // non-local is allowed 
  int nchareG      = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
  int ees_nonlocal = scProxy.ckLocalBranch()->cpcharmParaInfo->ees_nloc_on;
  if(ees_nonlocal==1){
    CkAssert(nchareG<=ngridc);
    if(thisIndex.y<=nchareG){
      gSpacePlaneProxy(thisIndex.x,thisIndex.y).startNLEes();
    }//endif
  }//endif
#endif

//============================================================================
// Send off the reduction unless you are skipping rho, whence you call doproduct

#ifdef _CP_DEBUG_RHO_OFF_
  if(thisIndex.x==0 && thisIndex.y==0){
    CkPrintf("EHART       = OFF FOR DEBUGGING\n");
    CkPrintf("EExt        = OFF FOR DEBUGGING\n");
    CkPrintf("EWALD_recip = OFF FOR DEBUGGING\n");
    CkPrintf("EEXC        = OFF FOR DEBUGGING\n");
    CkPrintf("EGGA        = OFF FOR DEBUGGING\n");
    CkPrintf("EEXC+EGGA   = OFF FOR DEBUGGING\n");
  }//endif
  bzero(data,ngrida*ngridb*sizeof(double));
  doProductThenFFT();
#else
  doReduction();
#endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================




//============================================================================
// No one else can use tmpdataR until I perform doReduction because 
// I will not be rolled out until I finish my scalar work.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::doReduction(){
//============================================================================

#ifdef _CP_DEBUG_STATER_VERBOSE_
    CkPrintf("In StateRSpacePlane[%d %d] doReduction %d\n", thisIndex.x, thisIndex.y,
                CmiMemoryUsage());
#endif

//============================================================================
// Perform the Reduction to get the density : vks holds psi^2 for us
// Return values for vks cannot hit this chare until the reduciton is complete.

   CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
   CkCallback cb(CkIndex_CP_Rho_RealSpacePlane::acceptDensity(0),
                 CkArrayIndex2D(thisIndex.y,0),rhoRealProxy.ckGetArrayID());

#ifndef CMK_OPTIMIZE    
   double StartTime=CmiWallTimer();
#endif

    FFTcache *fftcache  = fftCacheProxy.ckLocalBranch();
    double *data        = fftcache->tmpDataR;
    mcastGrp->contribute(ngrida*ngridb*sizeof(double),data,sumFastDoubleType,
                         cookie,cb);
    CmiNetworkProgress();

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(DoFFTContribute_, StartTime, CmiWallTimer());
#endif    
//============================================================================
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *   In this method, we receive vks from the density. We copy it out 
 *   because you don't want to work with charms squirely memory 
 *   and call the working doproduct. The copy is into cache scratch
 *   which is OK because we should not relinquish the proc to another
 *   chare until we have finsihed the working doproduct
 */
//============================================================================
void CP_State_RealSpacePlane::doProduct(ProductMsg *msg) {
//============================================================================

#ifdef _CP_DEBUG_STATER_VERBOSE_
  CkPrintf("In StateRSpacePlane[%d %d] doProd \n", thisIndex.x, thisIndex.y);
#endif

//============================================================================
// Unpack and check size, copy to cache temp, resume which calls doProduct
// without relinquishing control to other chares

  double *vks_in = msg->data;
  int mysize     = msg->datalen;
  CkAssert(mysize == ngrida*ngridb);
	
  FFTcache *fftcache  = fftCacheProxy.ckLocalBranch();
  double *Vks         = fftcache->tmpDataR;

  CmiMemcpy(Vks,vks_in,sizeof(double)*mysize);
  CmiNetworkProgress();
  
  RTH_Runtime_resume(run_thread); // this is scalar, we continue right on
                                  // as threaded loops calls working doproduct

//---------------------------------------------------------------------------
  }//endroutine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
* Really doing the product : called directly from the doProduct(msg) : I don't 
*                             give way to another chare so I can use the scratch 
*                             space in the Cache group
*/
//============================================================================
void CP_State_RealSpacePlane::doProductThenFFT() {
//============================================================================
// A little output under some circumstances

#ifdef _CP_DEBUG_STATER_VERBOSE_
  CkPrintf("In RealSpacePlane[%d %d] doProduct %d\n",
             thisIndex.x, thisIndex.y,CmiMemoryUsage());
#endif

#ifndef _CP_DEBUG_RHO_OFF_  
#ifdef _CP_DEBUG_VKS_RSPACE_
  if(thisIndex.x==0 && thisIndex.y == 0){
     FILE *fp = fopen("vks_real_y0_state0.out","w");
      for(int i=0;i<size;i++){
        fprintf(fp,"%g\n",vks[i]);
      }//endfor
     fclose(fp);
  }//endif
#endif    
#endif

//===================================================================
// Do the product and then the FFT

 //------------------------------------------------------------------
 // Start the timer
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

 //------------------------------------------------------------------
 // The product psi*vks
  double  *planeArrR  = rs.planeArrR;
  FFTcache *fftcache  = fftCacheProxy.ckLocalBranch();
  double *Vks         = fftcache->tmpDataR;

  int stride = sizeX/2+1;
  for(int i=0,i2=0;i<ngridb;i++,i2+=2){
    for(int j=i*ngrida;j<(i+1)*ngrida;j++){
      planeArrR[(j+i2)] *= Vks[j];
    }//endfor
  }//endfor
  CmiNetworkProgress();

 //------------------------------------------------------------------
 // The FFT
  int nplane_x      = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
  complex *planeArr = rs.planeArr;
  fftcache->doStpFFTRtoG_Rchare(planeArr,planeArrR,nplane_x,ngrida,ngridb);

 //------------------------------------------------------------------
 // End timer 
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(doRealBwFFT_, StartTime, CmiWallTimer());
#endif

//===================================================================
// Return to threaded loop which invokes the send routine

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::sendFPsiToGSP() {
//===================================================================
// Perform the transpose and then the blast off the final 1D-FFT

  int nchareG            = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
  int **tranpack         = scProxy.ckLocalBranch()->cpcharmParaInfo->index_tran_upack;
  int *nlines_per_chareG = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_per_chareG;
  complex *vks_on_state  = rs.planeArr;

 //------------------------------------------------------------------
  if (config.useMssInsGP){mssInstance.beginIteration();}

    for (int ic = 0; ic < nchareG; ic ++) { // chare arrays to which we will send

      int sendFFTDataSize = nlines_per_chareG[ic];
      GSIFFTMsg *msg      = new (sendFFTDataSize, 8 * sizeof(int)) GSIFFTMsg; 
      msg->size           = sendFFTDataSize;
      msg->offset         = thisIndex.y;    // z-index
      complex *data       = msg->data;

      if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.gsifftpriority+thisIndex.x*rs.planeSize[0];
      }//endif

      for(int i=0;i<sendFFTDataSize;i++){data[i] = vks_on_state[tranpack[ic][i]];}
      gproxy(thisIndex.x, ic).doIFFT(msg); // send the message
      CmiNetworkProgress();

    }//end for : chare sending

  if (config.useMssInsGP){mssInstance.endIteration();}
 //------------------------------------------------------------------

//===================================================================
// clean up the states

  if(config.conserveMemory){
     rs.destroy();
  }//endif

//===================================================================
//  CkPrintf("RSP [%d %d] sent back to GSP \n",thisIndex.x,thisIndex.y);
//============================================================================
   }//end routine : CP_State_RealSpacePlane::doProduct()
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
* Setting up the multicast trees for Gengbin's library 
*/
//============================================================================
void CP_State_RealSpacePlane::init(ProductMsg *msg){
    int i=1; 
    CkGetSectionInfo(cookie, msg);
    contribute(sizeof(int), &i, CkReduction::sum_int, 
	       CkCallback(CkIndex_main::doneInit(NULL),mainProxy));
    // do not delete nokeep message
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::ResumeFromSync(){
    if(config.useMssInsGP)
	ComlibResetProxy(&gproxy);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::printData() {
    char str[20];
    sprintf(str, "RSP%dx%d.out", thisIndex.x, thisIndex.y);
    FILE *outfile = fopen(str, "w");
    int ioff = 0;
    for (int i = 0; i < rs.planeSize[0]; i++) {
      for (int j = 0; j < rs.planeSize[1]; j++){
         fprintf(outfile, "%10.5lf", rs.planeArr[(j+ioff)].getMagSqr());
         fprintf(outfile, "\n");
      }
      ioff += rs.planeSize[0];
    }
    fclose(outfile);
}
//============================================================================
