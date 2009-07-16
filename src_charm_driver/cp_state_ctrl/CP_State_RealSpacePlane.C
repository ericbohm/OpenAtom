// Things to do
// dofft should do the fft
//#define RSVKS_BARRIER  
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

#include "CP_State_GSpacePlane.h"
#include "CP_State_Plane.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "main/groups.h"
#include "main/cpaimd.h"
#include "utility/util.h"
#include "charm++.h"

#include <iostream>
#include <fstream>
#include <cmath>
    using namespace std;



//============================================================================
extern CProxy_TimeKeeper                      TimeKeeperProxy;
extern CkVec <CProxy_AtomsGrp>                UatomsGrpProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>    UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>            UgSpaceDriverProxy;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>   UrhoRealProxy;
extern CkVec <CProxy_CP_State_RealSpacePlane> UrealSpacePlaneProxy;
extern CProxy_CPcharmParaInfoGrp              scProxy;
extern CProxy_main                            mainProxy;
extern CkVec <CProxy_CP_State_ParticlePlane>  UparticlePlaneProxy;
extern CkVec <CProxy_FFTcache>                UfftCacheProxy;
extern CProxy_InstanceController              instControllerProxy;
extern CkGroupID            mCastGrpId;
extern ComlibInstanceHandle mssInstance;

extern int    sizeX;
extern Config config;
extern CkReduction::reducerType sumFastDoubleType;
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

#ifdef RSVKS_BARRIER  // pause for every single chare to finish
    if(!(c->allVksDone())){
      RTH_Suspend(); // wait for broadcast that all vks is done  
    }//endif
#endif               //end pause
    c->thisProxy(c->thisIndex.x,c->thisIndex.y).doVksFFT(); // vks(r) arrives in doproduct(msg) which resumes
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
CP_State_RealSpacePlane::CP_State_RealSpacePlane( int gSpaceUnits, 
                  int realSpaceUnits, int _ngrida, int _ngridb, int _ngridc,
		  int _rfortime, int _rbacktime, UberCollection _instance)
  : thisInstance(_instance)
{
//============================================================================
//  ckout << "State R Space Constructor : "
//	<< thisIndex.x << " " << thisIndex.y << " " <<CkMyPe() << endl;
//============================================================================

    countProduct=0;
    count = 0;
    vksDone=false;
    rhoRsubplanes = config.rhoRsubplanes;
    numCookies=0;
    ngrida = _ngrida;
    ngridb = _ngridb;
    ngridc = _ngridc;
    csize = (ngrida/2 + 1)*ngridb; 
    rsize = (ngrida   + 2)*ngridb; ;
    iplane_ind  = thisIndex.y;
    forwardTimeKeep=_rfortime;
    backwardTimeKeep=_rbacktime;
    initRealStateSlab(&rs, ngrida, ngridb, ngridc, gSpaceUnits, 
                       realSpaceUnits, thisIndex.x, thisIndex.y);


    setMigratable(false);
    cookie= new CkSectionInfo[rhoRsubplanes];
    iteration = 0;
    run();

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);

  p|iplane_ind;
  p|iteration;
  p|rhoRsubplanes;
  p|ngrida;
  p|ngridb;
  p|ngridc;
  p|count;
  p|countProduct;
  p|csize;
  p|rsize;
  PUParray(p,cookie,rhoRsubplanes);
  p|gproxy;
  p|numCookies;
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
void CP_State_RealSpacePlane::acceptFFT(RSFFTMsg *msg) {
//============================================================================
#ifdef _CP_SUBSTEP_TIMING_
  if(forwardTimeKeep>0)
    {
      double rstart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&rstart,CkReduction::min_double, cb , forwardTimeKeep);
    }
#endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
    }
#endif

#ifdef RSVKS_BARRIER
  //reset each iteration
  vksDone=false;
#endif

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
      iteration++;
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

    FFTcache *fftcache  = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    int nplane_x        = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
    complex *planeArr   = rs.planeArr;
    double  *planeArrR  = rs.planeArrR;

#ifndef CMK_OPTIMIZE    
    double StartTime=CmiWallTimer();
#endif

    fftcache->doStpFFTGtoR_Rchare(planeArr,planeArrR,nplane_x,ngrida,ngridb,iplane_ind);

    fftcache->getCacheMem("CP_State_RealSpacePlane::doFFT");
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
  int ees_nonlocal = scProxy.ckLocalBranch()->cpcharmParaInfo->ees_nloc_on;
  if(ees_nonlocal==1 && config.launchNLeesFromRho==0){
    /*    if(config.nchareG>32)
      {
	fprintf(stderr,"[%d,%d]on %d nchareG is %d\n",thisIndex.x, thisIndex.y, CkMyPe(),config.nchareG);
	} */
    //    CkAssert(config.nchareG<=ngridc);
    int div    = (config.nchareG / ngridc);
    int rem    = (config.nchareG % ngridc);
    int ind    = thisIndex.y+ngridc;    
    if(div>1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("NchareG too large for ngridc for launch Scheme\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

    if(thisIndex.y<config.nchareG){
      UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y).startNonLocalEes(iteration);
    }//endif
    if((div==1) && (thisIndex.y<rem))
      {
	UgSpaceDriverProxy[thisInstance.proxyOffset](thisIndex.x,ind).startNonLocalEes(iteration);
      }
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
  bzero(planeArrR,(ngrida+2)*ngridb*sizeof(double));
  fftcache->freeCacheMem("CP_State_RealSpacePlane::doFFT");
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
// Grab some local pointer and do some checking.

   CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch(); 
   FFTcache *fftcache  = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
   double *data        = fftcache->tmpDataR;

#ifdef _CP_DEBUG_STATER_VERBOSE_
    CkPrintf("In StateRSpacePlane[%d %d] doReduction %d\n", thisIndex.x, thisIndex.y,
                CmiMemoryUsage());
#endif

#ifdef _NAN_CHECK_
  for(int i=0;i<ngrida*ngridb ;i++){
      if(isnan(data[i])!=0){
	CkPrintf("RS [%d %d] issuing nan at %d out of %d\n",
           thisIndex.x, thisIndex.y, i, ngridb*ngrida);
	CkAbort("RS nan in the fftcache");
      }
  }//endif
#endif

//============================================================================
// Perform the Reduction to get the density : vks holds psi^2 for us
// Return values for vks cannot hit this chare until the reduciton is complete.

#ifndef CMK_OPTIMIZE    
   double StartTime=CmiWallTimer();
#endif
  
  // Need loop of contribute calls, one for each nchareRhoSplit offset 
  // into data. Need a vector of cookies and callback functions

  int subSize = (ngridb/rhoRsubplanes);
  int subRem  = (ngridb % rhoRsubplanes);

  int off = 0;
  for(int subplane=0; subplane<rhoRsubplanes; subplane++){
    int dataSize = subSize*ngrida;
    if(subplane < subRem){dataSize += ngrida;}
    CkCallback cb(CkIndex_CP_Rho_RealSpacePlane::acceptDensity(0),
                  CkArrayIndex2D(thisIndex.y,subplane),UrhoRealProxy[thisInstance.proxyOffset].ckGetArrayID());
    mcastGrp->contribute(dataSize*sizeof(double),&(data[off]),
                         sumFastDoubleType,cookie[subplane],cb);
    off += dataSize;
  }//endfor : subplanes
  CkAssert(off==ngrida*ngridb);

  CmiNetworkProgress(); // yuck!

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(DoFFTContribute_, StartTime, CmiWallTimer());
#endif    

  fftcache->freeCacheMem("CP_State_RealSpacePlane::doReduction");
#ifdef _CP_SUBSTEP_TIMING_
  if(forwardTimeKeep>0)
    {
      double rend=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&rend,CkReduction::max_double, cb , forwardTimeKeep);
    }
#endif

//============================================================================
    }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *   In this method, we receive vks from the density. We apply this to psi to
 *   create psiVks, when all are here we call the working doProduct for
 *   FFTing.  This is a stream processing scheme.
 */
//============================================================================
void CP_State_RealSpacePlane::acceptProduct(ProductMsg *msg) {
//============================================================================
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
 
#ifdef _CP_DEBUG_STATER_VERBOSE_
  CkPrintf("In StateRSpacePlane[%d %d] doProd \n", thisIndex.x, thisIndex.y);
#endif
#ifdef _CP_SUBSTEP_TIMING_
  if(backwardTimeKeep>0)
    {
      double rstart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&rstart,CkReduction::min_double, cb , backwardTimeKeep);
    }
#endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->datalen ;i++){
      CkAssert(isnan(msg->data[i])==0);
  }
#endif

//============================================================================
//Unpack and check size, use message data for multiply, then resume
//which calls doProduct


  double *vks_tmp = msg->data;
  int mychare    = msg->idx;
  int mysize     = msg->datalen;
  int myindex    = msg->subplane; // subplaneindex of rhor
  CkAssert(mychare==thisIndex.y);

  int subSize    = (ngridb/rhoRsubplanes);
  int subRem     = (ngridb % rhoRsubplanes);
  int subAdd     = (myindex < subRem ? 1 : 0);
  int subMax     = (myindex < subRem ? myindex : subRem);
  int off        = (subSize*myindex + subMax)*(ngrida+2);
  int myNgridb   = subSize+subAdd;
  int size       = ngrida*myNgridb;
  CkAssert(mysize == size);

//============================================================================	

  double *psiVks      = rs.planeArrR;

  // multiply psi by vks to form psiVks
  for(int i=0,i2=off;i<myNgridb;i++,i2+=2){
    for(int j=i*ngrida;j<(i+1)*ngrida;j++){
      psiVks[(j+i2)] *= vks_tmp[j];
    }//endfor
  }//endfor
  CmiNetworkProgress();

//============================================================================	

  countProduct++;
  if(countProduct==rhoRsubplanes){
    countProduct=0;
#ifdef RSVKS_BARRIER
    int wehaveours=1;
    contribute(sizeof(int),&wehaveours,CkReduction::sum_int,
	       CkCallback(CkIndex_CP_State_RealSpacePlane::rdoneVks(NULL),UrealSpacePlaneProxy[thisInstance.proxyOffset]));
#endif

    RTH_Runtime_resume(run_thread); // this is scalar, we continue right on
                                    // as threaded loops calls do vksfft
  }//endif

//---------------------------------------------------------------------------
  }//endroutine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/*
 *        Do the fft of psi*vks stored in rho_rs.planeArrayR
 */
//============================================================================
void CP_State_RealSpacePlane::doVksFFT() {
//============================================================================
// A little output under some circumstances

#ifdef _CP_DEBUG_STATER_VERBOSE_
  CkPrintf("In RealSpacePlane[%d %d] doProduct %d\n",
             thisIndex.x, thisIndex.y,CmiMemoryUsage());
#endif

#ifndef _CP_DEBUG_RHO_OFF_  
#ifdef _CP_DEBUG_VKS_RSPACE_
  if(thisIndex.x==0 && thisIndex.y == 0){
     FILE *fp = fopen("psivks_real_y0_state0.out","w");
      for(int i=0;i<(ngrida+2)*ngridb;i++){
        fprintf(fp,"%g\n",rho_rs.planeArrayR[i]);
      }//endfor
     fclose(fp);
  }//endif
#endif    
#endif

//===================================================================
// Do the FFT of psi*vks

 //------------------------------------------------------------------
 // Start the timer
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

 //------------------------------------------------------------------
 // The FFT
  FFTcache *fftcache  = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  int nplane_x        = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
  complex *planeArr   = rs.planeArr;
  double *planeArrR   = rs.planeArrR;
  fftcache->doStpFFTRtoG_Rchare(planeArr,planeArrR,nplane_x,ngrida,ngridb,iplane_ind);

 //------------------------------------------------------------------
 // End timer 
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(doRealBwFFT_, StartTime, CmiWallTimer());
#endif

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

#ifdef USE_COMLIB
#ifdef OLD_COMMLIB
  if (config.useMssInsGP){mssInstance.beginIteration();}
#else
  //  if (config.useMssInsGP){ComlibBegin(gproxy);}
#endif
#endif

    for (int ic = 0; ic < nchareG; ic ++) { // chare arrays to which we will send

      int sendFFTDataSize = nlines_per_chareG[ic];
      GSIFFTMsg *msg      = new (sendFFTDataSize, 8 * sizeof(int)) GSIFFTMsg; 
      msg->size           = sendFFTDataSize;
      msg->offset         = thisIndex.y;    // z-index
      complex *data       = msg->data;

      if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.gsifftpriority+thisIndex.x*rs.ngrid_a;
      }//endif

      for(int i=0;i<sendFFTDataSize;i++){data[i] = vks_on_state[tranpack[ic][i]];}
      gproxy(thisIndex.x, ic).acceptIFFT(msg); // send the message
      CmiNetworkProgress();

    }//end for : chare sending

#ifdef USE_COMLIB
#ifdef OLD_COMMLIB
  if (config.useMssInsGP){mssInstance.endIteration();}
#else
  //  if (config.useMssInsGP){ComlibEnd(gproxy);}
#endif
#endif

 //------------------------------------------------------------------

//===================================================================
// clean up the states

  if(config.conserveMemory){
     rs.destroy();
  }//endif
#ifdef _CP_SUBSTEP_TIMING_
  if(backwardTimeKeep>0)
    {
      double rend=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&rend,CkReduction::max_double, cb , backwardTimeKeep);
    }
#endif

//============================================================================
   }//end routine : CP_State_RealSpacePlane::sendFPsiToGSP
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
* Setting up the multicast trees for Gengbin's library 
*/
//============================================================================
void CP_State_RealSpacePlane::init(ProductMsg *msg){
//============================================================================
// Do not delete msg. Its a nokeep.
//============================================================================
    gproxy = UgSpacePlaneProxy[thisInstance.proxyOffset];
    numCookies++;
    // based on where this came from, put it in the cookie vector
    CkGetSectionInfo(cookie[msg->subplane], msg);
    if(numCookies == config.rhoRsubplanes)
      contribute(sizeof(int), &numCookies, CkReduction::sum_int, 
	       CkCallback(CkIndex_InstanceController::doneInit(NULL),CkArrayIndex1D(thisInstance.proxyOffset),instControllerProxy), thisInstance.proxyOffset);
    // do not delete nokeep message
#ifdef USE_COMLIB
    if (config.useMssInsGP){
      ComlibAssociateProxy(&mssInstance,gproxy);
    }//endif
#endif

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::ResumeFromSync(){
#ifdef USE_COMLIB
    if(config.useMssInsGP)
	ComlibResetProxy(&gproxy);
#endif
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
    for (int i = 0; i < rs.ngrid_b; i++) {
      for (int j = 0; j < rs.ngrid_a; j++){
         fprintf(outfile, "%10.5lf", rs.planeArr[(j+ioff)].getMagSqr());
         fprintf(outfile, "\n");
      }
      ioff += rs.ngrid_a;
    }
    fclose(outfile);
}
//============================================================================



//============================================================================
//     All Rho RS objects have finished vks : Debugging only
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_State_RealSpacePlane::rdoneVks(CkReductionMsg *msg){
      delete msg;
      //let my ffts go!
      vksDone=true;
      RTH_Runtime_resume(run_thread);
  }
//============================================================================
