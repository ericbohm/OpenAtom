//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_RHartExc.C
 *
 *  This is a description of the "life" of a CP_Rho_RHartExc  object
 *
 * FIll in details here
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
#include "eesCache.h"
#include "CP_State_Plane.h"

#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CPcharmParaInfoGrp    scProxy;
extern CProxy_AtomsGrp              atomsGrpProxy;
extern CProxy_CP_Rho_RHartExt       rhoRHartExtProxy;
extern CProxy_CP_Rho_GHartExt       rhoGHartExtProxy;
extern CProxy_eesCache              eesCacheProxy;
extern CProxy_FFTcache              fftCacheProxy;

extern ComlibInstanceHandle commRHartGHartIns;

extern int    sizeX;
extern Config config;


#define _CP_RHART_VERBOSE_OFF_
//#define _CP_RHART_VERBOSE_


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *  This object just gets a rho message, computes GHartExt, and sends
 *  vks.  
 */
//============================================================================
CP_Rho_RHartExt::CP_Rho_RHartExt(int _ngrida, int _ngridb, int _ngridc, 
                                 int _ees_eext_on, int _natmTyp){
//============================================================================
// Set some pareameters

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  nplane_rho_x         = sim->nplane_rho_x;
  rhoRsubplanes        = config.rhoRsubplanes;

  CkAssert(nplane_rho_x >= rhoRsubplanes); // safety : should already be checked.

  ngrida      = _ngrida; 
  ngridb      = _ngridb; 
  ngridc      = _ngridc; 
  ees_eext_on = _ees_eext_on;
  natmTyp     = _natmTyp;

//============================================================================
// Initialize some variables

  countIntRtoG     = 0;
  countIntGtoR[0]  = 0;
  countIntGtoR[1]  = 0;
  countFFT[0]      = 0;
  countFFT[1]      = 0;
  nAtmTypRecv      = 0;
  registrationFlag = 0;
  launchFlag       = 0;
  iteration        = 0;
  iterAtmTyp       = 0;
  countDebug       = 0;

  // Parallelization of rho(x,y,z) by (y,z)
  int div      = (ngridb/rhoRsubplanes); 
  int rem      = (ngridb % rhoRsubplanes);
  int max      = (thisIndex.y < rem ? thisIndex.y : rem);
  myNgridb     = (thisIndex.y<rem ? div+1 : div);  // number of y values/lines of x
  myBoff       = div*thisIndex.y + max;            // offset into y 
  nptsB        =  ngrida*myNgridb;                 // size of plane without extra room
  nptsExpndB   = (ngrida+2)*myNgridb;              // extra memory for RealToComplex FFT

  // Parallelization after transpose : rho(gx,y,z) : parallelize gx and z
  int divb     = (nplane_rho_x/rhoRsubplanes);
  int remb     = (nplane_rho_x % rhoRsubplanes);
  int maxb     = (thisIndex.y < remb ? thisIndex.y : remb);
  myNplane_rho = (thisIndex.y<remb ? divb+1 : divb);// number of x values/lines of y
  myAoff       = divb*thisIndex.y + maxb;
  nptsA        = 2*myNplane_rho*ngridb;            // memory size for fft in doubles
  nptsExpndA   = 2*myNplane_rho*ngridb;            // memory size for fft in doubles

  csize        = 0;
  csizeInt     = 0;
  if(ees_eext_on==1){
    csize        = (ngrida/2+1)*myNgridb;  // complex variable size

 
    if(rhoRsubplanes==1){
      eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
      eesData->registerCacheRHart(thisIndex.x);
    }else{
      CkPrintf("rhart registration broken\n"); CkExit();
    }

    atmSFC      = (complex*) fftw_malloc(csize*sizeof(complex));
    atmSFR      = reinterpret_cast<double*> (atmSFC);
    atmForcC    = atmSFC;
    atmForcR    = atmSFR;

    atmEwdSFC   = (complex*) fftw_malloc(csize*sizeof(complex));
    atmEwdSFR   = reinterpret_cast<double*> (atmEwdSFC);
    atmEwdForcC = atmEwdSFC;
    atmEwdForcR = atmEwdSFR;

    if(rhoRsubplanes>1){
      csizeInt     = nptsA/2;
      atmSFCint      = (complex*) fftw_malloc(csizeInt*sizeof(complex));
      atmSFRint      = reinterpret_cast<double*> (atmSFC);
      atmForcCint    = atmSFCint;
      atmForcRint    = atmSFRint;

      atmEwdSFCint   = (complex*) fftw_malloc(csizeInt*sizeof(complex));
      atmEwdSFRint   = reinterpret_cast<double*> (atmEwdSFC);
      atmEwdForcCint = atmEwdSFCint;
      atmEwdForcRint = atmEwdSFRint;
    }//endif

    int i=1;
    CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(NULL),rhoRHartExtProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

  rhoGHartProxy_com    = rhoGHartExtProxy;;
  if (config.useRHartInsGHart){
    ComlibAssociateProxy(&commRHartGHartIns,rhoGHartProxy_com);          
  }//endif

  if(ees_eext_on==1 && rhoRsubplanes>1){
    if(thisIndex.y==0&&thisIndex.x==0){
       CkPrintf("Warning rhorhart: ees_eext_on under development\n");
     }//endif
  }//endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_RHartExt::~CP_Rho_RHartExt(){
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::pup(PUP::er &p){
//============================================================================

  ArrayElement2D::pup(p);

   p|countDebug;
   p|nplane_rho_x;
   p|rhoRsubplanes;
   p|registrationFlag;
   p|launchFlag;
   p|ngrida;
   p|ngridb;
   p|ngridc;
   p|ees_eext_on;
   p|natmTyp;
   p|nAtmTypRecv;
   p|countIntRtoG;
   p(countIntGtoR,2);
   p(countFFT,2);
   p|iteration;
   p|iterAtmTyp;
   p|csize;
   p|csizeInt;
   p|myNgridb;
   p|myBoff;
   p|nptsB;
   p|nptsExpndB;
   p|myNplane_rho;
   p|myAoff;
   p|nptsA;
   p|nptsExpndA;

   if (p.isUnpacking() && ees_eext_on==1) {
     atmSFC      = (complex*) fftw_malloc(csize*sizeof(complex));
     atmSFR      = reinterpret_cast<double*> (atmSFC);
     atmForcC    = atmSFC;
     atmForcR    = atmSFR;
     atmEwdSFC   = (complex*) fftw_malloc(csize*sizeof(complex));
     atmEwdSFR   = reinterpret_cast<double*> (atmEwdSFC);
     atmEwdForcC = atmEwdSFC;
     atmEwdForcR = atmEwdSFR;
     if(rhoRsubplanes>1){
       atmSFCint      = (complex*) fftw_malloc(csizeInt*sizeof(complex));
       atmSFRint      = reinterpret_cast<double*> (atmSFCint);
       atmForcCint    = atmSFC;
       atmForcRint    = atmSFR;
       atmEwdSFCint   = (complex*) fftw_malloc(csizeInt*sizeof(complex));
       atmEwdSFRint   = reinterpret_cast<double*> (atmEwdSFCint);
       atmEwdForcCint = atmEwdSFCint;
       atmEwdForcRint = atmEwdSFRint;
     }//endif
   }//endif

   if(ees_eext_on==1) {
     p((char*)atmSFC,csize*sizeof(complex));
     p((char*)atmEwdSFC,csize*sizeof(complex));
     if(rhoRsubplanes>1){
       p((char*)atmSFCint,csizeInt*sizeof(complex));
       p((char*)atmEwdSFCint,csizeInt*sizeof(complex));
     }//endif
   }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// Invoke by Rspace-density : Density has arrived in r-space and will soon arrive
//                            in g-space. Get moving RhartExt
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::startEextIter(){
//============================================================================
// Check for error

  if(atomsGrpProxy.ckLocalBranch()->iteration != iteration){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Flow of Control Error in RHartExtVks : atoms slow\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// This is a new time step : Increment time step counter, zero atm typ counter
//                           Launch if we are ready

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in start : %d\n",thisIndex.x,registrationFlag);
#endif

  if(iterAtmTyp==natmTyp){CkPrintf("%d signing off\n",thisIndex.x);CkExit();}

  iteration ++;
  iterAtmTyp  = 0; 
  nAtmTypRecv = 0;

  launchFlag=1;
  if(registrationFlag==1){computeAtmSF();}

//============================================================================
  }//end routine
//============================================================================


//==========================================================================
// Make sure everyone is registered on the 1st time step
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void CP_Rho_RHartExt::registrationDone(CkReductionMsg *msg) {
//==========================================================================

  int sum = ((int *)msg->getData())[0];
  delete msg;

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in reg : %d\n",thisIndex.x,sum);
#endif

  // use barrier as a debug tool
  if(iterAtmTyp==natmTyp){
    CkPrintf("HI, I am rhart chare %d exiting at %d\n",thisIndex.x,iterAtmTyp+1);
    CkExit();
  }//endif

  registrationFlag=1;
  if(launchFlag==1){computeAtmSF();}

  }
//==========================================================================



//============================================================================
// Start the real space part of the EES interpolation for atmSF(iatmTyp)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::computeAtmSF(){
//============================================================================

  iterAtmTyp++;
  if(iterAtmTyp>natmTyp){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Too many iterations RHartExtVks\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in computeatmsf at iter %d\n",thisIndex.x,iterAtmTyp);
#endif

//============================================================================
// Look up some constants and fill the r-space grid

  int myPlane      =  thisIndex.x;
  int mySubPlane   =  thisIndex.y;
  int itime        = iteration;

  eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
  if(iterAtmTyp==1){
    if(rhoRsubplanes==1){
      eesData->queryCacheRHart(myPlane,itime,iterAtmTyp);
    }else{
      CkPrintf("rhart query broken\n"); CkExit();
    }//endif
  }

  int *plane_index = eesData->RhoRHartData[myPlane].plane_index;
  int **igrid      = eesData->RhoRHartData[myPlane].igrid;
  double **mn      = eesData->RhoRHartData[myPlane].mn;

  AtomsGrp *ag     = atomsGrpProxy.ckLocalBranch(); // find me the local copy
  int natm         = ag->natm;

  if(rhoRsubplanes==1){
    CPLOCAL::eesPackGridRchare(natm,iterAtmTyp,atmSFR,myPlane,igrid,mn,plane_index);
  }else{
    CkPrintf("rhart pack broken\n"); CkExit();
  }//endif

//============================================================================
// FFT the result to G-space

  fftAtmSfRtoG();

//============================================================================
  }//end routine 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// FFT of SF(x,y,z,iatmTyp) -> SF(gx,gy,z,iatmTyp)
//============================================================================
void CP_Rho_RHartExt::fftAtmSfRtoG(){
//============================================================================
#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in fftsback at %d\n",thisIndex.x,iterAtmTyp);
#endif

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  

  if(rhoRsubplanes>1){
    fftcache->doEextFFTRxToGx_Rchare(atmSFC,atmSFR,nplane_rho_x,ngrida,myNgridb);
  }else{
    fftcache->doEextFFTRtoG_Rchare(atmSFC,atmSFR,nplane_rho_x,ngrida,ngridb);
  }//endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(doEextFFTRtoG_, StartTime, CmiWallTimer());    
#endif

  if(rhoRsubplanes>1){
    sendAtmSfRyToGy();   // double transpose method (yz ---> gx,z)
  }else{
    sendAtmSfRhoGHart();      // single transpose method (z ---> gx,gy)
  }//endif

//============================================================================
  }// end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// Double Transpose Fwd Send : A(gx,y,z) on the way to A(gx,gy,z)
//                             Send so that (y,z) parallelism is 
//                             switched to (gx,z)
//
//============================================================================
void CP_Rho_RHartExt::sendAtmSfRyToGy(){
//============================================================================
// Launch the communication

   CkAssert(rhoRsubplanes>1);
  //-----------------------------------------------------------------------------
  // Commlib launch : 

#ifdef _ERIC_SETS_UP_COMMLIB_
    if(config.useRInsRhoRP)    commRealInstanceRx.beginIteration();
#endif

   //-----------------------------------------------------------------------------
   // Send the data : I have myNgridB values of y  (gx,y) y=1...myNgridB and all gx
   //                 Send all the `y' I have for the gx range desired after transpose

    int stride = ngrida/2+1;
    int ix     = thisIndex.x;
    for(int ic = 0; ic < rhoRsubplanes; ic ++) { // chare arrays to which we will send

      int div     = (nplane_rho_x/rhoRsubplanes);   //parallelize gx
      int rem     = (nplane_rho_x % rhoRsubplanes);
      int add     = (ic < rem ? 1 : 0);
      int max     = (ic < rem ? ic : rem);
      int ist     = div*ic + max;        // start of gx desired by chare ic
      int iend    = ist + div + add;     // end   of gx desired by chare ic
      int size    = (iend-ist)*myNgridb; // data size

      int sendFFTDataSize = size;
      RhoGHartMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) RhoGHartMsg; 
      msg->size        = size;
      msg->offset      = myBoff;      // where the myNgridB y-lines start.
      msg->num         = myNgridb;    // number of y-lines I have. 
      msg->iter        = iterAtmTyp;
      complex *data    = msg->data;   // data

      if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.y;
      }//endif

      for(int i=ist,koff=0;i<iend;i++,koff+=myNgridb){
        for(int k=koff,ii=i;k<myNgridb+koff;k++,ii+=stride){
          data[k] = atmSFC[ii]; 
	}//endfor
      }//endfor

      rhoRHartExtProxy(ix,ic).recvAtmSfRyToGy(msg);

#ifdef _ERIC_SETS_UP_COMMLIB_
      rhoRHartExtProxy(ix,ic).recvAtmSfRyToGy(msg);
#endif

#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
    }//end for : chare sending

  //-----------------------------------------------------------------------------
  // Commlib stop

#ifdef _ERIC_SETS_UP_COMMLIB_
       if(config.useRInsRhoRP)    commRealInstanceRx.endIteration();
#endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// Double Transpose Fwd Recv : A(gx,y,z) on the way to A(gx,gy,z)
//                             Recv so that (y,z) parallel switched to (gx,z)
//
// Invoked natm_typ times per algorithm step : 
//
//============================================================================
void CP_Rho_RHartExt::recvAtmSfRyToGy(RhoGHartMsg *msg){
//============================================================================

  int size         = msg->size;  // msg size
  int iter         = msg->iter;  //
  int num          = msg->num;
  int offset       = msg->offset;
  complex *msgData = msg->data;

  CkAssert(size==myNplane_rho*num);
  CkAssert(rhoRsubplanes>1);
  CkAssert(iter==iterAtmTyp);

//============================================================================

  for(int js=0,j=offset;js<size;js+=num,j+=ngridb){
   for(int is=js,i=j;is<num+js;is++,i++){
     atmSFCint[i] = msgData[is];
   }//endfor
  }//endfor

  delete msg;

//============================================================================
// Do the Y fft, invoke communication

  countIntRtoG++;
  if(countIntRtoG==rhoRsubplanes){
    countIntRtoG = 0;
    FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
    fftcache->doEextFFTRyToGy_Rchare(atmSFCint,atmSFRint,nplane_rho_x,ngrida,myNgridb);
    sendAtmSfRhoGHart(); 
  }//endfor

//============================================================================
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Send SF(gx,gy,z,iatmTYP) to g-space whence the FFT will be completed
//============================================================================
void CP_Rho_RHartExt::sendAtmSfRhoGHart(){
//============================================================================

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in sendsf at %d\n",thisIndex.x,iterAtmTyp);
#endif

  CPcharmParaInfo *sim      = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int nchareG               = sim->nchareRhoGEext;

  int **tranpack            = sim->index_tran_upack_eext;
  int *nlines_per_chareG    = sim->nlines_per_chareRhoGEext;

  int ***tranupack_Y        = sim->index_tran_upack_eext_y;
  int **nlines_per_chareG_Y = sim->nline_send_eext_y;

//===================================================================
// Perform the transpose and then the blast off the final 1D-FFT

 //----------------------------------------------------------------
 // start commlib
  if(rhoRsubplanes==1){
     if(config.useRHartInsGHart){commRHartGHartIns.beginIteration();}
  }//endif

 //----------------------------------------------------------------
 // do the send

    int iy = thisIndex.y;
    for (int ic=0;ic<nchareG;ic++) { // chare arrays to which we will send

     //---------------------
     // malloc the message
      int sendFFTDataSize = nlines_per_chareG[ic];
      if(rhoRsubplanes!=1){sendFFTDataSize = nlines_per_chareG_Y[ic][iy];}
      RhoGHartMsg *msg    = new (sendFFTDataSize, 8 * sizeof(int)) RhoGHartMsg; 

     //----------------------
     // pack the message
      if(config.prioEextFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.rhogHartpriority+thisIndex.x*10;
      }//endif
      msg->iter           = iterAtmTyp;
      msg->size           = sendFFTDataSize;
      msg->offset         = thisIndex.x;    // c-plane-index
      msg->offsetGx       = thisIndex.y;    // gx parallelization index
      complex *data       = msg->data;
      if(rhoRsubplanes==1){      
        for(int i=0;i<sendFFTDataSize;i++){data[i] = atmSFC[tranpack[ic][i]];}
      }else{
        for(int i=0;i<sendFFTDataSize;i++){
          data[i] = atmSFCint[tranupack_Y[ic][iy][i]];
        }//endfor
      }//endif

     //-----------------
     // Send the message
      if(rhoRsubplanes==1){
        rhoGHartProxy_com(ic,0).recvAtmSFFromRhoRHart(msg); // send the message
      }else{
        rhoGHartExtProxy(ic,iy).recvAtmSFFromRhoRHart(msg); // send the message
      }//endif

#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
    }//end for : chare sending

 //----------------------------------------------------------------
 // stop commlib

  if(rhoRsubplanes==1){
    if(config.useRHartInsGHart){commRHartGHartIns.endIteration();}
  }//endif

//============================================================================
  }// end routine
//============================================================================


//============================================================================
// Hartree sends back atom forces from e-atm interation
// Depending on the flag, it is Ewald or e-atm interation
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::recvAtmForcFromRhoGHart(RhoRHartMsg *msg){
//============================================================================

  CPcharmParaInfo *sim      = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG               = sim->nchareRhoGEext;

  int **tranUnpack          = sim->index_tran_upack_eext;
  int *nlines_per_chareG    = sim->nlines_per_chareRhoGEext;

  int ***tranupack_Y        = sim->index_tran_upack_eext_y;
  int **nlines_per_chareG_Y = sim->nline_send_eext_y;
  int iy                    = thisIndex.y;
   
  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  int iter               = msg->iter;
  complex *partiallyFFTd = msg->data;

  int mySize;
  int csizenow;
  if(rhoRsubplanes==1){
     mySize   = nlines_per_chareG[Index];
     csizenow = csize;
  }else{
     mySize   = nlines_per_chareG_Y[Index][iy];
     csizenow = csizeInt;
  }//endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
  }
#endif

//============================================================================
// Perform some error checking

  countFFT[iopt]++;
  if (countFFT[iopt] > nchareG) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               countFFT[iopt],nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=mySize){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.1, %d != %d for rho rhart chare %d %d\n",size,mySize,
                  thisIndex.x,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(iter!=iterAtmTyp){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.2, %d != %d for rho rhart chare %d %d\n",iter,iterAtmTyp,
                  thisIndex.x,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message

  complex *data;
  if(rhoRsubplanes==1){
    switch(iopt){
      case 0 : data = atmSFC; break;
      case 1 : data = atmEwdSFC; break;
      default: CkAbort("impossible iopt");      break;
    }//end switch
  }else{
    switch(iopt){
      case 0 : data = atmSFCint; break;
      case 1 : data = atmEwdSFCint; break;
      default: CkAbort("impossible iopt");      break;
    }//end switch
  }//endif

  if(countFFT[iopt]==1){memset(data,0,sizeof(complex)*csizenow);}

  if(rhoRsubplanes==1){
    for(int i=0;i<size;i++){data[tranUnpack[Index][i]] = partiallyFFTd[i];}
  }else{
    for(int i=0;i<size;i++){
      data[tranupack_Y[Index][iy][i]] = partiallyFFTd[i];
    }//endfor
  }//endif

  delete msg;

//============================================================================
// When you have all the data : finish the FFT back to real space

  if (countFFT[iopt] == nchareG){
#ifdef _CP_RHART_VERBOSE_
    CkPrintf("HI, I am rhart chare %d in recvsf with opt %d at %d\n",
               thisIndex.x,iopt,iterAtmTyp);
#endif
    countFFT[iopt]=0;
    nAtmTypRecv++;
    fftAtmForcGtoR(iopt);
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Complete the FFT : sf(gx,gy,z) ->sf(x,y,z)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::fftAtmForcGtoR(int flagEwd){
//============================================================================

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in FFTAtmForc w %d at %d \n",
        thisIndex.x,flagEwd,iterAtmTyp);
#endif

  double *dataR;
  complex *dataC;  
  if(rhoRsubplanes==1){
    switch(flagEwd){
      case 0: dataR=atmSFR;    dataC=atmSFC;    break;
      case 1: dataR=atmEwdSFR; dataC=atmEwdSFC; break;
      default: CkAbort("impossible iopt");      break;
    }//endif
  }else{
    switch(flagEwd){
      case 0: dataR=atmSFRint;    dataC=atmSFCint;    break;
      case 1: dataR=atmEwdSFRint; dataC=atmEwdSFCint; break;
      default: CkAbort("impossible iopt");            break;
    }//endif
  }//endif

//============================================================================

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
  if(rhoRsubplanes==1){  
    fftcache->doEextFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb);
    computeAtmForc(flagEwd);
  }else{
    fftcache->doEextFFTGyToRy_Rchare(dataC,dataR,myNplane_rho,ngrida,ngridb);
    sendAtmForcGxToRx(flagEwd);
  }//endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(doEextFFTGtoR_, StartTime, CmiWallTimer());    
#endif

//============================================================================ 
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::sendAtmForcGxToRx(int iopt){
//============================================================================
  CkAssert(rhoRsubplanes>1);

  complex *FFTresult;
  switch(iopt){   
    case 0 : FFTresult = atmSFCint;      break;
    case 1 : FFTresult = atmEwdSFCint;   break;
    default: CkAbort("impossible iopt"); break;
  }//end switch

//============================================================================
// Launch the communication

  //-----------------------------------------------------------------------------
  // Commlib launch : 

#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){   
       case 0 : if(config.useRInsRhoRP)    commRealInstanceRx.beginIteration();    break;
       case 1 : if(config.useRInsIGXRhoRP) commRealIGXInstanceRx.beginIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch
#endif

   //-----------------------------------------------------------------------------
   // Send the data : I have myNgridB values of y  (gx,y) y=1...myNgridB and all gx
   //                 Send all the `y' I have for the gx range desired after transpose

    int ix = thisIndex.x;
    for(int ic = 0; ic < rhoRsubplanes; ic ++) { // chare arrays to which we will send
      int div     = (ngridb/rhoRsubplanes);   //parallelize y
      int rem     = (ngridb % rhoRsubplanes);
      int add     = (ic < rem ? 1 : 0);
      int max     = (ic < rem ? ic : rem);
      int ist     = div*ic + max;        // start of y desired by chare ic
      int iend    = ist + div + add;     // end   of y desired by chare ic
      int size    = (iend-ist)*myNplane_rho; // data size : send all the gx I have

      int sendFFTDataSize = size;
      RhoGHartMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) RhoGHartMsg; 
      msg->size        = size;
      msg->iopt        = iopt;
      msg->iter        = iterAtmTyp;
      msg->offset      = myAoff;       // where the gx-lines I have start.
      msg->num         = myNplane_rho; // number of gx-lines I have. 
      complex *data    = msg->data;    // data

      if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.y;
      }//endif

      for(int i=ist,koff=0;i<iend;i++,koff+=myNplane_rho){ 
        for(int k=koff,ii=i;k<myNplane_rho+koff;k++,ii+=ngridb){
          data[k] = FFTresult[ii]; 
	}//endfor
      }//endfor

      rhoRHartExtProxy(ix,ic).recvAtmForcGxToRx(msg);

#ifdef _ERIC_SETS_UP_COMMLIB_
      switch(iopt){
        case 0 : rhoGProxy_com(ic,0).acceptRhoGradVksGxToRx(msg);break;
        case 1 : rhoGProxyIGX_com(ic,0).acceptRhoGradVksGxToRx(msg); break;
        default: CkAbort("impossible iopt");break;
      }//end switch
#endif

#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
    }//end for : chare sending

  //-----------------------------------------------------------------------------
  // Commlib stop

#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){
       case 0 : if(config.useRInsRhoRP)    commRealInstanceRx.endIteration();    break;
       case 1 : if(config.useRInsIGXRhoRP) commRealIGXInstanceRx.endIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch
#endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::recvAtmForcGxToRx(RhoGHartMsg *msg){
//============================================================================

  int size         = msg->size;  // msg size
  int iopt         = msg->iopt;  
  int iter         = msg->iter;  
  int num          = msg->num;   // number of lines along `a' sent
  int offset       = msg->offset; // offset into lines along `a'
  complex *msgData = msg->data;

  CkAssert(size==myNgridb*num);
  CkAssert(rhoRsubplanes>1);
  CkAssert(iter!=iterAtmTyp);

  complex *dataC;
  double  *dataR;
  switch(iopt){
    case 0: dataR=atmSFR;    dataC=atmSFC;    break;
    case 1: dataR=atmEwdSFR; dataC=atmEwdSFC; break;
    default: CkAbort("impossible iopt");      break;
  }//end switch

//============================================================================
// Unpack the message 

  countIntGtoR[iopt]++;
  if(countIntGtoR[iopt]==1){bzero(dataR,sizeof(double)*nptsExpndB);}
  int stride = (ngrida/2+1);
  for(int js=0,j=offset;js<myNgridb;js+=num,j+=stride){
   for(int is=js,i=j;is<num+js;is++,i++){
     dataC[i] = msgData[is];
   }//endfor
  }//endfor

  delete msg;

//============================================================================
// Do the FFT when you have all the parts

  if(countIntGtoR[iopt]==rhoRsubplanes){
    countIntGtoR[iopt]=0;
    FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
    fftcache->doEextFFTGxToRx_Rchare(dataC,dataR,nplane_rho_x,ngrida,myNgridb);
    computeAtmForc(iopt);
  }//endif

//============================================================================
  }//end routine
//============================================================================



//============================================================================
// Get forces on atoms from Ewald or electrons
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::computeAtmForc(int flagEwd){
//============================================================================
#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in computeAtmForc w %d at %d %d %d\n",
             thisIndex.x,flagEwd,iterAtmTyp,nAtmTypRecv,natmTyp);
#endif
//============================================================================
// Check for errors

  // Ewald total SF  arrives on the last iteration only .
  // In parallel it can arrive before the last atmtyp SF
  if( (flagEwd==1) && (iterAtmTyp!=natmTyp)){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("You can't have flagEwd==1 unless you are on the last step\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  // If we are not on the last iteration when ewald can show up
  // the iteration # must match the number of SF received  
  if( (flagEwd==0) && (iterAtmTyp!=nAtmTypRecv) && (iterAtmTyp!=natmTyp)){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Atm SF Recv and atm SF calc out of sync %d\n",thisIndex.x);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  // the number of iteration <=natmTyp and you recieve natmtyp SF plus Ewald total
  if( iterAtmTyp>natmTyp || nAtmTypRecv>natmTyp+1 ){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Too much action for rhoRhart\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// set the variables to evaluate atom forces

  int myPlane      = thisIndex.x;
  int itime        = iteration;
  double *data;    
  if(flagEwd==0){data=atmSFR;}else{data=atmEwdSFR;}

  // you have already queried for this step:
  eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
  int *plane_index = eesData->RhoRHartData[myPlane].plane_index;
  int **igrid      = eesData->RhoRHartData[myPlane].igrid;
  double **dmn_x   = eesData->RhoRHartData[myPlane].dmn_x;
  double **dmn_y   = eesData->RhoRHartData[myPlane].dmn_y;
  double **dmn_z   = eesData->RhoRHartData[myPlane].dmn_z;

  AtomsGrp *ag         = atomsGrpProxy.ckLocalBranch(); // find me the local copy
  FastAtoms *fastAtoms = &(ag->fastAtoms);
  int natm             = ag->natm;

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

   if(rhoRsubplanes==1){
     CPLOCAL::eesAtmForceRchare(natm,fastAtoms,nAtmTypRecv,igrid,dmn_x,dmn_y,dmn_z,
                                plane_index,data,myPlane,flagEwd);
   }else{
    CkPrintf("rhart atmforc broken\n"); CkExit();
   }//endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(eesAtmForcR_, StartTime, CmiWallTimer());    
#endif

//============================================================================
// If you aren't done, start another round.

 if(iterAtmTyp<natmTyp){computeAtmSF();}

//============================================================================
// If your are done, Tell rho real and reset yourself for the next time step

#define _CP_DEBUG_STOP_RHART_OFF

#ifndef _CP_DEBUG_STOP_RHART_
 if(iterAtmTyp==natmTyp && nAtmTypRecv==(natmTyp+1)){
#ifdef _CP_RHART_VERBOSE_
   CkPrintf("HI, I am rhart chare %d in computeAtmForc w %d at %d done\n",
             thisIndex.x,flagEwd,iterAtmTyp);
#endif
   CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
   int index = (thisIndex.x % sim->sizeZ);
   rhoRealProxy(index,thisIndex.y).RHartReport();
   iterAtmTyp  = 0;
   nAtmTypRecv = 0;
 }//endif
#endif

//============================================================================
// If your are done and debugging, exit

#ifdef _CP_DEBUG_STOP_RHART_
 if(iterAtmTyp==natmTyp && nAtmTypRecv==(natmTyp+1)){
   CkPrintf("HI, I am rhart chare %d in computeAtmForc w %d at %d contrib\n",
              thisIndex.x,flagEwd,iterAtmTyp);
   int i=1;
   CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(NULL),rhoRHartExtProxy);
   contribute(sizeof(int),&i,CkReduction::sum_int,cb);
 }//endif
#endif

//============================================================================
  }// end routine
//============================================================================


//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Glenn's Rhart exit 
//============================================================================
void CP_Rho_RHartExt::exitForDebugging(){
  countDebug++;  
  if(countDebug==(rhoRsubplanes*ngridc)){
    countDebug=0;
    CkPrintf("I am in the exitfordebuging rhoreal puppy. Bye-bye\n");
    CkExit();
  }//endif
}
//============================================================================
