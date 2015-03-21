//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_RHartExt.C
 *
 *  At the start of the program, the constructor CP_Rho_RHartExt is
 *  called, FFTs set up and registration with the atoms cache is
 *  performed. This chare is only needed for the N log N EES Ext method.
 *
 *  CP_Rho_RHartExt is started by a message sent from a launch
 *  point chare array selected by the user as a configurable
 *  parameter.
 *
 *  When the launch command arrives, the chare pops the atom cache for
 *  an EES approximate to the 1st atom type structure factor in real space.
 *  This is partly ffted to atmsf(gx,gy,z) then sent back to
 *  GHartExt. GHartEext does its dance and sends back data that can be
 *  used to compute the atom forces from the Ext energy of this atom
 *  type. The fft is completed and forces on the atoms added. The 
 *
 */ 
//============================================================================

#include "charm++.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "main/AtomsCache.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "main/eesCache.h"
#include "cp_state_ctrl/CP_State_Plane.h"

#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern CkVec <CProxy_CP_Rho_RealSpacePlane> UrhoRealProxy;
extern CkVec <CProxy_AtomsCache>              UatomsCacheProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>       UrhoRHartExtProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>       UrhoGHartExtProxy;
extern CkVec <CProxy_eesCache>              UeesCacheProxy;
extern CkVec <CProxy_FFTcache>              UfftCacheProxy;

extern ComlibInstanceHandle commRHartGHartIns;

extern int    sizeX;
extern Config config;


//#define _CP_RHART_VERBOSE_

/** \file CP_Rho_RHartExt.C
 * @addtogroup Density
 @{
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *  This object just gets a rho message, computes GHartExt, and sends
 *  vks.  
 */
//============================================================================
CP_Rho_RHartExt::CP_Rho_RHartExt(int _ngrida, int _ngridb, int _ngridc, 
    int _ees_eext_on, int _natmTyp, 
    UberCollection _instance) 
: thisInstance(_instance)
{
  //============================================================================
  // Set some pareameters

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  nplane_rho_x         = sim->nplane_rho_x;
  listSubFlag          = sim->listSubFlag;
  rhoRsubplanes        = config.rhoRsubplanes;

  CkAssert(nplane_rho_x >= rhoRsubplanes); // safety : should already be checked.

  ngrida      = _ngrida; 
  ngridb      = _ngridb; 
  ngridc      = _ngridc; 
  ees_eext_on = _ees_eext_on;
  natmTypTot  = _natmTyp;

  nchareHartAtmT = config.nchareHartAtmT;
  int div        = (natmTypTot/nchareHartAtmT); 
  int rem        = (natmTypTot % nchareHartAtmT);
  int max        = (thisIndex.z < rem ? thisIndex.z : rem);
  natmTyp        = (thisIndex.z<rem ? div+1 : div);  
  atmTypoff      = div*thisIndex.z + max;         

  iplane_ind  = thisIndex.y*ngridc + thisIndex.x;

  //============================================================================
  // Compute messages sizes and zero message counters

  int nchareG = sim->nchareRhoGEext;
  if(rhoRsubplanes>1){
    recvCountFromGHartExt = 0;
    for(int i=0;i<nchareG;i++){
      if(sim->nline_send_eext_y[i][thisIndex.y]>0){recvCountFromGHartExt++;}
    }//endfor
  }else{
    recvCountFromGHartExt=nchareG;
  }//endif

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

  //============================================================================
  // Parallelization 

  //------------------------------------------------------
  // rho(x,y,z) by (y,z)
  div      = (ngridb/rhoRsubplanes); 
  rem      = (ngridb % rhoRsubplanes);
  max      = (thisIndex.y < rem ? thisIndex.y : rem);
  myNgridb     = (thisIndex.y<rem ? div+1 : div);  // number of y values/lines of x
  myBoff       = div*thisIndex.y + max;            // offset into y 
  nptsB        =  ngrida*myNgridb;                 // size of plane without extra room
  nptsExpndB   = (ngrida+2)*myNgridb;              // extra memory for RealToComplex FFT

  //------------------------------------------------------
  // Parallelization after transpose : rho(gx,y,z) : parallelize gx and z
  if(rhoRsubplanes>1){
    myNplane_rho = sim->numSubGx[thisIndex.y];
  }else{
    myNplane_rho = nplane_rho_x;
  }//endif
  nptsA        = 2*myNplane_rho*ngridb;            // memory size for fft in doubles
  nptsExpndA   = 2*myNplane_rho*ngridb;            // memory size for fft in doubles


  //---------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *  Register in the atom cache etc.
 */
//============================================================================
void CP_Rho_RHartExt::init(){
  //============================================================================
  // Malloc data sets : Register in the cache

  csize        = 0;
  csizeInt     = 0;
  if(ees_eext_on==1){
    csize        = (ngrida/2+1)*myNgridb;  // complex variable size

    eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    eesData->registerCacheRHart(thisIndex.x);

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
    //    CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(NULL),UrhoRHartExtProxy[thisInstance.proxyOffset]);
    CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(NULL),thisProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

  //============================================================================
  // Set up proxies

  rhoGHartProxy_com    = UrhoGHartExtProxy[thisInstance.proxyOffset];;
#ifdef USE_COMLIB
  if (config.useRHartInsGHart){
    ComlibAssociateProxy(commRHartGHartIns,rhoGHartProxy_com);          
  }//endif
#endif

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

  ArrayElement3D::pup(p);
  p|listSubFlag;
  p|countDebug;
  p|nplane_rho_x;
  p|rhoRsubplanes;
  p|registrationFlag;
  p|launchFlag;
  p|ngrida;
  p|ngridb;
  p|ngridc;
  p|iplane_ind;
  p|ees_eext_on;
  p|natmTypTot;
  p|atmTypoff;
  p|natmTyp;
  p|nchareHartAtmT;
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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Invoke by Rspace-density : Density has arrived in r-space and will soon arrive
 *                            in g-space. Get moving RhartExt
 */
//============================================================================
void CP_Rho_RHartExt::startEextIter(){
  //============================================================================
  // Check for error

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt = sim->cp_min_opt;

  iteration ++;
  // the atoms haven't moved yet
  // TODO: We'll need to make something like this for BOMD (keep iteration counters in sync)
  if(cp_min_opt==0){
    if(UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration != iteration-1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error in GHartExtVks : atoms slow %d %d\n",
          UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration,iteration);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

  //============================================================================
  // This is a new time step : Increment time step counter, zero atm typ counter
  //                           Launch if we are ready

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d %d in start : %d on %d\n",
      thisIndex.x,thisIndex.y,registrationFlag,CkMyPe());
#endif

  if(iterAtmTyp==natmTyp){CkPrintf("%d signing off\n",thisIndex.x);CkExit();}

  iterAtmTyp  = 0; 
  nAtmTypRecv = 0;

  launchFlag=1;
  if(registrationFlag==1){computeAtmSF();}

  //============================================================================
}//end routine
//============================================================================


//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Make sure everyone is registered on the 1st time step
 */
//==========================================================================
void CP_Rho_RHartExt::registrationDone(CkReductionMsg *msg) {
  //==========================================================================

  int sum = ((int *)msg->getData())[0];
  delete msg;

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d %d in reg : %d on %d\n",
      thisIndex.x,thisIndex.y,sum,CkMyPe());
#endif

  // use barrier as a debug tool
  if(iterAtmTyp==natmTyp){
    CkPrintf("HI, I am rhart chare %d exiting at %d\n",thisIndex.x,iterAtmTyp+1);
    CkExit();
  }//endif

  registrationFlag=1;
  if(launchFlag==1){computeAtmSF();}

}//end routine
//==========================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Start the real space part of the EES interpolation foratmSF(iatmTyp)
 */
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
  CkPrintf("HI, I am rhart chare %d %d %d in computeatmsf at iter %d on %d\n",
      thisIndex.x,thisIndex.y,thisIndex.z, iterAtmTyp,CkMyPe());
#endif

  //============================================================================
  // Look up some constants and fill the r-space grid

  int myPlane        = thisIndex.x;
  int mySubplane     = thisIndex.y;
  int itime          = iteration;
  int iterAtmTypFull = iterAtmTyp+atmTypoff;

  eesCache *eesData= UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  if(iterAtmTyp==1){eesData->queryCacheRHart(myPlane,itime,iterAtmTyp);}

  int *plane_index = eesData->RhoRHartData[myPlane]->plane_index;
  int **nSub       = eesData->RhoRHartData[myPlane]->nSub;
  int ***igrid     = eesData->RhoRHartData[myPlane]->igrid;
  double ***mn     = eesData->RhoRHartData[myPlane]->mn;

  AtomsCache *ag     = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch(); // find me the local copy
  int natm         = ag->natm;

  CPLOCAL::eesPackGridRchare(natm,iterAtmTypFull,atmSFR,myPlane,
      mySubplane,igrid,mn,
      plane_index,nSub,myNgridb);

#ifdef _CP_RHART_VERBOSE_
  int off = 0;
  for(int b = 0; b < myNgridb; b++) {
    for(int a = 0; a < 2*(ngrida/2+1); a++) {
      CkPrintf("%d %d %d %lf\n", myPlane, b, a, atmSFR[off]);
      off++;
    }
  }
#endif

  //============================================================================
  // FFT the result to G-space

  fftAtmSfRtoG();

  //============================================================================
}//end routine 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * FFT of SF(x,y,z,iatmTyp) -> SF(gx,gy,z,iatmTyp)
 **/
//============================================================================
void CP_Rho_RHartExt::fftAtmSfRtoG(){
  //============================================================================
#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d %d %d in fftsback at %d on %d\n",
      thisIndex.x,thisIndex.y,thisIndex.z, iterAtmTyp,CkMyPe());
#endif

  //==========================================================================
  // Do the FFT
#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif    

  FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  

  if(rhoRsubplanes>1){
    fftcache->doEextFFTRxToGx_Rchare(atmSFC,atmSFR,nplane_rho_x,ngrida,myNgridb,iplane_ind);
  }else{
    fftcache->doEextFFTRtoG_Rchare(atmSFC,atmSFR,nplane_rho_x,ngrida,ngridb,iplane_ind);
  }//endif

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(doEextFFTRtoG_, StartTime, CmiWallTimer());    
#endif

  //==========================================================================
  // Do the communication if you are not debugging

  if(rhoRsubplanes>1){
    sendAtmSfRyToGy();   // double transpose method (yz ---> gx,z)
  }else{
#ifdef DEBUG_INT_TRANS_FWD
    char name[100];
    sprintf(name,"partFFTGxGyZ%d.out.%d.%d",rhoRsubplanes,thisIndex.x,thisIndex.y);
    FILE *fp = fopen(name,"w");
    for(int ix =0;ix<nplane_rho_x;ix++){
      for(int iy =0;iy<ngridb;iy++){
        int i = iy*(ngrida/2+1) + ix;
        fprintf(fp,"%d %d : %g %g\n",iy,ix,atmSFC[i].re,atmSFC[i].im);
      }//endfor
    }//endof
    fclose(fp);
    UrhoRHartExtProxy[thisInstance.proxyOffset](0,0,0).exitForDebugging();
#else
#ifdef _GLENN_DBG_KPT_
    char name[100];
    sprintf(name,"sfAtmTypGxGyZ_inR.p%d.t%d.out",thisIndex.x,iterAtmTyp);
    FILE *fp = fopen(name,"w");
    for(int ix =0;ix<nplane_rho_x;ix++){
      for(int iy =0;iy<ngridb;iy++){
        int i = iy*(ngrida/2+1) + ix;
        fprintf(fp,"%d %d : %g %g\n",iy,ix,atmSFC[i].re,atmSFC[i].im);
      }//endfor
    }//endof
    fclose(fp);
#endif
    sendAtmSfRhoGHart(); // single transpose method (z ---> gx,gy)
#endif
  }//endif

  //============================================================================
}// end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Double Transpose Fwd Send : A(gx,y,z) on the way to A(gx,gy,z)
 *                             Send so that (y,z) parallelism is 
 *                             switched to (gx,z)
 **/
//============================================================================
void CP_Rho_RHartExt::sendAtmSfRyToGy(){
  //============================================================================
  // Launch the communication

  CkAssert(rhoRsubplanes>1);
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int **listSubGx      = sim->listSubGx;
  int  *numSubGx       = sim->numSubGx;

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
    int num   = numSubGx[ic]; // number of gx values chare ic wants
    int size  = num*myNgridb; // num*(all the y's i have)

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

    if(listSubFlag==1){
      for(int i=0,koff=0;i<num;i++,koff+=myNgridb){
        for(int k=koff,ii=listSubGx[ic][i];k<myNgridb+koff;k++,ii+=stride){
          data[k] = atmSFC[ii];  // all y's of this gx
        }//endfor
      }//endfor
    }else{
      int nst=listSubGx[ic][0];
      for(int i=0,ist=nst,koff=0;i<num;i++,koff+=myNgridb,ist++){
        for(int k=koff,ii=ist;k<myNgridb+koff;k++,ii+=stride){
          data[k] = atmSFC[ii];  // all y's of this gx
        }//endfor
      }//endfor
    }//endif

    UrhoRHartExtProxy[thisInstance.proxyOffset](ix,ic,thisIndex.z).recvAtmSfRyToGy(msg);

#ifdef _ERIC_SETS_UP_COMMLIB_
    UrhoRHartExtProxy[thisInstance.proxyOffset](ix,ic,thisIndex.z).recvAtmSfRyToGy(msg);
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
/**
 * Double Transpose Fwd Recv : A(gx,y,z) on the way to A(gx,gy,z)
 *                             Recv so that (y,z) parallel switched to (gx,z)
 *
 * Invoked natm_typ times per algorithm step : 
 */
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

  //============================================================================

  for(int js=0,j=offset;js<size;js+=num,j+=ngridb){
    for(int is=js,i=j;is<num+js;is++,i++){
      atmSFCint[i] = msgData[is];
    }//endfor
  }//endfor

  delete msg;

  //============================================================================
  // Do the Y fft, invoke communication if not debugging

  countIntRtoG++;
  if(countIntRtoG==rhoRsubplanes){

    countIntRtoG = 0;
    FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif
    fftcache->doEextFFTRyToGy_Rchare(atmSFCint,atmSFRint,myNplane_rho,
        ngrida,ngridb,iplane_ind);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(doEextFFTRytoGy_, StartTime, CmiWallTimer());    
#endif

#ifndef DEBUG_INT_TRANS_FWD
    sendAtmSfRhoGHart(); 
#endif
  }//endif

  //============================================================================
  // Debugging 

#ifdef DEBUG_INT_TRANS_FWD
  if(countIntRtoG==rhoRsubplanes){
    CPcharmParaInfo *sim =  CPcharmParaInfo::get();
    int **listSubGx = sim->listSubGx;
    int ic          = thisIndex.y;
    char name[100];
    sprintf(name,"partFFTGxGyZ%d.out.%d.%d",rhoRsubplanes,thisIndex.x,thisIndex.y);
    FILE *fp = fopen(name,"w");
    for(int ix =0;ix<myNplane_rho;ix++){
      for(int iy =0;iy<ngridb;iy++){
        int i = ix*ngridb + iy;
        fprintf(fp,"%d %d : %g %g\n",iy,listSubGx[ic][ix],atmSFCint[i].re,atmSFCint[i].im);
      }//endfor
    }//endof
    fclose(fp);
    UrhoRHartExtProxy[thisInstance.proxyOffset](0,0,0).exitForDebugging();
  }//endif
#endif

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Send EESatmSF(gx,gy,z,iatmTYP) to g-space whence the FFT will be
 * completed
 **/
//============================================================================
void CP_Rho_RHartExt::sendAtmSfRhoGHart(){
  //============================================================================

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d %d %d in sendsf at %d on %d\n",
      thisIndex.x,thisIndex.y,thisIndex.z,iterAtmTyp,CkMyPe());
#endif

  CPcharmParaInfo *sim      = CPcharmParaInfo::get();
  int nchareG               = sim->nchareRhoGEext;

  int **tranpack            = sim->index_tran_upack_eext;
  int *nlines_per_chareG    = sim->nlines_per_chareRhoGEext;

  int ***tranupack_Y        = sim->index_tran_upack_eext_y;
  int **nlines_per_chareG_Y = sim->nline_send_eext_y;

  //===================================================================
  // Perform the transpose and then the blast off the final 1D-FFT

  //----------------------------------------------------------------
  // start commlib
#ifdef USE_COMLIB
  if(rhoRsubplanes==1){
#ifdef OLD_COMMLIB
    if(config.useRHartInsGHart){commRHartGHartIns.beginIteration();}
#else
    if(config.useRHartInsGHart){ComlibBegin(rhoGHartProxy_com,0);}
#endif
  }//endif
#endif

  //----------------------------------------------------------------
  // do the send

  int iy = thisIndex.y;
  for (int ic=0;ic<nchareG;ic++) { // chare arrays to which we will send

    //---------------------
    // malloc the message
    int sendFFTDataSize = nlines_per_chareG[ic];
    if(rhoRsubplanes!=1){sendFFTDataSize = nlines_per_chareG_Y[ic][iy];}
    if(sendFFTDataSize>0)
    {
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
        rhoGHartProxy_com(ic,thisIndex.z).recvAtmSFFromRhoRHart(msg); // send the message
      }else{
        UrhoGHartExtProxy[thisInstance.proxyOffset](ic,thisIndex.z).recvAtmSFFromRhoRHart(msg); // send the message
      }//endif
    }
  }//end for : chare sending

  //----------------------------------------------------------------
  // stop commlib
#ifdef USE_COMLIB
  if(rhoRsubplanes==1){
#ifdef OLD_COMMLIB
    if(config.useRHartInsGHart){commRHartGHartIns.endIteration();}
#else
    if(config.useRHartInsGHart){ComlibBegin(rhoGHartProxy_com,1);}
#endif
  }//endif
#endif

  //============================================================================
}// end routine
//============================================================================


//============================================================================
// 
// 
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Hartree sends back atom forces from ext interation
 * Depending on the flag, it is Ewald or e-atm interation
 */
//============================================================================
void CP_Rho_RHartExt::recvAtmForcFromRhoGHart(RhoRHartMsg *msg){
  //============================================================================

  CPcharmParaInfo *sim      = CPcharmParaInfo::get();
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
  if (countFFT[iopt] > 	recvCountFromGHartExt) {
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
    CkPrintf("Dude.2, %d != %d for rho rhart chare %d %d %d\n",iter,iterAtmTyp,
        thisIndex.x,Index,iopt);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  if((iopt==1) && (iterAtmTyp!=natmTyp) ){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Dude.3, %d != %d for rho rhart chare %d %d : %d %d\n",iter,natmTyp,
        thisIndex.x,iy,Index,iopt);
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

  if (countFFT[iopt] == recvCountFromGHartExt){
#ifdef _CP_RHART_VERBOSE_
    CkPrintf("HI, I am rhart chare %d %d %d in recvsf with opt %d at %d on %d\n",
        thisIndex.x,thisIndex.y, thisIndex.z,iopt,iterAtmTyp,CkMyPe());
#endif
    countFFT[iopt]=0;
    if(rhoRsubplanes==1){nAtmTypRecv++;}
    fftAtmForcGtoR(iopt);
  }//endif

  //============================================================================
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
    Start the FFT : eesExtAtmForc(gx,gy,z) -> eesExtAtmForc(x,y,z)
 */
//============================================================================
void CP_Rho_RHartExt::fftAtmForcGtoR(int flagEwd){
  //============================================================================

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d %d %d in FFTAtmForc w %d at %d on %d\n",
      thisIndex.x,thisIndex.y,thisIndex.z,flagEwd,iterAtmTyp,CkMyPe());
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

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif    

  FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
  if(rhoRsubplanes==1){  
    fftcache->doEextFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb,iplane_ind);
    computeAtmForc(flagEwd);
  }else{
    fftcache->doEextFFTGyToRy_Rchare(dataC,dataR,myNplane_rho,ngrida,ngridb,iplane_ind);
    sendAtmForcGxToRx(flagEwd);
  }//endif

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(doEextFFTGtoR_, StartTime, CmiWallTimer());    
#endif

  //============================================================================ 
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *   Send out the partially transformed EES atom force for Rhosubplanes>1
 */
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
    msg->offset      = thisIndex.y;  // my chare index
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

    //HEY is thisIndex.z right here?
    UrhoRHartExtProxy[thisInstance.proxyOffset](ix,ic,thisIndex.z).recvAtmForcGxToRx(msg);

#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){
      case 0 : UrhoGProxy[thisInstance.proxyOffset]_com(ic,0).acceptRhoGradVksGxToRx(msg);break;
      case 1 : UrhoGProxy[thisInstance.proxyOffset]IGX_com(ic,0).acceptRhoGradVksGxToRx(msg); break;
      default: CkAbort("impossible iopt");break;
    }//end switch
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
/**
 *  Recieve the partially transformed EES atom structure factor for Rhosubplanes>1
 */
//============================================================================
void CP_Rho_RHartExt::recvAtmForcGxToRx(RhoGHartMsg *msg){
  //============================================================================

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int **listSubGx      = sim->listSubGx;
  int  *numSubGx       = sim->numSubGx;

  int size         = msg->size;  // msg size
  int iopt         = msg->iopt;  
  int iter         = msg->iter;  
  int num          = msg->num;   // number of lines along `a' sent
  int offset       = msg->offset; // offset into lines along `a'
  complex *msgData = msg->data;

  CkAssert(rhoRsubplanes>1);
  CkAssert(size==myNgridb*num);
  if(iterAtmTyp!=iter){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("iter!=iterAtmTyp recvatmforcgxtorx : %d %d %d : %d %d %d\n",
        thisIndex.x,thisIndex.y,thisIndex.z,iterAtmTyp,iter,iopt);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

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

  int stride = ngrida/2+1;
  if(listSubFlag==1){
    for(int js=0,j=0;js<size;js+=num,j++){
      int jj = j*stride;
      for(int is=js,i=0;is<(num+js);is++,i++){
        dataC[(listSubGx[offset][i]+jj)] = msgData[is];
      }//endfor
    }//endfor
  }else{
    int nst = listSubGx[offset][0];
    for(int js=0,j=0;js<size;js+=num,j++){
      int jj = j*stride+nst;
      for(int is=js,i=jj;is<(num+js);is++,i++){
        dataC[i] = msgData[is];
      }//endfor
    }//endfor
  }//endif
  delete msg;

  //============================================================================
  // Do the FFT when you have all the parts

  if(countIntGtoR[iopt]==rhoRsubplanes){
    countIntGtoR[iopt]=0;
    FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif
    fftcache->doEextFFTGxToRx_Rchare(dataC,dataR,nplane_rho_x,ngrida,myNgridb,iplane_ind);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(doEextFFTGxtoRx_, StartTime, CmiWallTimer());    
#endif

    nAtmTypRecv++;
    computeAtmForc(iopt);
  }//endif

  //============================================================================
}//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Get forces on atoms from Ewald or electrons by EES method
 */
//============================================================================
void CP_Rho_RHartExt::computeAtmForc(int flagEwd){
  //============================================================================
#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d %d %d in computeAtmForc.1 w %d at %d %d %d on %d\n",
      thisIndex.x,thisIndex.y, thisIndex.z, flagEwd,iterAtmTyp,nAtmTypRecv,natmTyp,CkMyPe());
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

  int myPlane        = thisIndex.x;
  int mySubplane     = thisIndex.y;
  int iterAtmTypFull = iterAtmTyp+atmTypoff;

  double *data;    
  if(flagEwd==0){data=atmSFR;}else{data=atmEwdSFR;}

  // you have already queried for this step:
  eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  AtomsCache *ag       = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch(); // find me the local copy

  CkAssert(eesData->allowedRhoRHartChares[myPlane]!=0);
  int *plane_index   = eesData->RhoRHartData[myPlane]->plane_index;
  int **nSub         = eesData->RhoRHartData[myPlane]->nSub;
  int ***igrid       = eesData->RhoRHartData[myPlane]->igrid;
  double ***dmn_x    = eesData->RhoRHartData[myPlane]->dmn_x;
  double ***dmn_y    = eesData->RhoRHartData[myPlane]->dmn_y;
  double ***dmn_z    = eesData->RhoRHartData[myPlane]->dmn_z;

  FastAtoms *fastAtoms = &(ag->fastAtoms);
  int natm             = ag->natm;

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif    
#ifdef _CP_RHART_VERBOSE_
  int off = 0;
  char myFileName[MAX_CHAR_ARRAY_LENGTH];
  sprintf(myFileName, "atmsf_%d_%d_%d.out", thisIndex.x, thisIndex.y, nAtmTypRecv);
  FILE *fp = fopen(myFileName,"w");
  for(int b = 0; b < myNgridb; b++) {
    for(int a = 0; a < 2*(ngrida/2+1); a++) {
      fprintf(fp, "%d %d %d %lf\n", myPlane, b, a, data[off]);
      off++;
    }
  }
  fclose(fp);
#endif

  CPLOCAL::eesAtmForceRchare(natm,fastAtoms,iterAtmTypFull,igrid,
      dmn_x,dmn_y,dmn_z, plane_index,nSub,data,
      myPlane,mySubplane,flagEwd);
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesAtmForcR_, StartTime, CmiWallTimer());    
#endif

  //============================================================================
  // If you aren't done, start another round.

  if(iterAtmTyp<natmTyp){computeAtmSF();}

  //============================================================================
  // If your are done, Tell rho real and reset yourself for the next time step

  int nrecvMax=natmTyp;
  if(thisIndex.z==0){nrecvMax++;}

#define _CP_DEBUG_STOP_RHART_OFF

#ifndef _CP_DEBUG_STOP_RHART_
  if(iterAtmTyp==natmTyp && nAtmTypRecv==nrecvMax){
#ifdef _CP_RHART_VERBOSE_
    CkPrintf("HI, I am rhart chare %d %d %d in computeAtmForc.2 w %d at %d done on %d\n",
        thisIndex.x,thisIndex.y,thisIndex.z,flagEwd,iterAtmTyp,CkMyPe());
#endif
    CPcharmParaInfo *sim   = CPcharmParaInfo::get();
    int index = (thisIndex.x % sim->sizeZ);
    UrhoRealProxy[thisInstance.proxyOffset](index,thisIndex.y).RHartReport();
    iterAtmTyp  = 0;
    nAtmTypRecv = 0;
  }//endif
#endif

  //============================================================================
  // If your are done and debugging, exit

#ifdef _CP_DEBUG_STOP_RHART_
  if(iterAtmTyp==natmTyp && nAtmTypRecv==nrecvMax){
    CkPrintf("HI, I am rhart chare %d in computeAtmForc w %d at %d contrib\n",
        thisIndex.x,flagEwd,iterAtmTyp);
    int i=1;
    CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(NULL),UrhoRHartExtProxy[thisInstance.proxyOffset]);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif
#endif

  //============================================================================
}// end routine
//============================================================================


//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Glenn's RhartExt exit  - handy, dandy debugging exit call
 */
//============================================================================
void CP_Rho_RHartExt::exitForDebugging(){
  countDebug++;  
  if(countDebug==(rhoRsubplanes*ngridc*nchareHartAtmT)){
    countDebug=0;
    CkPrintf("I am in the exitfordebuging rhorhartext puppy. Bye-bye\n");
    CkExit();
  }//endif
}
//============================================================================
/*@}*/
