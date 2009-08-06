/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_RealSpacePlane.C
 * This is the description of the "life" of a CP_Rho_RealSpacePlane object.
 *
 * At the start of the program, the constructor CP_Rho_RealSpacePlane() is called.
 * 
 * The CP_State_RealSpacePlanes send the squared magnitudes of the psi_r values 
 * using the acceptData() method. The squared magnitudes are summed across states.
 * A copy of this data is made, inverse fft is done in the z and x directions 
 * and sent to rhoGDensity. The other copy is processed using CP_exc_calc. 
 * Then the CP_Rho_RealSpacePlane object waits for a reply from the RhoGDensity 
 * object.
 *
 * The reply from RhoGDensity comes in the form of the method 
 * acceptDensityForSumming(). The data obtained from this reply is taken and
 * forward fft in z and x directions is performed. The resulting data is 
 * summed with the result of CP_exc_calc. The sum is sent to the 
 * CP_State_RealSpacePlane objects.
 */
//============================================================================

#include "charm++.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "main/groups.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_ctrl/CP_State_Plane.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================
extern CProxy_TimeKeeper                 TimeKeeperProxy;
extern CkVec <CProxy_CP_State_RealSpacePlane>    UrealSpacePlaneProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
extern CkVec <CProxy_CP_Rho_GSpacePlane>         UrhoGProxy;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>       UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;

extern ComlibInstanceHandle commRealInstance;
extern ComlibInstanceHandle commRealIGXInstance;
extern ComlibInstanceHandle commRealIGYInstance;
extern ComlibInstanceHandle commRealIGZInstance;
extern ComlibInstanceHandle mcastInstance;
extern CkGroupID            mCastGrpId;

extern Config    config;
extern int       nstates;

bool is_pow2(int );


//#define  _CP_RHO_RSP_VERBOSE_OFF

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the states
// Performs lots of operations to get exc, eext energies and vks
//
//============================================================================
CP_Rho_RealSpacePlane::CP_Rho_RealSpacePlane(int xdim, bool _useCommlib, 
                                             int _ees_eext_on,int _ngridcEext,
					     int _rhokeeperid, 
					     UberCollection _instance) :
  thisInstance(_instance)
//============================================================================
   {//begin routine
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
    CkPrintf("[%d %d] RhoR constructs \n",thisIndex.x, thisIndex.y);
#endif

//============================================================================
// Get parameters from the globals/groups

    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;

    ngridcEext           = _ngridcEext;    // FFT sizes for rharteext
    ngrida               = sim->sizeX;     // FFT sizes for rho
    ngridb               = sim->sizeY;
    ngridc               = sim->sizeZ;
    nplane_rho_x         = sim->nplane_rho_x; // gx_max 
    listSubFlag          = sim->listSubFlag;

    iplane_ind           = thisIndex.y*ngridc + thisIndex.x;

    rhoRsubplanes        = config.rhoRsubplanes;
    CkAssert(nplane_rho_x >= rhoRsubplanes); // safety : should already be checked.

    cp_grad_corr_on      = sim->cp_grad_corr_on;
    ees_eext_on          = _ees_eext_on;
    rhoKeeperId          = _rhokeeperid;

    double vol           = sim->vol;
    int numFFT           = config.numFFTPoints;
    FFTscale     = 1.0/((double)numFFT);  // these are based on the full size
    volumeFactor = vol*FFTscale;
    probScale    = 1.0/vol;

//============================================================================
// Compute number of messages to be received 

    int nchareRhoG=sim->nchareRhoG;
    if(rhoRsubplanes>1){
	recvCountFromGRho = 0;
	for(int i=0;i<nchareRhoG;i++){
	  if(sim->nline_send_rho_y[i][thisIndex.y]>0){recvCountFromGRho++;}
        }//endfor
    }else{
	recvCountFromGRho=nchareRhoG;
    }//endif

    int nchareRhoGEext=sim->nchareRhoGEext;
    if(rhoRsubplanes>1){
       recvCountFromGHartExt = 0;
       for(int i=0;i<nchareRhoGEext;i++){
         if(sim->nline_send_eext_y[i][thisIndex.y]>0)recvCountFromGHartExt++;
       }//endfor
    }else{
      recvCountFromGHartExt=nchareRhoGEext;
    }//endif

//============================================================================
// Parallelization 

   // Before transpose : rho(x,y,z) : parallelize y and z
    int div      = (ngridb/rhoRsubplanes); 
    int rem      = (ngridb % rhoRsubplanes);
    int max      = (thisIndex.y < rem ? thisIndex.y : rem);
    myNgridb     = (thisIndex.y<rem ? div+1 : div);  // number of y values/lines of x
    myBoff       = div*thisIndex.y + max;            // offset into y 
    nptsB        =  ngrida*myNgridb;                 // size of plane without extra room
    nptsExpndB   = (ngrida+2)*myNgridb;              // extra memory for RealToComplex FFT

   // After transpose : rho(gx,y,z) : parallelize gx and z
    if(rhoRsubplanes>1){
      myNplane_rho = sim->numSubGx[thisIndex.y];
    }else{
      myNplane_rho = nplane_rho_x;
    }//endif
    nptsA        = 2*myNplane_rho*ngridb;            // memory size fft in doubles
    nptsExpndA   = 2*myNplane_rho*ngridb;            // memory size fft in doubles

//============================================================================    
// Set up the data class : mallocs rho,gradrho, etc.

    initRhoRealSlab(&rho_rs,ngrida,myNgridb,ngridc,myNplane_rho,ngridb,
                    thisIndex.x,thisIndex.y,rhoRsubplanes);

//============================================================================
// Initialize counters, set booleans.myTime

    myTime          = 0;
    countDebug      = 0; // does nothing in the working code.
    countWhiteByrd  = 0;
    doneGradRhoVks  = 0;
    countRHart      = 0;
    countFFTRyToGy  = 0;

    countRHartValue = 1; if(thisIndex.x<(ngridcEext-rho_rs.sizeZ)){countRHartValue=2;}
    countRHartValue*=(config.nchareHartAtmT);

    doneHartVks     = false;
    doneRHart       = false;
    doneWhiteByrd   = false;
    for(int i=0;i<5;i++){countGradVks[i]=0;}//ctrls bkc-transpose rho(gx,gy,z):gxy->gx/z
    for(int i=0;i<5;i++){countIntGtoR[i]=0;}//ctrls bkc-int-transpose rho(gx,y,z):gx/z->yz
    for(int i=0;i<5;i++){countIntRtoG[i]=0;}//ctrls fwd-int-transpose rho(gx,y,z):yz->gx/z


//============================================================================
// Migration

    usesAtSync = CmiTrue;
   //    if(config.lbdensity){
   //      setMigratable(true);
   //    }else{
    setMigratable(false);
   //    }//endif

//============================================================================
   }//end routine
//============================================================================

/** post constructor initialization */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::init(){

// make sections in the realSpacePlane array. These will be used when 
// computing real-space densities and multicasting v_ks values 

    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    CkArrayIndexMax *elems   = new CkArrayIndexMax[nstates];

    int j;
    // section i has al the portions with all 
    CkArrayIndex2D idx(0, thisIndex.x);
    if(is_pow2(nstates)){
	for (j = 0; j < nstates; j++) {
	  idx.index[0] = j^((thisIndex.x+thisIndex.y)%nstates);
	    elems[j] = idx;
	}//endfor
    }else{
	for (j = 0; j < nstates; j++) {
	    idx.index[0] = (j+thisIndex.x+thisIndex.y)%nstates;
	    elems[j] = idx;
	}//endfor
    }//endif

    realSpaceSectionProxy = CProxySection_CP_State_RealSpacePlane::
        ckNew(UrealSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(), elems, nstates);

    realSpaceSectionCProxy = CProxySection_CP_State_RealSpacePlane::
        ckNew(UrealSpacePlaneProxy[thisInstance.proxyOffset].ckGetArrayID(), elems, nstates);

    realSpaceSectionProxy.ckDelegate
        (CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());

    mcastGrp->setSection(realSpaceSectionProxy);

#ifdef USE_COMLIB
    ComlibAssociateProxy(&mcastInstance,realSpaceSectionCProxy);
#endif

    delete [] elems;   

    ProductMsg *dummyProductMessage = new (0) ProductMsg;    
    // inform realspace element of this section proxy.
    dummyProductMessage->subplane=thisIndex.y;
    realSpaceSectionProxy.init(dummyProductMessage);

    rhoGProxy_com    = UrhoGProxy[thisInstance.proxyOffset];
    rhoGProxyIGX_com = UrhoGProxy[thisInstance.proxyOffset];
    rhoGProxyIGY_com = UrhoGProxy[thisInstance.proxyOffset];
    rhoGProxyIGZ_com = UrhoGProxy[thisInstance.proxyOffset];

#ifdef USE_COMLIB
    if (config.useRInsRhoGP) 
	ComlibAssociateProxy(&commRealInstance,rhoGProxy_com);          
    if (config.useRInsIGXRhoGP) 
	ComlibAssociateProxy(&commRealIGXInstance,rhoGProxyIGX_com);          
    if (config.useRInsIGYRhoGP) 
	ComlibAssociateProxy(&commRealIGYInstance,rhoGProxyIGY_com);          
    if (config.useRInsIGZRhoGP) 
	ComlibAssociateProxy(&commRealIGZInstance,rhoGProxyIGZ_com);          
#endif

}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_RealSpacePlane::~CP_Rho_RealSpacePlane(){
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//  Pup my variables for migration
//============================================================================
void CP_Rho_RealSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);
  p|listSubFlag;
  p|myTime;
  p|nplane_rho_x;
  p|ngrida;
  p|myNgridb;
  p|myNplane_rho;
  p|ngridb;
  p|ngridc;
  p|iplane_ind;
  p|nptsExpndA;
  p|nptsExpndB;
  p|nptsA;
  p|nptsB;
  p|ees_eext_on;
  p|ngridcEext;
  p|cp_grad_corr_on;
  p|FFTscale;        
  p|volumeFactor;        
  p|probScale;             
  p|rhoGHelpers;
  p|doneGradRhoVks;
  p|countWhiteByrd;
  p|doneWhiteByrd;
  p|doneHartVks;
  p|countDebug;
  p|rhoRsubplanes;

  PUParray(p,countGradVks,5);
  PUParray(p,countIntGtoR,5);
  PUParray(p,countIntGtoR,5);

  rho_rs.pup(p);   // pup your data class

 // Pupping Proxies???
  p|realSpaceSectionProxy;
  p|realSpaceSectionCProxy;
  if(p.isUnpacking()){
    rhoGProxy_com = UrhoGProxy[thisInstance.proxyOffset];
    rhoGProxyIGX_com = UrhoGProxy[thisInstance.proxyOffset];
    rhoGProxyIGY_com = UrhoGProxy[thisInstance.proxyOffset];
    rhoGProxyIGZ_com = UrhoGProxy[thisInstance.proxyOffset];

#ifdef USE_COMLIB
    if (config.useRInsRhoGP) 
	ComlibAssociateProxy(&commRealInstance,rhoGProxy_com);          
    if (config.useRInsIGXRhoGP) 
	ComlibAssociateProxy(&commRealIGXInstance,rhoGProxyIGX_com);          
    if (config.useRInsIGYRhoGP) 
	ComlibAssociateProxy(&commRealIGYInstance,rhoGProxyIGY_com);          
    if (config.useRInsIGZRhoGP) 
	ComlibAssociateProxy(&commRealIGZInstance,rhoGProxyIGZ_com);          
#endif

    }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Data comes from StateRspacePlane once an algorithm step.
/** 
 * Here the density from all the states is added up. The data from all the
 * states is received via an array section reduction. Nothing happens in this
 * chare until the density arrives.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity(CkReductionMsg *msg) {
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
    CkPrintf("RhoReal accepting Density %d %d %d\n",
              thisIndex.x,thisIndex.y,CkMyPe());
#endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->getSize()/sizeof(double) ;i++){
      CkAssert(isnan(((double*) msg->getData())[i])==0);
  }//endif
#endif
#ifdef _CP_SUBSTEP_TIMING_
  if(rhoKeeperId>0){
      double rhostart=CmiWallTimer();
      CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL),0,TimeKeeperProxy);
      contribute(sizeof(double),&rhostart,CkReduction::min_double, cb ,rhoKeeperId);
  }//endif
#endif

//============================================================================
// Set the flags : you are not done unless certain conditions apply.

   myTime++;
   doneHartVks   = false;
   doneWhiteByrd = false;
   doneRHart     = false;
   if(cp_grad_corr_on==0){doneWhiteByrd = true;}
   if(ees_eext_on==0)    {doneRHart     = true;}

#ifdef _CP_DEBUG_HARTEEXT_OFF_
   doneHartVks    = true;
   doneRHart      = true;
#endif

//============================================================================
// Unpack into spread out form and delete the message

   double *realValues = (double *) msg->getData(); 
   double *density    = rho_rs.density;
   CkAssert(msg->getSize() == nptsB * sizeof(double));

   rho_rs.uPackScaleGrow(density,realValues,probScale);

   delete msg;

//============================================================================
// If debugging, generate output!

#ifdef _CP_DEBUG_RHOR_RHO_
      char myFileName[MAX_CHAR_ARRAY_LENGTH];
      sprintf(myFileName, "Rho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
      FILE *fp = fopen(myFileName,"w");
        for (int i = 0; i <nptsB; i++){
          fprintf(fp,"%g\n",realValues[i]*probScale);
        }//endfor
      fclose(fp);
#endif

//============================================================================
// Compute the exchange correlation energy (density no-grad part)

  energyComputation();

//============================================================================
// 2nd Launch real-space external-hartree and the G-space non-local
// The energy comp through RhoG is he more expensive critical path.
   launchEextRNlG();


//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// The density is here : Launch ees NL and ees Ext routines
//
// Do this once an algorithm step
//
//============================================================================
void CP_Rho_RealSpacePlane::launchEextRNlG() {
//============================================================================
// Launch the external energy computation in r-space :
//   rhart has the same rhoRsubplanes for simplicity.

#ifndef _CP_DEBUG_HARTEEXT_OFF_
  if(ees_eext_on==1){
    int div    = (ngridcEext / ngridc);
    int rem    = (ngridcEext % ngridc);
    int ind    = thisIndex.x+ngridc;
    if(div>1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Eext Grid size too large for launch Scheme\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
#ifdef _CP_RHO_RSP_VERBOSE_
    CkPrintf("HI, I am r-rho chare %d lauchning %d : ngrids %d %d : %d\n",
             thisIndex.x,thisIndex.x,ngridcEext,ngridc,rem);
#endif
    for(int j=0;j<config.nchareHartAtmT;j++){
      UrhoRHartExtProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y,j).startEextIter();
      if(thisIndex.x<rem){
#ifdef _CP_RHO_RSP_VERBOSE_
	CkPrintf("HI, I am r-rho chare %d also lauchning %d\n",thisIndex.x,ind);
#endif
	UrhoRHartExtProxy[thisInstance.proxyOffset](ind,thisIndex.y,j).startEextIter();
      }//endif : the launch
    }//endfor : atmTyp parallism
  }//endif : Launch is needed
#endif

//============================================================================
// Launch nonlocal g space if it wasn't done in RS
//  Spread the launch over all the rhoRchares you can.

  CPcharmParaInfo *sim  = (scProxy.ckLocalBranch ())->cpcharmParaInfo;

  if(sim->ees_nloc_on==1 && config.launchNLeesFromRho==1 && sim->natm_nl>0){ 

      CkAssert(rho_rs.sizeZ>=config.nchareG);
      if(thisIndex.x<config.nchareG){
         int nstates = config.nstates; 
         int div     = (nstates/rhoRsubplanes);
         int rem     = (nstates % rhoRsubplanes);
         int add     = (thisIndex.y < rem ? 1 : 0);
         int max     = (thisIndex.y < rem ? thisIndex.y : rem);
         int ist     = div*thisIndex.y + max;
         int iend    = ist + div + add;
          for(int ns=ist;ns<iend;ns++){
           //           CkPrintf("RhoRP[%d,%d] triggering NL %d %d \n",
           //                    thisIndex.x, thisIndex.y, ns, thisIndex.x);
           CkAssert(ns<config.nstates);
           //           CkAssert(thisIndex.x<32);
           UgSpaceDriverProxy[thisInstance.proxyOffset](ns,thisIndex.x).startNonLocalEes(myTime);
         }//endfor
       }//endif

  }//endif : launch the non-local ees

//----------------------------------------------------------------------------
   }//end routine 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * Compute one part of the EXC energy using PINY CP_exc_calc.
 * This is done once an algorithm step
 */
//============================================================================
void CP_Rho_RealSpacePlane::energyComputation(){
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
   CkPrintf("In RhoRealSpacePlane[%d] energyComp %d %d\n",thisIndex.x, 
            CkMyPe(),CmiMemoryUsage());
#endif

//============================================================================

   double *density  = rho_rs.density;
   double *Vks      = rho_rs.Vks;
   int nf1          = ngrida;
   int nf2          = myNgridb;
   int nf3          = ngridc;
   int npts         = config.numFFTPoints;  // total number of points
   double *exc_ret  = &(rho_rs.exc_ret);
   double *muxc_ret = &(rho_rs.muxc_ret);

//============================================================================
// Perform exchange correlation computation (no grad corr here).

   CPXCFNCTS::CP_exc_calc(npts,nf1,nf2,nf3,density,Vks,exc_ret,muxc_ret,config.nfreq_xcfnctl);

#ifdef CMK_BLUEGENEL
   CmiNetworkProgress();
#endif

   if(cp_grad_corr_on==0){
     double exc[2];
     exc[0]=rho_rs.exc_ret;
     exc[1]=0.0;
     contribute(2*sizeof(double),exc,CkReduction::sum_double);
   }//endif

//============================================================================
// Invoke FFT to take rho(r) to rho(g) : do not over-write the density!!

   fftRhoRtoRhoG();

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// 1) Perform FFT of density: Single Transpose method rho(x,y,z) ---> rho(gx,gy,z)
//                            Double Transpose method rho(x,y,z) ---> rho(gx,y,z)
//
// 2) launch the real space part of the non-local 
//
//    Routine invoked once an algorithm step.
//
//============================================================================
void CP_Rho_RealSpacePlane::fftRhoRtoRhoG(){
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("In RhoRealSpacePlane[%d %d] FFT_RSpacetoGSpace %d %d\n",thisIndex.x, 
                   thisIndex.y, CkMyPe(),CmiMemoryUsage());
#endif

//============================================================================
// FFT myself back to G-space part way

  FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
  double  *density   = rho_rs.density;  // we need to save the density and vks.
  double  *dataR     = rho_rs.rhoIRX;   // rhoirx is around doing nothing now
  complex *dataC     = rho_rs.rhoIRXC;  // so we can use it to store the FFT

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

  rho_rs.uPackScale(dataR,density,volumeFactor);  // can't avoid the CmiMemcpy-scaling
  if(rhoRsubplanes>1){
    fftcache->doRhoFFTRxToGx_Rchare(dataC,dataR,nplane_rho_x,ngrida,myNgridb,iplane_ind);
  }else{
    fftcache->doRhoFFTRtoG_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb,iplane_ind);
  }//endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(RhoRtoGFFT_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_BLUEGENEL
  CmiNetworkProgress();
#endif


//============================================================================
// Launch the transpose :go directly to rhog if the 1 transpose method is on.
//                     : do an internal transpose if the 2 transpose method is on.

//#define DEBUG_INT_TRANS_FWD
  int iopt = 0;
  if(rhoRsubplanes>1){
    sendPartlyFFTRyToGy(iopt);    // double transpose method (yz ---> gx,z)
  }else{ 
#ifndef DEBUG_INT_TRANS_FWD
    sendPartlyFFTtoRhoG(iopt);    // single transpose method (z ---> gx,gy)
#else
    char name[100];
    sprintf(name,"partFFTGxGyZ%d.out.%d.%d",rhoRsubplanes,thisIndex.x,thisIndex.y);
    FILE *fp = fopen(name,"w");
    for(int ix =0;ix<nplane_rho_x;ix++){
      for(int iy =0;iy<ngridb;iy++){
        int i = iy*(ngrida/2+1) + ix;
        fprintf(fp,"%d %d : %g %g\n",iy,ix,dataC[i].re,dataC[i].im);
      }//endfor
    }//endof
    fclose(fp);
    UrhoRealProxy[thisInstance.proxyOffset](0,0).exitForDebugging();
#endif
  }//endif

//============================================================================
// Launch non-local real space FFT : allow NL to advance after density works a bit
// Now that we have previously done our bit for the critical path through rhog

  launchNLRealFFT();

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 *  Launch ees-nonlocal real here. This routine can be called from any
 *  other routine except gradCorr(). GradCorr comp is optional (PINY option)
 *  e.g. it is not always computed.
 */
//============================================================================
void CP_Rho_RealSpacePlane::launchNLRealFFT(){
//============================================================================
// Tell NLeesR its ok to compute its FFT. Otherwise we get no overlap
//   Spread the launch over all the RhoR chares

  CPcharmParaInfo *sim  = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  if(sim->ees_nloc_on==1){ 
    CkAssert(rho_rs.sizeZ>=sim->ngrid_nloc_c);;
    if(thisIndex.x < sim->ngrid_nloc_c){
        int nstates =  config.nstates; 
        int div     = (nstates/rhoRsubplanes);
        int rem     = (nstates % rhoRsubplanes);
        int add     = (thisIndex.y < rem ? 1 : 0);
        int max     = (thisIndex.y < rem ? thisIndex.y : rem);
        int ist     = div*thisIndex.y + max;
        int iend    = ist + div + add;
        for(int ns=ist;ns<iend;ns++){
          CkAssert(ns<config.nstates);
          UrealParticlePlaneProxy[thisInstance.proxyOffset](ns,thisIndex.x).launchFFTControl(myTime);
	  if(ns%4)
	    CmiNetworkProgress();
        }//endfor
    }//endif
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// Double Transpose Fwd Send : A(gx,y,z) on the way to A(gx,gy,z)
//                             Send so that (y,z) parallelism is 
//                             switched to (gx,z)
//
//  Invoked 4 times per algorithm step : case 0   density(gx,y,z)
//                                     : cast 1-3 gradients(gx,y,z)
//
//============================================================================
void CP_Rho_RealSpacePlane::sendPartlyFFTRyToGy(int iopt){
//============================================================================

    CPcharmParaInfo *sim  = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    int **listSubGx       = sim->listSubGx;
    int  *numSubGx        = sim->numSubGx;

    CkAssert(rhoRsubplanes>1);
    complex *FFTresult;
    switch(iopt){   
      case 0 : FFTresult = rho_rs.rhoIRXC;  break; //Scr variable for rho send
      case 1 : FFTresult = rho_rs.rhoIRXC;  break;
      case 2 : FFTresult = rho_rs.rhoIRYC;  break;
      case 3 : FFTresult = rho_rs.rhoIRZC;  break;
      default: CkAbort("impossible iopt");  break;
    }//end switch

//============================================================================
// Launch the communication

  //-----------------------------------------------------------------------------
  // Commlib launch : 

#ifdef USE_COMLIB
#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){
       case 0 : if(config.useRInsRhoRP)    commRealInstanceRx.beginIteration();    break;
       case 1 : if(config.useRInsIGXRhoRP) commRealIGXInstanceRx.beginIteration(); break;
       case 2 : if(config.useRInsIGYRhoRP) commRealIGYInstanceRx.beginIteration(); break;
       case 3 : if(config.useRInsIGZRhoRP) commRealIGZInstanceRx.beginIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch
#endif
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
      RhoGSFFTMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) RhoGSFFTMsg; 
      msg->size        = size;
      msg->iopt        = iopt;
      msg->offset      = myBoff;      // where the myNgridB y-lines start.
      msg->num         = myNgridb;    // number of y-lines I have. 
      complex *data    = msg->data;   // data

      if(config.prioFFTMsg){
          CkSetQueueing(msg, CK_QUEUEING_IFIFO);
          *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.y;
      }//endif

      if(listSubFlag==1){
        for(int i=0,koff=0;i<num;i++,koff+=myNgridb){
          for(int k=koff,ii=listSubGx[ic][i];k<myNgridb+koff;k++,ii+=stride){
            data[k] = FFTresult[ii];  // all y's of this gx
          }//endfor
        }//endfor
      }else{
        int nst=listSubGx[ic][0];
        for(int i=0,ist=nst,koff=0;i<num;i++,koff+=myNgridb,ist++){
          for(int k=koff,ii=ist;k<myNgridb+koff;k++,ii+=stride){
            data[k] = FFTresult[ii];  // all y's of this gx
          }//endfor
        }//endfor
      }//endif

      switch(iopt){
        case 0 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksRyToGy(msg);      break;
        case 1 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksRyToGy(msg);      break;
        case 2 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksRyToGy(msg);      break;
        case 3 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksRyToGy(msg);      break;
        default: CkAbort("impossible iopt");break;
      }//end switch

#ifdef _ERIC_SETS_UP_COMMLIB_
      switch(iopt){
        case 0 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksRyToGy(msg);      break;
        default: CkAbort("impossible iopt");break;
      }//end switch
#endif

#ifdef CMK_BLUEGENEL
      CmiNetworkProgress();
#endif
    }//end for : chare sending

  //-----------------------------------------------------------------------------
  // Commlib stop

#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){
       case 0 : if(config.useRInsRhoRP)    commRealInstanceRx.endIteration();    break;
       case 1 : if(config.useRInsIGXRhoRP) commRealIGXInstanceRx.endIteration(); break;
       case 2 : if(config.useRInsIGYRhoRP) commRealIGYInstanceRx.endIteration(); break;
       case 3 : if(config.useRInsIGZRhoRP) commRealIGZInstanceRx.endIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch
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
// Invoked 4 times per algorithm step : case 0   rho(gx,y,z)
//                                    : cast 1-3 gradRho(gx,y,z)
//
//============================================================================
void CP_Rho_RealSpacePlane::acceptRhoGradVksRyToGy(RhoGSFFTMsg *msg){
//============================================================================

  int size         = msg->size;  // msg size
  int iopt         = msg->iopt;  //
  int num          = msg->num;
  int offset       = msg->offset;
  complex *msgData = msg->data;

  CkAssert(size==myNplane_rho*num);
  CkAssert(rhoRsubplanes>1);

  complex *dataC;
  switch(iopt){
    case  0: dataC   = rho_rs.rhoIRXCint;  break;
    case  1: dataC   = rho_rs.rhoIRXCint;  break;
    case  2: dataC   = rho_rs.rhoIRYCint;  break;
    case  3: dataC   = rho_rs.rhoIRZCint;  break;
    default: CkAbort("Impossible option\n"); break;
  }//endif

  for(int js=0,j=offset;js<size;js+=num,j+=ngridb){
   for(int is=js,i=j;is<num+js;is++,i++){
     dataC[i] = msgData[is];
   }//endfor
  }//endfor

  delete msg;

  countIntRtoG[iopt]++;
  if(countIntRtoG[iopt]==rhoRsubplanes){
    countIntRtoG[iopt]=0;
    fftRhoRyToGy(iopt);
  }//endfor

//============================================================================
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Double Transpose Fwd FFT : A(gx,y,z) -> A(gx,gy,z)
//
//  Invoked 4 times per algorithm step : 
//               case 0   rho(gx,y,z)    -> rho(gx,gy,z)
//               cast 1-3 gradRho(gx,y,z)-> gradRho(gx,gy,z)
//
//============================================================================
void CP_Rho_RealSpacePlane::fftRhoRyToGy(int iopt){
//============================================================================
// FFT myself back to G-space part way : e.g. along gy here

  CkAssert(rhoRsubplanes>1);
  complex *dataC;
  double *dataR;
  switch(iopt){
    case 0: dataC = rho_rs.rhoIRXCint; dataR = rho_rs.rhoIRXint; break;
    case 1: dataC = rho_rs.rhoIRXCint; dataR = rho_rs.rhoIRXint; break;
    case 2: dataC = rho_rs.rhoIRYCint; dataR = rho_rs.rhoIRYint; break;
    case 3: dataC = rho_rs.rhoIRZCint; dataR = rho_rs.rhoIRZint; break;
    default: CkAbort("Impossible option\n"); break;
  }//endif

  FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
  fftcache->doRhoFFTRyToGy_Rchare(dataC,dataR,myNplane_rho,ngrida,ngridb,iplane_ind);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(doRhoFFTRytoGy_, StartTime, CmiWallTimer());    
#endif

//============================================================================
// Send chunk to RhoGDensity 

  int igo=0;
  if(iopt>=1 && iopt <= 3){countFFTRyToGy++; igo=1;}

#ifndef DEBUG_INT_TRANS_FWD
  if(config.rhoGToRhoRMsgComb==0 || iopt==0){sendPartlyFFTtoRhoG(iopt);}
  if(config.rhoGToRhoRMsgComb==1 && countFFTRyToGy==3 && igo==1){
    countFFTRyToGy=0;
    sendPartlyFFTtoRhoGall();
  }//endif
#else
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int **listSubGx = sim->listSubGx;
  int ic          = thisIndex.y;
  CkPrintf("%d %d : %d\n",thisIndex.x,thisIndex.y,myNplane_rho);
  char name[100];
  sprintf(name,"partFFTGxGyZ%d.out.%d.%d",rhoRsubplanes,thisIndex.x,thisIndex.y);
  FILE *fp = fopen(name,"w");
  for(int ix =0;ix<myNplane_rho;ix++){
   for(int iy =0;iy<ngridb;iy++){
     int i = ix*ngridb + iy;
     fprintf(fp,"%d %d : %g %g\n",iy,listSubGx[ic][ix],dataC[i].re,dataC[i].im);
   }//endfor
  }//endof
  fclose(fp);
  UrhoRealProxy[thisInstance.proxyOffset](0,0).exitForDebugging();
#endif

//============================================================================
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// The Tranpose to G-space : A(gx,gy,z) on the way to A(gx,gy,gz)
//                           Change parallel by gx,z to parallel by {gx,gy}
//                           We switch chare arrays here from RhoR to RhoG
//
//  Invoked 4 times per algorithm step : 
//               case 0   send rho(gx,gy,z)     -> rho(gx,gy,z) in Rhog
//               cast 1-3 send gradRho(gx,gy,z) -> gradRho(gx,gy,z) in Rhog
//
// Send the partly FFTed array A(gx,gy,z), to rhoGSpacePlane :  
// You then wait while RHOG works, sends back RhoIRalpha which you whitebyrdize
// You send back the whitebyrdized RhoIRalpha puppies to RhoGspacePlane.
// Whilst all this is going on, HartG is churning and will send you another
// part of Vks.  All this requires keeping density, vks, div_rho_alpha and
// vkshart memory available at all times to receive messages.
//
//============================================================================
void CP_Rho_RealSpacePlane::sendPartlyFFTtoRhoG(int iopt){
//============================================================================
// Local pointers and variables

    CPcharmParaInfo *sim        = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    int nchareRhoG              = sim->nchareRhoG;
    int **tranpack_rho          = sim->index_tran_upack_rho;
    int *nlines_per_chareRhoG   = sim->nlines_per_chareRhoG;
    int ***tranupack_rhoY       = sim->index_tran_upack_rho_y;
    int **nlines_per_chareRhoGY = sim->nline_send_rho_y;

    complex *FFTresult;
    if(rhoRsubplanes==1){
     switch(iopt){   
      case 0 : FFTresult = rho_rs.rhoIRXC;  break; //Scr variable for rho send
      case 1 : FFTresult = rho_rs.rhoIRXC;  break;
      case 2 : FFTresult = rho_rs.rhoIRYC;  break;
      case 3 : FFTresult = rho_rs.rhoIRZC;  break;
      default: CkAbort("impossible iopt");  break;
     }//end switch
    }else{
     switch(iopt){   
      case 0 : FFTresult = rho_rs.rhoIRXCint;  break; //Scr variable for rho send
      case 1 : FFTresult = rho_rs.rhoIRXCint;  break;
      case 2 : FFTresult = rho_rs.rhoIRYCint;  break;
      case 3 : FFTresult = rho_rs.rhoIRZCint;  break;
      default: CkAbort("impossible iopt");  break;
     }//end switch
    }//endif

//============================================================================
// Commlib launch

#ifdef USE_COMLIB
    if(rhoRsubplanes==1){
      switch(iopt){
#ifdef OLD_COMMLIB
      case 0 : if(config.useRInsRhoGP)    commRealInstance.beginIteration();    break;
      case 1 : if(config.useRInsIGXRhoGP) commRealIGXInstance.beginIteration(); break;
      case 2 : if(config.useRInsIGYRhoGP) commRealIGYInstance.beginIteration(); break;
      case 3 : if(config.useRInsIGZRhoGP) commRealIGZInstance.beginIteration(); break;
#else
      case 0 : if(config.useRInsRhoGP)    ComlibBegin(rhoGProxy_com);    break;
      case 1 : if(config.useRInsIGXRhoGP) ComlibBegin(rhoGProxyIGX_com); break;
      case 2 : if(config.useRInsIGYRhoGP) ComlibBegin(rhoGProxyIGY_com); break;
      case 3 : if(config.useRInsIGZRhoGP) ComlibBegin(rhoGProxyIGZ_com); break;
#endif
       default: CkAbort("impossible iopt");break;
     }//end switch
    }//endif
#endif

//============================================================================
// Send the data

    int iy = thisIndex.y;
    for(int ic = 0; ic < nchareRhoG; ic ++) { // chare arrays to which we will send

     //---------------------------
     //malloc the message
      int sendFFTDataSize = nlines_per_chareRhoG[ic];
      if(rhoRsubplanes!=1){sendFFTDataSize = nlines_per_chareRhoGY[ic][iy];}
      if(sendFFTDataSize>0)
	{
	  RhoGSFFTMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) RhoGSFFTMsg; 

	  //---------------------------
	  //Pack the message
	  msg->size        = sendFFTDataSize;
	  msg->iopt        = iopt;
	  msg->offset      = thisIndex.x;    // z-index
	  msg->offsetGx    = thisIndex.y;    // gx parallelization index
	  complex *data    = msg->data;
	  if(config.prioFFTMsg){
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	    *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.x;
	  }//endif

	  if(rhoRsubplanes==1){
	    for(int i=0;i<sendFFTDataSize;i++){
	      data[i] = FFTresult[tranpack_rho[ic][i]];
	    }//endfor
	  }else{
	    for(int i=0;i<sendFFTDataSize;i++){
	      data[i] = FFTresult[tranupack_rhoY[ic][iy][i]];
	    }//endfor
	  }//endif

	  //---------------------------
	  // Send the message
	  if(rhoRsubplanes==1){
	    switch(iopt){
	    case 0 : rhoGProxy_com(ic,0).acceptRhoData(msg);      break;
	    case 1 : rhoGProxyIGX_com(ic,0).acceptWhiteByrd(msg); break;
	    case 2 : rhoGProxyIGY_com(ic,0).acceptWhiteByrd(msg); break;
	    case 3 : rhoGProxyIGZ_com(ic,0).acceptWhiteByrd(msg); break;
	    default: CkAbort("impossible iopt");break;
	    }//end switch
	  }else{
	    switch(iopt){
	    case 0 : UrhoGProxy[thisInstance.proxyOffset](ic,0).acceptRhoData(msg);   break;
	    case 1 : UrhoGProxy[thisInstance.proxyOffset](ic,0).acceptWhiteByrd(msg); break;
	    case 2 : UrhoGProxy[thisInstance.proxyOffset](ic,0).acceptWhiteByrd(msg); break;
	    case 3 : UrhoGProxy[thisInstance.proxyOffset](ic,0).acceptWhiteByrd(msg); break;
	    default: CkAbort("impossible iopt");break;
	    }//end switch
	  }//endif
	}//end if nonzero
#ifdef CMK_BLUEGENEL
      if(ic%4==0)
	CmiNetworkProgress();
#endif

    }//end for : chare sending

//============================================================================
// Commlib stop

#ifdef USE_COMLIB
    if(rhoRsubplanes==1){
      switch(iopt){
#ifdef OLD_COMMLIB
      case 0 : if(config.useRInsRhoGP)    commRealInstance.endIteration(); break;
      case 1 : if(config.useRInsIGXRhoGP) commRealIGXInstance.endIteration(); break;
      case 2 : if(config.useRInsIGYRhoGP) commRealIGYInstance.endIteration(); break;
      case 3 : if(config.useRInsIGZRhoGP) commRealIGZInstance.endIteration(); break;
#else
      case 0 : if(config.useRInsRhoGP)    ComlibEnd(rhoGProxy_com);    break;
      case 1 : if(config.useRInsIGXRhoGP) ComlibEnd(rhoGProxyIGX_com); break;
      case 2 : if(config.useRInsIGYRhoGP) ComlibEnd(rhoGProxyIGY_com); break;
      case 3 : if(config.useRInsIGZRhoGP) ComlibEnd(rhoGProxyIGZ_com); break;
#endif
      default: CkAbort("impossible iopt");break;
      }//end switch
    }//endif
#endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::sendPartlyFFTtoRhoGall(){
//============================================================================
// Local pointers and variables

    CPcharmParaInfo *sim        = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    int nchareRhoG              = sim->nchareRhoG;
    int **tranpack_rho          = sim->index_tran_upack_rho;
    int *nlines_per_chareRhoG   = sim->nlines_per_chareRhoG;
    int ***tranupack_rhoY       = sim->index_tran_upack_rho_y;
    int **nlines_per_chareRhoGY = sim->nline_send_rho_y;

    complex *FFTresultX;
    complex *FFTresultY;
    complex *FFTresultZ;
    if(rhoRsubplanes==1){
      FFTresultX = rho_rs.rhoIRXC;
      FFTresultY = rho_rs.rhoIRYC;
      FFTresultZ = rho_rs.rhoIRZC;
    }else{
      FFTresultX = rho_rs.rhoIRXCint;
      FFTresultY = rho_rs.rhoIRYCint;
      FFTresultZ = rho_rs.rhoIRZCint;
    }//endif

//============================================================================
// Send the data

    int iy = thisIndex.y;
    for(int ic = 0; ic < nchareRhoG; ic ++) { // chare arrays to which we will send

     //---------------------------
     //malloc the message
      int sendFFTDataSize = nlines_per_chareRhoG[ic];
      if(rhoRsubplanes!=1){sendFFTDataSize = nlines_per_chareRhoGY[ic][iy];}

      if(sendFFTDataSize>0){
        //---------------------------
        //Pack the message
	  RhoGSFFTMsg *msg = new (3*sendFFTDataSize, 8 * sizeof(int)) RhoGSFFTMsg; 
	  msg->size        = sendFFTDataSize;
	  msg->iopt        = 1;
	  msg->offset      = thisIndex.x;    // z-index
	  msg->offsetGx    = thisIndex.y;    // gx parallelization index
	  complex *data    = msg->data;
	  if(config.prioFFTMsg){
	    CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	    *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.x;
	  }//endif

	  if(rhoRsubplanes==1){
	    for(int i=0,ii=0;ii<sendFFTDataSize;i+=3,ii++){
              int j = tranpack_rho[ic][ii];
	      data[i]   = FFTresultX[j];
	      data[i+1] = FFTresultY[j];
	      data[i+2] = FFTresultZ[j];
	    }//endfor
	  }else{
	    for(int i=0,ii=0;ii<sendFFTDataSize;i+=3,ii++){
              int j = tranupack_rhoY[ic][iy][ii];
	      data[i]   = FFTresultX[j];
	      data[i+1] = FFTresultY[j];
	      data[i+2] = FFTresultZ[j];
	    }//endfor
	  }//endif
	 //---------------------------
	 // Send the message
          UrhoGProxy[thisInstance.proxyOffset](ic,0).acceptWhiteByrdAll(msg);
      }//end if : nonzero msg
#ifdef CMK_BLUEGENEL
      if(ic%4==0){CmiNetworkProgress();}
#endif
   }//end for : send to rhog(ic)

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *
 * Accept transpose data from RhoG : receive grad_rho(gy,gx,z)
 *
 * Invoked 3 times per algorithm step : once for each grad_rho
 * 
 * Memory required is : rho_igx,rho_igy,rho_igz so stuff can come in any order
 *                    : density and vks are needed later so no reusing for you.
 *                    : VksHart can also arrive at any time and cannot be used
 *                      here.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptGradRhoVks(RhoRSFFTMsg *msg){
//============================================================================
// Unpack the message

  CPcharmParaInfo *sim        = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG                 = sim->nchareRhoG;
  int **tranUnpack            = sim->index_tran_upack_rho;
  int *nlines_per_chareG      = sim->nlines_per_chareRhoG;
  int ***tranupack_rhoY       = sim->index_tran_upack_rho_y;
  int **nlines_per_chareRhoGY = sim->nline_send_rho_y;
  int iy                      = thisIndex.y;

  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;

  int mySize;
  int nptsExpnd;
  if(rhoRsubplanes==1){
     mySize = nlines_per_chareG[Index];
     nptsExpnd = nptsExpndB;
  }else{
     mySize = nlines_per_chareRhoGY[Index][iy];
     nptsExpnd = nptsExpndA;
  }//endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
  }
#endif
#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("Data from RhoG arriving at RhoR : %d %d %d %d\n",
           thisIndex.x,thisIndex.y,iopt,countGradVks[iopt]);
#endif

//============================================================================
// Perform some error checking

  countGradVks[iopt]++;
  if (countGradVks[iopt] > recvCountFromGRho) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               countGradVks[iopt],nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=mySize){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.1, %d != %d for rho chare %d %d %d\n",size,mySize,
                  thisIndex.y,Index,iopt);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(1> iopt || iopt >3){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Wrong option in rhoR \n",iopt);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message

  complex *dataC;
  double  *dataR;
  if(rhoRsubplanes==1){
   switch(iopt){
     case 1 : dataC = rho_rs.rhoIRXC; dataR = rho_rs.rhoIRX; break;
     case 2 : dataC = rho_rs.rhoIRYC; dataR = rho_rs.rhoIRY; break;
     case 3 : dataC = rho_rs.rhoIRZC; dataR = rho_rs.rhoIRZ; break;
     default: CkAbort("impossible iopt");break;
   }//end switch
  }else{
   switch(iopt){
     case 1 : dataC = rho_rs.rhoIRXCint; dataR = rho_rs.rhoIRXint; break;
     case 2 : dataC = rho_rs.rhoIRYCint; dataR = rho_rs.rhoIRYint; break;
     case 3 : dataC = rho_rs.rhoIRZCint; dataR = rho_rs.rhoIRZint; break;
     default: CkAbort("impossible iopt");break;
   }//end switch
  }//endif

  // you must zero because messages don't fill the plane
  if(countGradVks[iopt]==1){bzero(dataR,sizeof(double)*nptsExpnd);}

  if(rhoRsubplanes==1){
    for(int i=0;i<size;i++){ 
      dataC[tranUnpack[Index][i]] = partiallyFFTd[i]*probScale;
    }//endfor
  }else{
    for(int i=0;i<size;i++){
      dataC[tranupack_rhoY[Index][iy][i]] = partiallyFFTd[i]*probScale;
    }//endfor
  }//endif

  delete msg;

//============================================================================
// When you have all the data : finish the FFT back to real space

  if (countGradVks[iopt] == recvCountFromGRho){

    countGradVks[iopt]=0;
    if(rhoRsubplanes==1){doneGradRhoVks++;}

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

    FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
    if(rhoRsubplanes==1){
      fftcache->doRhoFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb,iplane_ind);
    }else{
      fftcache->doRhoFFTGyToRy_Rchare(dataC,dataR,myNplane_rho,ngrida,ngridb,iplane_ind);
      sendPartlyFFTGxToRx(iopt);
    }//endif

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(fwFFTGtoRnot0_, StartTime, CmiWallTimer());    
#endif

  }//endif : I captured a divRho

//============================================================================
// When you have rhoiRX,rhoiRY,rhoiRZ and Vks invoke gradient correction

  if(doneGradRhoVks==3 && rhoRsubplanes==1){ 
    doneGradRhoVks = 0;
    GradCorr();         // if rhosubplanes>1 you have a transpose 
                        // and another fft to do before you can GradCorr);
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *
 * Accept transpose data from RhoG : receive grad_rho(gy,gx,z)
 *
 * Invoked 3 times per algorithm step : once for each grad_rho
 * 
 * Memory required is : rho_igx,rho_igy,rho_igz so stuff can come in any order
 *                    : density and vks are needed later so no reusing for you.
 *                    : VksHart can also arrive at any time and cannot be used
 *                      here.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptGradRhoVksAll(RhoRSFFTMsg *msg){
//============================================================================
// Unpack the message

  CPcharmParaInfo *sim        = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG                 = sim->nchareRhoG;
  int **tranUnpack            = sim->index_tran_upack_rho;
  int *nlines_per_chareG      = sim->nlines_per_chareRhoG;
  int ***tranupack_rhoY       = sim->index_tran_upack_rho_y;
  int **nlines_per_chareRhoGY = sim->nline_send_rho_y;
  int iy                      = thisIndex.y;

  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;

  int mySize;
  int nptsExpnd;
  if(rhoRsubplanes==1){
     mySize = nlines_per_chareG[Index];
     nptsExpnd = nptsExpndB;
  }else{
     mySize = nlines_per_chareRhoGY[Index][iy];
     nptsExpnd = nptsExpndA;
  }//endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
  }
#endif

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("Data from RhoG arriving at RhoR : %d %d %d %d\n",
           thisIndex.x,thisIndex.y,iopt,countGradVks[iopt]);
#endif

//============================================================================
// Perform some error checking

  countGradVks[iopt]++;

  if (countGradVks[iopt] > recvCountFromGRho) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               countGradVks[iopt],nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=mySize){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.1, %d != %d for rho chare %d %d %d\n",size,mySize,
                  thisIndex.y,Index,iopt);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(iopt!=1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Wrong option in rhoR \n",iopt);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message

  complex *dataCX,*dataCY,*dataCZ;
  double  *dataRX,*dataRY,*dataRZ;
  if(rhoRsubplanes==1){
     dataCX = rho_rs.rhoIRXC;    dataRX = rho_rs.rhoIRX;
     dataCY = rho_rs.rhoIRYC;    dataRY = rho_rs.rhoIRY;
     dataCZ = rho_rs.rhoIRZC;    dataRZ = rho_rs.rhoIRZ;
  }else{
     dataCX = rho_rs.rhoIRXCint; dataRX = rho_rs.rhoIRXint;
     dataCY = rho_rs.rhoIRYCint; dataRY = rho_rs.rhoIRYint;
     dataCZ = rho_rs.rhoIRZCint; dataRZ = rho_rs.rhoIRZint;
  }//endif

  // you must zero because messages don't fill the plane
  if(countGradVks[iopt]==1){
     bzero(dataRX,sizeof(double)*nptsExpnd);
     bzero(dataRY,sizeof(double)*nptsExpnd);
     bzero(dataRZ,sizeof(double)*nptsExpnd);
  }//endif

  if(rhoRsubplanes==1){
    for(int i=0,ii=0;ii<size;i+=3,ii++){ 
      int j = tranUnpack[Index][ii];
      dataCX[j] = partiallyFFTd[i]*probScale;
      dataCY[j] = partiallyFFTd[i+1]*probScale;
      dataCZ[j] = partiallyFFTd[i+2]*probScale;
    }//endfor
  }else{
    for(int i=0,ii=0;ii<size;i+=3,ii++){ 
      int j = tranupack_rhoY[Index][iy][ii];
      dataCX[j] = partiallyFFTd[i]*probScale;
      dataCY[j] = partiallyFFTd[i+1]*probScale;
      dataCZ[j] = partiallyFFTd[i+2]*probScale;
    }//endfor
  }//endif

  delete msg;

//============================================================================
// When you have all the data : finish the FFT back to real space

  if (countGradVks[iopt] == recvCountFromGRho){

    countGradVks[iopt]=0;
    if(rhoRsubplanes==1){doneGradRhoVks+=3;}

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

    FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
    if(rhoRsubplanes==1){
      fftcache->doRhoFFTGtoR_Rchare(dataCX,dataRX,nplane_rho_x,ngrida,ngridb,iplane_ind);
      fftcache->doRhoFFTGtoR_Rchare(dataCY,dataRY,nplane_rho_x,ngrida,ngridb,iplane_ind);
      fftcache->doRhoFFTGtoR_Rchare(dataCZ,dataRZ,nplane_rho_x,ngrida,ngridb,iplane_ind);
    }else{
      fftcache->doRhoFFTGyToRy_Rchare(dataCX,dataRX,myNplane_rho,ngrida,ngridb,iplane_ind);
      fftcache->doRhoFFTGyToRy_Rchare(dataCY,dataRY,myNplane_rho,ngrida,ngridb,iplane_ind);
      fftcache->doRhoFFTGyToRy_Rchare(dataCZ,dataRZ,myNplane_rho,ngrida,ngridb,iplane_ind);
      sendPartlyFFTGxToRx(1);
      sendPartlyFFTGxToRx(2);
      sendPartlyFFTGxToRx(3);
    }//endif

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(fwFFTGtoRnot0_, StartTime, CmiWallTimer());    
#endif

  }//endif : I captured a divRho

//============================================================================
// When you have rhoiRX,rhoiRY,rhoiRZ and Vks invoke gradient correction

  if(doneGradRhoVks==3 && rhoRsubplanes==1){ 
    doneGradRhoVks = 0;
    GradCorr();         // if rhosubplanes>1 you have a transpose 
                        // and another fft to do before you can GradCorr);
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// Double Transpose Bck Send :  A(gx,y,z) on the way to A(x,y,z)
//                              Send so (gx,z) parallel -> (y,z) parallel 
//
// Invoked 5 times an algorithm step:
//             case 1-3: Gradients
//             case  0 : VksWhiteByrd
//             case  4 : VksHarteext
//
//============================================================================
void CP_Rho_RealSpacePlane::sendPartlyFFTGxToRx(int iopt){
//============================================================================
  CkAssert(rhoRsubplanes>1);

    complex *FFTresult;
    switch(iopt){   
      case 0 : FFTresult = rho_rs.rhoIRXCint; break;
      case 1 : FFTresult = rho_rs.rhoIRXCint; break;
      case 2 : FFTresult = rho_rs.rhoIRYCint; break;
      case 3 : FFTresult = rho_rs.rhoIRZCint; break;
      case 4 : FFTresult = rho_rs.VksHartCint; break;
      default: CkAbort("impossible iopt");  break;
    }//end switch

//============================================================================
// Launch the communication

  //-----------------------------------------------------------------------------
  // Commlib launch : 

#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){
       case 0 : if(config.useRInsRhoRP)    commRealInstanceRx.beginIteration();    break;
       case 1 : if(config.useRInsIGXRhoRP) commRealIGXInstanceRx.beginIteration(); break;
       case 2 : if(config.useRInsIGYRhoRP) commRealIGYInstanceRx.beginIteration(); break;
       case 3 : if(config.useRInsIGZRhoRP) commRealIGZInstanceRx.beginIteration(); break;
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
      RhoGSFFTMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) RhoGSFFTMsg; 
      msg->size        = size;
      msg->iopt        = iopt;
      msg->offset      = thisIndex.y;   // my chare index
      msg->num         = myNplane_rho;  // number of gx-lines I have. 
      complex *data    = msg->data;     // data

      if(config.prioFFTMsg){
          CkSetQueueing(msg, CK_QUEUEING_IFIFO);
          *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.y;
      }//endif

      for(int i=ist,koff=0;i<iend;i++,koff+=myNplane_rho){ 
        for(int k=koff,ii=i;k<myNplane_rho+koff;k++,ii+=ngridb){
          data[k] = FFTresult[ii]; 
        }//endfor
      }//endfor

      switch(iopt){
        case 0 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksGxToRx(msg);break;
        case 1 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksGxToRx(msg); break;
        case 2 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksGxToRx(msg); break;
        case 3 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksGxToRx(msg); break;
        case 4 : UrhoRealProxy[thisInstance.proxyOffset](ix,ic).acceptRhoGradVksGxToRx(msg); break;
        default: CkAbort("impossible iopt");break;
      }//end switch

#ifdef _ERIC_SETS_UP_COMMLIB_
      switch(iopt){
        case 0 : UrhoGProxy[thisInstance.proxyOffset]_com(ic,0).acceptRhoGradVksGxToRx(msg);break;
        case 1 : UrhoGProxy[thisInstance.proxyOffset]IGX_com(ic,0).acceptRhoGradVksGxToRx(msg); break;
        case 2 : UrhoGProxy[thisInstance.proxyOffset]IGY_com(ic,0).acceptRhoGradVksGxToRx(msg); break;
        case 3 : UrhoGProxy[thisInstance.proxyOffset]IGZ_com(ic,0).acceptRhoGradVksGxToRx(msg); break;
        default: CkAbort("impossible iopt");break;
      }//end switch
#endif

#ifdef CMK_BLUEGENEL
      if(ic%3==0)
       CmiNetworkProgress();
#endif
    }//end for : chare sending

  //-----------------------------------------------------------------------------
  // Commlib stop

#ifdef _ERIC_SETS_UP_COMMLIB_
    switch(iopt){
       case 0 : if(config.useRInsRhoRP)    commRealInstanceRx.endIteration();    break;
       case 1 : if(config.useRInsIGXRhoRP) commRealIGXInstanceRx.endIteration(); break;
       case 2 : if(config.useRInsIGYRhoRP) commRealIGYInstanceRx.endIteration(); break;
       case 3 : if(config.useRInsIGZRhoRP) commRealIGZInstanceRx.endIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch
#endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// Double Transpose Bck Recv :  A(gx,y,z) on the way to A(x,y,z)
//                              Recv (gx,z) parallel -> (y,z) parallel 
//
// Invoked 5 times an algorithm step:
//             case 1-3: Gradients
//             case  0 : VksWhiteByrd
//             case  4 : VksHarteext
//
//============================================================================
void CP_Rho_RealSpacePlane::acceptRhoGradVksGxToRx(RhoGSFFTMsg *msg){
//============================================================================

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int **listSubGx      = sim->listSubGx;
  int  *numSubGx       = sim->numSubGx;

  int size         = msg->size;  // msg size
  int iopt         = msg->iopt;  
  int num          = msg->num;    // number of lines along `a' sent
  int offset       = msg->offset; // chare array that sent the data
  complex *msgData = msg->data;

  CkAssert(size==myNgridb*num);
  CkAssert(rhoRsubplanes>1);

  complex *dataC;
  double  *dataR;
  switch(iopt){
    case  0: dataC = rho_rs.densityC; dataR = rho_rs.density; break;
    case  1: dataC = rho_rs.rhoIRXC;  dataR = rho_rs.rhoIRX; break;
    case  2: dataC = rho_rs.rhoIRYC;  dataR = rho_rs.rhoIRY; break;
    case  3: dataC = rho_rs.rhoIRZC;  dataR = rho_rs.rhoIRZ; break;
    case  4: dataC = rho_rs.VksHartC; dataR = rho_rs.VksHart; break;
    default: CkAbort("Impossible option\n"); break;
  }//endif

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

  int done = 0;
  if(countIntGtoR[iopt]==rhoRsubplanes){
    done = 1;
    countIntGtoR[iopt]=0;
    FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    fftcache->doRhoFFTGxToRx_Rchare(dataC,dataR,nplane_rho_x,ngrida,myNgridb,iplane_ind);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(doRhoFFTGxtoRx_, StartTime, CmiWallTimer());    
#endif

  }//endif

//============================================================================
// When you have completed an FFT, you have some choices of what to do next

  if(done == 1){ 

   //---------------------------------------------------
   // If you have finished a gradient, increment counter.
   // If you have completed all gradients : Gradcorr.

    if(1 <= iopt && iopt <=3){
      doneGradRhoVks++;
      if(doneGradRhoVks==3){
       doneGradRhoVks=0;
       GradCorr();
      }//endif
    }//endif

   //---------------------------------------------------
   // if you have finished the whiteByrd : add it to vks

    if(iopt==0){addWhiteByrdVks();}

   //---------------------------------------------------
   // if you have finished the HartEext : add it to vks

    if(iopt==4){addHartEextVks();}

  }//endif : you have completed an fft

//============================================================================
  }//end routine
//============================================================================




//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// The gradient of the density is now completed. You can compute the 
// GGA-DFT functional now. Density is now available to be used as scratch.
//
// Invoked once per algorithm step
//
//============================================================================
void CP_Rho_RealSpacePlane::GradCorr(){
//============================================================================

  if(cp_grad_corr_on==0){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Don't come in the grad corr routines when\n");
     CkPrintf("gradient corrections are off\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

   double *density            = rho_rs.density;
   double *Vks                = rho_rs.Vks;
   double *rhoIRX             = rho_rs.rhoIRX;
   double *rhoIRY             = rho_rs.rhoIRY;
   double *rhoIRZ             = rho_rs.rhoIRZ;

   int nf1                    = ngrida;
   int nf2                    = myNgridb;
   int nf3                    = ngridc;
   int npts         = config.numFFTPoints;  // total number of points
   double *exc_gga_ret        = &(rho_rs.exc_gga_ret);

//============================================================================

#ifdef _CP_DEBUG_RHOR_VKSC_
    char myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "BGradRho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
    FILE *fp = fopen(myFileName,"w");
      for (int i = 0; i <ngridb*ngrida; i++){
        fprintf(fp,"%g %g %g %g\n",rho_rs.rhoIRX[i],rho_rs.rhoIRY[i],rho_rs.rhoIRZ[i],
                                   rho_rs.Vks[i]);
      }//endfor
    fclose(fp);
#endif

//============================================================================
// Compute the gradient corrected functional : Density is toast after this.

    rho_rs.exc_gga_ret = 0.0;
#define GGA_ON
#ifdef GGA_ON

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
    CPXCFNCTS::CP_getGGAFunctional(npts,nf1,nf2,nf3,density,rhoIRX,rhoIRY,rhoIRZ,
                                   Vks,thisIndex.x,exc_gga_ret,config.nfreq_xcfnctl);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(GradCorrGGA_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_BLUEGENEL
     CmiNetworkProgress();
#endif

#endif // GGA ON

//============================================================================
// Reduce the exchange correlation energy

   double exc[2];
   exc[0]=rho_rs.exc_ret;
   exc[1]=rho_rs.exc_gga_ret;
   contribute(2*sizeof(double),exc,CkReduction::sum_double);

//============================================================================
// output

#ifdef _CP_DEBUG_RHOR_VKSD_
    myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "AGradRho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
    fp = fopen(myFileName,"w");
      for (int i = 0; i <ngridb*ngrida; i++){
        fprintf(fp,"%g %g %g %g\n",rho_rs.rhoIRX[i],rho_rs.rhoIRY[i],rho_rs.rhoIRZ[i],
                                   rho_rs.Vks[i]);
      }//endfor
    fclose(fp);
#endif

//============================================================================
// Start the white bird puppy : back fft of rhoirx, rhoiry, rhoirz

    whiteByrdFFT();

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//
// The white-bird term : First fwfft redefined delrho(r) to delrho(g)
// then send to RhoGspacePlane. RhoGspacePlane sends you back back another term.
// After this routine, rhoIRX, rhoIRY and rhoIRZ are `free'.
//
// Invoked once per algorithm step
//
//
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::whiteByrdFFT(){
//============================================================================
// Constants and pointers

   FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
   double  *rhoIRX    = rho_rs.rhoIRX;
   double  *rhoIRY    = rho_rs.rhoIRY;
   double  *rhoIRZ    = rho_rs.rhoIRZ;
   complex *rhoIRXC   = rho_rs.rhoIRXC;
   complex *rhoIRYC   = rho_rs.rhoIRYC;
   complex *rhoIRZC   = rho_rs.rhoIRZC;

   int ioptx          = 1;
   int iopty          = 2;
   int ioptz          = 3;

//============================================================================
// I) rhoIRX : Scale, Real to complex FFT, perform FFT, transpose

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

  rho_rs.scale(rhoIRX,FFTscale);

  if(rhoRsubplanes==1){
    fftcache->doRhoFFTRtoG_Rchare(rhoIRXC,rhoIRX,nplane_rho_x,ngrida,ngridb,iplane_ind);
    if(config.rhoGToRhoRMsgComb==0){sendPartlyFFTtoRhoG(ioptx);}
  }else{
    fftcache->doRhoFFTRxToGx_Rchare(rhoIRXC,rhoIRX,nplane_rho_x,ngrida,myNgridb,iplane_ind);
    sendPartlyFFTRyToGy(ioptx);// transpose and do the y-gy fft
  }//endif

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(WhiteByrdFFTX_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_BLUEGENEL
  CmiNetworkProgress();
#endif

//============================================================================
// II) rhoIRY : Scale, real to complex FFT, perform FFT, transpose

#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  rho_rs.scale(rhoIRY,FFTscale);

  if(rhoRsubplanes==1){
    fftcache->doRhoFFTRtoG_Rchare(rhoIRYC,rhoIRY,nplane_rho_x,ngrida,ngridb,iplane_ind);
    if(config.rhoGToRhoRMsgComb==0){sendPartlyFFTtoRhoG(iopty);}
  }else{
    fftcache->doRhoFFTRxToGx_Rchare(rhoIRYC,rhoIRY,nplane_rho_x,ngrida,myNgridb,iplane_ind);
    sendPartlyFFTRyToGy(iopty); // transpose and do the y-gy fft
  }//endif

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(WhiteByrdFFTY_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_BLUEGENEL
  CmiNetworkProgress();
#endif

//============================================================================
// III) rhoIRZ : Scale, real to complex FFT, perform FFT, transpose

#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  rho_rs.scale(rhoIRZ,FFTscale);
  if(rhoRsubplanes==1){
    fftcache->doRhoFFTRtoG_Rchare(rhoIRZC,rhoIRZ,nplane_rho_x,ngrida,ngridb,iplane_ind);
    if(config.rhoGToRhoRMsgComb==0){sendPartlyFFTtoRhoG(ioptz);}
  }else{
    fftcache->doRhoFFTRxToGx_Rchare(rhoIRZC,rhoIRZ,nplane_rho_x,ngrida,myNgridb,iplane_ind);
    sendPartlyFFTRyToGy(ioptz);// transpose and do the y-gy fft
  }//endif

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(WhiteByrdFFTZ_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_BLUEGENEL
  CmiNetworkProgress();
#endif

//============================================================================
// Send all 3 components at once

  if(rhoRsubplanes==1 && config.rhoGToRhoRMsgComb==1){
     sendPartlyFFTtoRhoGall();
  }//endif

//============================================================================
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// The white bird vks correction is returned from RhoG : VksW(gx,gy,z)
//                This routine recvs the transpose {gx,gy} to (gx,z)
//
// Invoked once per algorithm step
//
// After this routine, VksW(gx,gy,z) -> VksW(x,y,gz) and is added vksTot();
//
//============================================================================
void CP_Rho_RealSpacePlane::acceptWhiteByrd(RhoRSFFTMsg *msg){
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("WhiteByrd Data from RhoG arriving at RhoR : %d %d\n",
           thisIndex.x,thisIndex.y);
#endif

//============================================================================
// Local Pointers

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoG;
  int **tranUnpack       = sim->index_tran_upack_rho;
  int *nlines_per_chareG = sim->nlines_per_chareRhoG;
  int ***tranupack_rhoY       = sim->index_tran_upack_rho_y;
  int **nlines_per_chareRhoGY = sim->nline_send_rho_y;
  int iy                      = thisIndex.y;
 
  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;

  int mySize;
  int nptsExpnd;
  if(rhoRsubplanes==1){
     mySize    = nlines_per_chareG[Index];
     nptsExpnd = nptsExpndB;
  }else{
     mySize    = nlines_per_chareRhoGY[Index][iy];
     nptsExpnd = nptsExpndA;
  }//endif

//============================================================================
// Perform some error checking

  countWhiteByrd++;
  if (countWhiteByrd > recvCountFromGRho) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               countWhiteByrd,nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=mySize){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.2, %d != %d for rho chare %d %d\n",size,mySize,
               thisIndex.y,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message : The density is free for use as scr

  double *dataR;
  complex *dataC;
  if(rhoRsubplanes==1){
    dataR = rho_rs.density; 
    dataC = rho_rs.densityC; 
  }else{
    dataR = rho_rs.rhoIRXint;
    dataC = rho_rs.rhoIRXCint;
  }//endif

  // zero because input data does not fill the plane
  if(countWhiteByrd==1){bzero(dataR,sizeof(double)*nptsExpnd);}

  if(rhoRsubplanes==1){
    for(int i=0;i<size;i++){ 
      dataC[tranUnpack[Index][i]] = partiallyFFTd[i];
    }//endfor
  }else{
    for(int i=0;i<size;i++){
      dataC[tranupack_rhoY[Index][iy][i]] = partiallyFFTd[i];
    }//endfor
  }//endif

  delete msg;

//============================================================================
// When you have all the messages, do the last fft, and add in the correction

  if(countWhiteByrd == recvCountFromGRho){
    countWhiteByrd=0;

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
    if(rhoRsubplanes==1){
      fftcache->doRhoFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb,iplane_ind);
      addWhiteByrdVks();
    }else{
      fftcache->doRhoFFTGyToRy_Rchare(dataC,dataR,myNplane_rho,ngrida,ngridb,iplane_ind);
      sendPartlyFFTGxToRx(0);
    }
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(PostByrdfwFFTGtoR_, StartTime, CmiWallTimer());    
#endif

  }//endif : communication from rhog has all arrived safely

//============================================================================
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Add the VksWhiteByrd to VksTot : Set the done flag.
//
//============================================================================
void CP_Rho_RealSpacePlane::addWhiteByrdVks(){
//============================================================================
// Add the whitebyrd contrib to vks

    int nptsExpnd  = nptsExpndB;
    double *dataR  = rho_rs.density; 
    double *Vks    = rho_rs.Vks;
    for(int i=0;i<nptsExpnd;i++){Vks[i] -= dataR[i];}

//============================================================================
// Our we done yet?

    doneWhiteByrd = true;
    doMulticastCheck();

//============================================================================
   }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *
 * Accept hartExt transpose data : receive VksHartEext(gx,gy,z) gx,z is parallel.
 * Since the message can come at any time memory, VksHart has to be ready to
 * receive it.
 *
 * Invoked once per algorithm step.
 *
 * After the routine, VksHartEext(gx,gy,z)--> VksHarteext(x,y,z) and is added
 *                    to Vkstot(x,y,z)
 *
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptHartVks(RhoHartRSFFTMsg *msg){
//============================================================================
// Local pointers

  CPcharmParaInfo *sim        = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG                 = sim->nchareRhoG;
  int **tranUnpack            = sim->index_tran_upack_rho;
  int *nlines_per_chareG      = sim->nlines_per_chareRhoG;
  int ***tranupack_rhoY       = sim->index_tran_upack_eext_ys;
  int **nlines_per_chareRhoGY = sim->nline_send_eext_y;
  int iy                      = thisIndex.y;
   
  int size               = msg->size; 
  int IndexS             = msg->index;
  int Index              = msg->senderBigIndex;
  int istrt_lines        = msg->senderStrtLine;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;

  double  *dataR;
  complex *dataC;
  int nptsExpnd;
  int mySize;
  if(rhoRsubplanes==1){
    dataR     = rho_rs.VksHart;
    dataC     = rho_rs.VksHartC;
    nptsExpnd = nptsExpndB;
    mySize    = size;
  }else{
    dataR     = rho_rs.VksHartint;
    dataC     = rho_rs.VksHartCint;
    nptsExpnd = nptsExpndA;
    mySize    = nlines_per_chareRhoGY[IndexS][iy];
  }//endif
  
  CkAssert(size==mySize);

//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("Data from RhoG arriving at RhoR : %d %d %d %d\n",
           thisIndex.x,thisIndex.y,iopt,countGradVks[iopt]);
#endif

  CkAssert(iopt==0);
  if(size!=mySize){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.2, %d != %d for rho chare %d %d : %d\n",size,mySize,
               thisIndex.y,IndexS,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// Copy out the hart-eext contrib to vks : iopt==0

  countGradVks[iopt]++; 
  if(countGradVks[iopt]==1){bzero(dataR,sizeof(double)*nptsExpnd);}
   
  if(rhoRsubplanes==1){
    for(int i=0,j=istrt_lines;i<size;i++,j++){
      dataC[tranUnpack[Index][j]] = partiallyFFTd[i];
    }//endfor
  }else{
    for(int i=0;i<size;i++){
      dataC[tranupack_rhoY[IndexS][iy][i]] = partiallyFFTd[i];
    }//endfor
  }//endif

  delete msg;  

//============================================================================
// fft the puppy if you've got it all : only atmtyp index=1 of ghart sends rho_real

  if (countGradVks[iopt] == recvCountFromGHartExt){
      countGradVks[iopt]=0;

      FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif
      if(rhoRsubplanes==1){
        fftcache->doRhoFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb,iplane_ind);
        addHartEextVks();
      }else{
        fftcache->doRhoFFTGyToRy_Rchare(dataC,dataR,myNplane_rho,ngrida,ngridb,iplane_ind);
        sendPartlyFFTGxToRx(4);
      }//endif

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(fwFFTGtoR0_, StartTime, CmiWallTimer());    
#endif
  }//endif : communication from rhog 

  
//============================================================================
  }//end routine 
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Add the VksHartEext to VksTot : Set the done flag.
//
//============================================================================
void CP_Rho_RealSpacePlane::addHartEextVks(){
//============================================================================
// Add the whitebyrd contrib to vks

    int nptsExpnd  = nptsExpndB;
    double *dataR  = rho_rs.VksHart; 
    double *Vks    = rho_rs.Vks;
    for(int i=0;i<nptsExpnd;i++){Vks[i]+=dataR[i];}

//============================================================================
// Our we done yet?

    doneHartVks = true;
    doMulticastCheck();

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Under ees-eext Rhart chare reports its completion :  Set the done flag.
//
//============================================================================
void CP_Rho_RealSpacePlane::RHartReport(){
  countRHart++;
  if(countRHart==countRHartValue){
      doneRHart=true;
      doMulticastCheck();
      countRHart=0;
  }//endif
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// If all the parts of exc-eext-hart are done, invoke blast of vks to states
//
//============================================================================
  void CP_Rho_RealSpacePlane::doMulticastCheck(){
//============================================================================

  if(doneWhiteByrd && doneRHart && doneHartVks){doMulticast();}

//============================================================================
  }//end routine
//============================================================================



//============================================================================
// Send vks back to the states
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::doMulticast(){
//============================================================================

  if(!doneWhiteByrd || !doneHartVks || !doneRHart){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Flow of Control Error : Attempting to rho multicast\n");
    CkPrintf("without harteext or gradcorr (whitebyrd) \n");
    CkPrintf("without rhodensity or rhart \n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

 // overkill on the resetting
  countRHart=0;
  doneWhiteByrd    = false;
  doneHartVks      = false;
  doneRHart        = false;
    
//============================================================================
// Send vks back to the states in real space

   int dataSize    = ngrida*myNgridb;
   double *Vks     = rho_rs.Vks;

   if ((config.useGMulticast+config.useCommlibMulticast)!=1) {
     CkAbort("No multicast strategy\n");
   }//endif

   if (config.useGMulticast || config.useCommlibMulticast) {

      ProductMsg *msg = new (dataSize, 0) ProductMsg;

      // ADD new index here for Y
      msg->idx        = thisIndex.x;
      msg->datalen    = dataSize;
      msg->hops       = 0;
      msg->subplane   = thisIndex.y;
      double *dataR   = msg->data;

      rho_rs.uPackShrink(dataR,Vks); // down pack Vks for the send
//#define _CP_DEBUG_RHOR_VKSE_
#ifdef _CP_DEBUG_RHOR_VKSE_
       char myFileName[MAX_CHAR_ARRAY_LENGTH];
       sprintf(myFileName, "vks_rho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
       FILE *fp = fopen(myFileName,"w");
       for (int j = 0, iii=0; j <myNgridb; j++){
       for (int i = 0; i <ngrida; i++,iii++){
         fprintf(fp,"%d %d %g\n",i,j+myBoff,dataR[iii]);
       }}//endfor
       fclose(fp);
       UrhoRealProxy[thisInstance.proxyOffset](0,0).exitForDebugging();
#else
       /*
#ifdef OLD_COMMLIB
      if(config.useCommlibMulticast){mcastInstance.beginIteration();}
#else
      if(config.useCommlibMulticast){ComlibBegin(realSpaceSectionCProxy);}
#endif
       */
      if(config.useCommlibMulticast){
        realSpaceSectionCProxy.acceptProduct(msg);
      }else{
        realSpaceSectionProxy.acceptProduct(msg);
      }//enddif
      /*
#ifdef OLD_COMMLIB
     if(config.useCommlibMulticast){mcastInstance.endIteration();}
#else
     if(config.useCommlibMulticast){ComlibEnd(realSpaceSectionCProxy);}
#endif
      */
#ifdef _CP_SUBSTEP_TIMING_
     if(rhoKeeperId>0)
       {
	 double rhoend=CmiWallTimer();
	 contribute(sizeof(double), &rhoend, CkReduction::max_double,  CkCallback(CkIndex_TimeKeeper::collectEnd(NULL),0,TimeKeeperProxy),rhoKeeperId);
       }
#endif

#endif
   }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Glenn's RhoReal exit 
//============================================================================
void CP_Rho_RealSpacePlane::exitForDebugging(){
  countDebug++;  
  if(countDebug==(rhoRsubplanes*ngridc)){
    countDebug=0;
    CkPrintf("I am in the exitfordebuging rhoreal puppy. Bye-bye\n");
    CkExit();
  }//endif
}
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::ResumeFromSync(){
  /*
    if(config.useCommlibMulticast)
        ComlibResetSectionProxy(&realSpaceSectionCProxy);
    if(config.useRInsRhoGP)
        ComlibResetProxy(&UrhoGProxy[thisInstance.proxyOffset]_com);
  */
}
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//! return tru if input is power of 2
//============================================================================
bool is_pow2(int input){
    int y=0;
    for(int x=0;x<32;x++){
      y = 1<<x;
      if(y==input){return true;}
    }//endfor
    return false;
}
//============================================================================
