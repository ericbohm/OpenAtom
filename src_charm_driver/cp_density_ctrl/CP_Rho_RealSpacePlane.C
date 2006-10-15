//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_RealSpacePlane
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
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "../../include/debug_flags.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CP_State_RealParticlePlane realParticlePlaneProxy;
extern CProxy_CP_Rho_RealSpacePlane   rhoRealProxy;
extern CProxy_CP_Rho_GSpacePlane      rhoGProxy;
extern CProxy_CPcharmParaInfoGrp      scProxy;
extern CProxy_FFTcache                fftCacheProxy;
extern CProxy_CP_Rho_RHartExt         rhoRHartExtProxy;
extern CProxy_CP_State_GSpacePlane    gSpacePlaneProxy;

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
RTH_Routine_locals(CP_Rho_RealSpacePlane,run)
//---------------------------------------------------------------------------
RTH_Routine_code(CP_Rho_RealSpacePlane,run) {
//---------------------------------------------------------------------------
  while(1) {

    // 1st entry point is the constructor which invokes run() to put you here.
    RTH_Suspend(); 
    // 2nd entry point is acceptDensity(msg *) data from statereal
    c->acceptDensity();
    c->launchNLRealFFT();
    RTH_Suspend(); 
    // could insert the launchFFT for nonlocal real 

    // 3rd entry point is acceptGradRhoVks(RhoRSFFTMsg *) 
    if(c->cp_grad_corr_on!=0){
       c->GradCorr();
       RTH_Suspend(); 
    }//endif

    //    c->launchNLRealFFT(); candidate location

    // 4th entry point is acceptWhiteByrd(RhoRSFFTMsg *) || acceptHartVks.
    // Whichever arrives LAST calls resume (the last shall be first).
    c->doneRhoReal=true;
    if(c->doneRHart)
      c->doMulticast();

  } //end while not done
//--------------------------------------------------------------------------
   } RTH_Routine_end(CP_Rho_RealSpacePlane,run)
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::run () {
  run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_Rho_RealSpacePlane,run),this);
  RTH_Runtime_resume(run_thread);
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the states,
// adds them up and remembers them for the next stage in computation
//
//============================================================================
CP_Rho_RealSpacePlane::CP_Rho_RealSpacePlane(int xdim, size2d yzdim,bool _useCommlib, 
                                             int _ees_eext_on,int _ngridcEext)
//============================================================================
   {//begin routine
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
    CkPrintf("[%d %d] RhoR constructs \n",thisIndex.x, thisIndex.y);
#endif

//============================================================================
// Malloc Memory, Intialize counters, set constants

   // Get parameters from the globals/groups
    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    cp_grad_corr_on      = sim->cp_grad_corr_on;
    nplane_rho_x         = sim->nplane_rho_x;
    double vol           = sim->vol;
    int numFFT           = config.numFFTPoints;
    rhoGHelpers          = config.rhoGHelpers;

   // Set up the data class
    initRhoRealSlab(&rho_rs,xdim,yzdim[0],yzdim[1],thisIndex.x,thisIndex.y);

   // Set up constants
    ees_eext_on  = _ees_eext_on;
    ngridcEext   = _ngridcEext;
    ngrida       = rho_rs.sizeX;
    ngridb       = rho_rs.sizeY;
    ngridc       = rho_rs.sizeZ;
    npts         =  ngrida*ngridb;
    nptsExpnd    = (ngrida+2)*ngridb;
    
    FFTscale     = 1.0/((double)numFFT);
    volumeFactor = vol*FFTscale;
    probScale    = 1.0/vol;

   // Initialize counters, set booleans.
    countWhiteByrd = 0;
    doneGradRhoVks = 0;
    countRHart     = 0;
    countRHartValue = 1;
    if( thisIndex.x   < ngridcEext - rho_rs.sizeZ)
      countRHartValue=2;
    doneHartVks    = true;
    doneWhiteByrd  = true;
    doneRHart      = false;
    doneRhoReal    = false;
    for(int i=0;i<4;i++){countGradVks[i]=0;}

//============================================================================
// make sections in the realSpacePlane array. These will be used when 
// computing real-space densities and multicasting v_ks values 

    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    CkArrayIndexMax *elems   = new CkArrayIndexMax[nstates];

    int j;
    // section i has al the portions with all 
    CkArrayIndex2D idx(0, thisIndex.x);
    if(is_pow2(nstates)){
	for (j = 0; j < nstates; j++) {
	    idx.index[0] = j^(thisIndex.x%nstates);
	    elems[j] = idx;
	}//endfor
    }else{
	for (j = 0; j < nstates; j++) {
	    idx.index[0] = (j+thisIndex.x)%nstates;
	    elems[j] = idx;
	}//endfor
    }//endif


    realSpaceSectionProxy = CProxySection_CP_State_RealSpacePlane::
        ckNew(realSpacePlaneProxy.ckGetArrayID(), elems, nstates);

    realSpaceSectionCProxy = CProxySection_CP_State_RealSpacePlane::
        ckNew(realSpacePlaneProxy.ckGetArrayID(), elems, nstates);

    realSpaceSectionProxy.ckDelegate
        (CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());

    mcastGrp->setSection(realSpaceSectionProxy);

    ComlibAssociateProxy(&mcastInstance,realSpaceSectionCProxy);

    delete [] elems;    

    ProductMsg *dummyProductMessage = new (0) ProductMsg;    
    // inform realspace element of this section proxy.
    dummyProductMessage->subplane=thisIndex.y;
    realSpaceSectionProxy.init(dummyProductMessage);

    rhoGProxy_com    = rhoGProxy;
    rhoGProxyIGX_com = rhoGProxy;
    rhoGProxyIGY_com = rhoGProxy;
    rhoGProxyIGZ_com = rhoGProxy;
    if (config.useRInsRhoGP) 
	ComlibAssociateProxy(&commRealInstance,rhoGProxy_com);          
    if (config.useRInsIGXRhoGP) 
	ComlibAssociateProxy(&commRealIGXInstance,rhoGProxyIGX_com);          
    if (config.useRInsIGYRhoGP) 
	ComlibAssociateProxy(&commRealIGYInstance,rhoGProxyIGY_com);          
    if (config.useRInsIGZRhoGP) 
	ComlibAssociateProxy(&commRealIGZInstance,rhoGProxyIGZ_com);          

//============================================================================
// Migration

    usesAtSync = CmiTrue;
   //    if(config.lbdensity){
   //      setMigratable(true);
   //    }else{
    setMigratable(false);
   //    }//endif

//============================================================================
// put yourself into the threaded loop

    run();

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_RealSpacePlane::~CP_Rho_RealSpacePlane(){
}
//============================================================================


//============================================================================
//  Pup my variables for migration
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);
  p|nplane_rho_x;
  p|ngrida;
  p|ngridb;
  p|ngridc;
  p|nptsExpnd;
  p|npts;
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

  PUParray(p,countGradVks,5);
  rho_rs.pup(p);   // pup your data class

 // Pupping Proxies???
  p|realSpaceSectionProxy;
  p|realSpaceSectionCProxy;
  if(p.isUnpacking()){
    rhoGProxy_com = rhoGProxy;
    rhoGProxyIGX_com = rhoGProxy;
    rhoGProxyIGY_com = rhoGProxy;
    rhoGProxyIGZ_com = rhoGProxy;
    if (config.useRInsRhoGP) 
	ComlibAssociateProxy(&commRealInstance,rhoGProxy_com);          
    if (config.useRInsIGXRhoGP) 
	ComlibAssociateProxy(&commRealIGXInstance,rhoGProxyIGX_com);          
    if (config.useRInsIGYRhoGP) 
	ComlibAssociateProxy(&commRealIGYInstance,rhoGProxyIGY_com);          
    if (config.useRInsIGZRhoGP) 
	ComlibAssociateProxy(&commRealIGZInstance,rhoGProxyIGZ_com);          

    }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
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
//============================================================================
// Set the flags

   if(cp_grad_corr_on!=0){
     doneWhiteByrd  = false;
   }else{
     doneWhiteByrd  = true;
   }//endif
   doneHartVks      = false;
   doneRHart        = false;
   doneRhoReal      = false;

   if(ees_eext_on==1){doneRhart = true;}
#ifdef _CP_DEBUG_HARTEEXT_OFF_
   doneHartVks    = true;
   doneRhart      = true;
#endif



#ifdef _NAN_CHECK_
  for(int i=0;i<msg->getSize()/sizeof(double) ;i++)
    {
      CkAssert(isnan(((double*) msg->getData())[i])==0);
    }
#endif

//============================================================================
// Unpack into spread out form and delete the message

    double *realValues = (double *) msg->getData(); 
    double *density    = rho_rs.density;
    CkAssert(msg->getSize() == ngridb * ngrida * sizeof(double));

    rho_rs.uPackScaleGrow(density,realValues,probScale);

    delete msg;

//============================================================================

#ifdef _CP_DEBUG_RHOR_RHO_
      char myFileName[MAX_CHAR_ARRAY_LENGTH];
      sprintf(myFileName, "Rho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
      FILE *fp = fopen(myFileName,"w");
        for (int i = 0; i <ngridb*ngrida; i++){
          fprintf(fp,"%g\n",realValues[i]*probScale);
        }//endfor
      fclose(fp);
#endif

//============================================================================
// the threaded loop will throw you down to the next routine

    RTH_Runtime_resume(run_thread);

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Now that the density is here : 1) Launch ees routines 2) compute EEXC
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity() {
//============================================================================
// Launch the external energy computation in r-space
 
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
    rhoRHartExtProxy(thisIndex.x,0).startEextIter();
    if(thisIndex.x<rem){
#ifdef _CP_RHO_RSP_VERBOSE_
      CkPrintf("HI, I am r-rho chare %d also lauchning %d\n",thisIndex.x,ind);
#endif
      rhoRHartExtProxy(ind,0).startEextIter();
    }//endif
#endif
    // launch nonlocal g space if it wasn't done in RS

    CPcharmParaInfo *sim      = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    if(sim->ees_nloc_on==1 && config.launchNLeesFromRho)
      { // kick off the NLeesG from here

	if(thisIndex.x<config.nchareG)
	  for(int ns=0;ns<config.nstates;ns++)
	    {
	      //	    CkPrintf("RhoRP[%d,%d] triggering NL %d %d \n",thisIndex.x, thisIndex.y, thisIndex.x, nchg);
	      gSpacePlaneProxy(ns,thisIndex.x).startNLEes(false);
	    }

      }
  }//endif


//============================================================================
// Compute the exchange correlation energy (density no-grad part)

  energyComputation();

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 *  launch nonlocal real no matter what
 *  this can be shifted to anywhere (outside grad corr which is optional)
 *
 */
//============================================================================
void CP_Rho_RealSpacePlane::launchNLRealFFT()
{
  CPcharmParaInfo *sim      = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  if(sim->ees_nloc_on==1)
    { // tell NLeesR its ok to compute
      if(thisIndex.x< sim->ngrid_nloc_c)
	for(int ns=0;ns<config.nstates;ns++)
	  {
	    realParticlePlaneProxy(ns,thisIndex.x).launchFFTControl();
	  }
      
    }
  //----------------------------------------------------------------------------
} //end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * Compute one part of the EXC energy using PINY CP_exc_calc
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
   int nf2          = ngridb;
   int nf3          = ngridc;
   double *exc_ret  = &(rho_rs.exc_ret);
   double *muxc_ret = &(rho_rs.muxc_ret);

//============================================================================
// Perform exchange correlation computation (no grad corr here).

   CPXCFNCTS::CP_exc_calc(npts,nf1,nf2,nf3,density,Vks,exc_ret,muxc_ret);

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

  if(cp_grad_corr_on==0){
    double exc[2];
    exc[0]=rho_rs.exc_ret;
    exc[1]=0.0;
    contribute(2*sizeof(double),exc,CkReduction::sum_double);
  }//endif

//============================================================================
// Invoke FFT to take rho(r) to rho(g)

   fftRhoRtoRhoG();

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *
 * Perform fft on density to get rho(g), part way
 *   FFT Rho(x,y,z) ----> Rho(gx,gy,z) : store result in rhoIRX which is idle.
 *
 */
//============================================================================
void CP_Rho_RealSpacePlane::fftRhoRtoRhoG(){
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("In RhoRealSpacePlane[%d %d] FFT_RSpacetoGSpace %d %d\n",thisIndex.x, 
                   thisIndex.y, CkMyPe(),CmiMemoryUsage());
#endif

//============================================================================
// FFT myself back to G-space part way

  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
  double  *density   = rho_rs.density;  // we need to save the density and vks.
  double  *dataR     = rho_rs.rhoIRX;   // rhoirx is around doing nothing now
  complex *dataC     = rho_rs.rhoIRXC;  // so we can use it to store the FFT

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

  // SPLIT this 2D XY FFT to just X

  rho_rs.uPackScale(dataR,density,volumeFactor);  // can't avoid the memcpy-scaling
  fftcache->doRhoFFTRtoG_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(RhoRtoGFFT_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

   // Add transpose and communicate

   // FFT along Y
    
   /// now you can send to RhoG

//============================================================================
// Send chunk to RhoGDensity to complete the FFT

 int iopt = 0;
 sendPartlyFFTtoRhoG(iopt);

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// Send the partly FFTed density, Rho(gx,gy,z), to rhoGSpacePlane :  
// You then wait while RHOG works, sends back RhoIRalpha which you whitebyrdize
// You send back the whitebyrdized RhoIRalpha puppies to RhoGspacePlane.
// Whilst all this is going on, HartG is churning and will send you another
// part of Vks.  All this requires keeping density, vks, div_rho_alpha and
// vkshart memory available at all times to receive messages.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::sendPartlyFFTtoRhoG(int iopt){
//============================================================================

    CPcharmParaInfo *sim      = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    int nchareRhoG            = sim->nchareRhoG;
    int **tranpack_rho        = sim->index_tran_upack_rho;
    int *nlines_per_chareRhoG = sim->nlines_per_chareRhoG;

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
  // Commlib launch
    switch(iopt){
       case 0 : if(config.useRInsRhoGP)    commRealInstance.beginIteration();    break;
       case 1 : if(config.useRInsIGXRhoGP) commRealIGXInstance.beginIteration(); break;
       case 2 : if(config.useRInsIGYRhoGP) commRealIGYInstance.beginIteration(); break;
       case 3 : if(config.useRInsIGZRhoGP) commRealIGZInstance.beginIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch

   //-----------------------------------------------------------------------------
   // Send the data
    for(int ic = 0; ic < nchareRhoG; ic ++) { // chare arrays to which we will send

      int sendFFTDataSize = nlines_per_chareRhoG[ic];
      RhoGSFFTMsg *msg = new (sendFFTDataSize, 8 * sizeof(int)) RhoGSFFTMsg; 
 
      msg->size        = sendFFTDataSize;
      msg->iopt        = iopt;
      msg->offset      = thisIndex.x;    // z-index
      complex *data    = msg->data;
      if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.rhogpriority+thisIndex.x;
      }//endif

      for(int i=0;i<sendFFTDataSize;i++){
        data[i] = FFTresult[tranpack_rho[ic][i]];
      }//endfor

      switch(iopt){
        case 0 : rhoGProxy_com(ic,0).acceptRhoData(msg);      break;
        case 1 : rhoGProxyIGX_com(ic,0).acceptWhiteByrd(msg); break;
        case 2 : rhoGProxyIGY_com(ic,0).acceptWhiteByrd(msg); break;
        case 3 : rhoGProxyIGZ_com(ic,0).acceptWhiteByrd(msg); break;
        default: CkAbort("impossible iopt");break;
      }//end switch

#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
    }//end for : chare sending

  //-----------------------------------------------------------------------------
  // Commlib stop
    switch(iopt){
       case 0 : if(config.useRInsRhoGP)    commRealInstance.endIteration(); break;
       case 1 : if(config.useRInsIGXRhoGP) commRealIGXInstance.endIteration(); break;
       case 2 : if(config.useRInsIGYRhoGP) commRealIGYInstance.endIteration(); break;
       case 3 : if(config.useRInsIGZRhoGP) commRealIGZInstance.endIteration(); break;
       default: CkAbort("impossible iopt");break;
    }//end switch

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
/**
 *
 * accept tranpose data : receive grad_rho(z,gy,gx) from RHOGSPACE
 * Memory required is : rho_igx,rho_igy,rho_igz so stuff can come in any order
 *                    : density and vks are needed later so no reusing for you.
 *                    : VksHart can also arrive at any time and cannot be used
 *                      here.
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptGradRhoVks(RhoRSFFTMsg *msg){
//============================================================================
// Unpack the message

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoG;
  int **tranUnpack       = sim->index_tran_upack_rho;
  int *nlines_per_chareG = sim->nlines_per_chareRhoG;

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
    }
#endif
   
  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("Data from RhoG arriving at RhoR : %d %d %d %d\n",
	   thisIndex.x,thisIndex.y,iopt,countGradVks[iopt]);
#endif

//============================================================================
// Perform some error checking

  countGradVks[iopt]++;
  if (countGradVks[iopt] > nchareG) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               countGradVks[iopt],nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=nlines_per_chareG[Index]){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.1, %d != %d for rho chare %d %d %d\n",size,nlines_per_chareG[Index],
                  thisIndex.y,Index,iopt);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message

  complex *dataC;
  double  *dataR;
  switch(iopt){
     case 1 : dataC = rho_rs.rhoIRXC; dataR = rho_rs.rhoIRX; break;
     case 2 : dataC = rho_rs.rhoIRYC; dataR = rho_rs.rhoIRY; break;
     case 3 : dataC = rho_rs.rhoIRZC; dataR = rho_rs.rhoIRZ; break;
     default: CkAbort("impossible iopt");break;
  }//end switch

  // you must zero because messages don't fill the plane
  if(countGradVks[iopt]==1){bzero(dataR,sizeof(double)*nptsExpnd);}

  // scaling here should be OK saving a few flops
  for(int i=0;i<size;i++){dataC[tranUnpack[Index][i]] = partiallyFFTd[i]*probScale;}

  delete msg;

//============================================================================
// When you have all the data : finish the FFT back to real space

  if (countGradVks[iopt] == nchareG){

    countGradVks[iopt]=0;
    doneGradRhoVks++;

    FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    // split this 2D XY to just Y

    fftcache->doRhoFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(fwFFTGtoRnot0_, StartTime, CmiWallTimer());    
#endif

    // transpose and communicate

    // fft along X

  }//endif : I captured a divRho

//============================================================================
// When you have rhoiRX,rhoiRY,rhoiRZ and Vks invoke gradient correction

  if(doneGradRhoVks==3){ 
    doneGradRhoVks=0;
    RTH_Runtime_resume(run_thread);
  }//endif

//----------------------------------------------------------------------------
   }//end routine
//============================================================================



//============================================================================
// The gradient of the density is now completed. You can compute the 
// GGA-DFT functional now. Density is now available to be used as scratch.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
   int nf2                    = ngridb;
   int nf3                    = ngridc;
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

  // revise NF1 NF2 NF3
    CPXCFNCTS::CP_getGGAFunctional(npts,nf1,nf2,nf3,density,rhoIRX,rhoIRY,rhoIRZ,
                                   Vks,thisIndex.x,exc_gga_ret);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(GradCorrGGA_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_VERSION_BLUEGENE
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

  // Candidate location for launch of NL Real, but only if you have grad corr on

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
// The white-bird term : First fwfft redefined delrho(r) to delrho(g)
// then send to RhoGspacePlane. RhoGspacePlane sends you back back another term.
// After this routine, rhoIRX, rhoIRY and rhoIRZ are `free'.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::whiteByrdFFT(){
//============================================================================

   FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
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
  // revise this 2D XY to just X
  fftcache->doRhoFFTRtoG_Rchare(rhoIRXC,rhoIRX,nplane_rho_x,ngrida,ngridb);


  // transpose communicate

  // FFT along Y
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(WhiteByrdFFTX_, StartTime, CmiWallTimer());    
#endif
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

  sendPartlyFFTtoRhoG(ioptx);

//============================================================================
// II) rhoIRY : Scale, real to complex FFT, perform FFT, transpose

#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  rho_rs.scale(rhoIRY,FFTscale);
  // revise this 2D XY to just X
  fftcache->doRhoFFTRtoG_Rchare(rhoIRYC,rhoIRY,nplane_rho_x,ngrida,ngridb);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(WhiteByrdFFTY_, StartTime, CmiWallTimer());    
#endif
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

  // transpose communicate

  // FFT along Y


  sendPartlyFFTtoRhoG(iopty);

//============================================================================
// III) rhoIRZ : Scale, real to complex FFT, perform FFT, transpose


#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  rho_rs.scale(rhoIRZ,FFTscale);
  // split this 2D XY to just X
  fftcache->doRhoFFTRtoG_Rchare(rhoIRZC,rhoIRZ,nplane_rho_x,ngrida,ngridb);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(WhiteByrdFFTZ_, StartTime, CmiWallTimer());    
#endif
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif
  // transpose communicate

  // FFT along Y

  sendPartlyFFTtoRhoG(ioptz);

//============================================================================
   }//end routine
//============================================================================



//============================================================================
// The white bird force correction : finish fft, add to vks
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptWhiteByrd(RhoRSFFTMsg *msg){
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("WhiteByrd Data from RhoG arriving at RhoR : %d %d\n",
	   thisIndex.x,thisIndex.y);
#endif

//============================================================================

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoG;
  int **tranUnpack       = sim->index_tran_upack_rho;
  int *nlines_per_chareG = sim->nlines_per_chareRhoG;
   
  int size               = msg->size; 
  int Index              = msg->senderIndex;
  complex *partiallyFFTd = msg->data;

//============================================================================
// Perform some error checking

  countWhiteByrd++;
  if (countWhiteByrd > nchareG) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               countWhiteByrd,nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=nlines_per_chareG[Index]){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.2, %d != %d for rho chare %d %d\n",size,nlines_per_chareG[Index],
                  thisIndex.y,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message : The density is free for use as scr

  double  *dataR = rho_rs.density; 
  complex *dataC = rho_rs.densityC; 

  // zero because input data does not fill the plane
  if(countWhiteByrd==1){bzero(dataR,sizeof(double)*nptsExpnd);}

  for(int i=0;i<size;i++){dataC[tranUnpack[Index][i]] = partiallyFFTd[i];}

  delete msg;

//============================================================================
// When you have all the messages, do the last fft, and add in the correction
// and resume

  if (countWhiteByrd == nchareG){

    countWhiteByrd=0;

    FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
    double *Vks        = rho_rs.Vks;
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    // Split this XY to just Y
    fftcache->doRhoFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb);
    
    // transpose and communicate

    // FFT along X

    for(int i=0;i<nptsExpnd;i++){Vks[i] -= dataR[i];}

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(PostByrdfwFFTGtoR_, StartTime, CmiWallTimer());    
#endif

#ifdef _CP_DEBUG_RHOR_VKSE_
    char myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "WaRho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
    FILE *fp = fopen(myFileName,"w");
    for (int i = 0; i <ngridb*ngrida; i++){
        fprintf(fp,"%g\n",rho_rs.Vks[i]);
    }//endfor
#endif

    doneWhiteByrd=true;
    if(doneHartVks){
      RTH_Runtime_resume(run_thread);
    }//endif

  }//endif : communication from rhog has all arrived safely

//============================================================================
   }//end routine
//============================================================================



//============================================================================
/**
 *
 * Accept hartExt tranpose data : receive VksHartEext(gx,gy,z) z is parallel.
 * Since the message can come at any time memory, VksHart has to be ready to
 * receive it.
 *
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptHartVks(RhoHartRSFFTMsg *msg){
//============================================================================
// Local pointers

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoG;
  int **tranUnpack       = sim->index_tran_upack_rho;
  int *nlines_per_chareG = sim->nlines_per_chareRhoG;
   
  int size               = msg->size; 
  int Index              = msg->senderBigIndex;
  int istrt_lines        = msg->senderStrtLine;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;

  CkAssert(iopt==0);
  countGradVks[iopt]++;

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("Data from RhoG arriving at RhoR : %d %d %d %d\n",
	   thisIndex.x,thisIndex.y,iopt,countGradVks[iopt]);
#endif

//============================================================================
// Copy out the hart-eext contrib to vks : iopt==0

  double  *dataR = rho_rs.VksHart;
  complex *dataC = rho_rs.VksHartC;

  // you must zero because partiallyFFT does not fill the plane
  if(countGradVks[iopt]==1){bzero(dataR,sizeof(double)*nptsExpnd);}

  for(int i=0,j=istrt_lines;i<size;i++,j++){
    dataC[tranUnpack[Index][j]] = partiallyFFTd[i];
  }//endfor

  delete msg;  

//============================================================================
// fft the puppy if you've got it all

  if (countGradVks[iopt] == nchareG*rhoGHelpers){
      countGradVks[iopt]=0;

      FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
      double *vks        = rho_rs.Vks;
#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif
      // split this 2D XY to just Y
      fftcache->doRhoFFTGtoR_Rchare(dataC,dataR,nplane_rho_x,ngrida,ngridb);

      // tranpose and communicate 

      // FFT along X

      for(int i=0;i<nptsExpnd;i++){vks[i]+=dataR[i];}

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(fwFFTGtoR0_, StartTime, CmiWallTimer());    
#endif
      doneHartVks = true;
      if(doneWhiteByrd){
	  RTH_Runtime_resume(run_thread);
      }//endif : check if whitebyrd is done

  // Candidate location for launch of NL Real

  }//endif : communication from rhog 


  
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

  if(!doneWhiteByrd || !doneHartVks|| !doneRhoReal || !doneRHart){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Flow of Control Error : Attempting to rho multicast\n");
    CkPrintf("without harteext or gradcorr (whitebyrd) \n");
    CkPrintf("without rhodensity or rhart \n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  countRHart=0;
  doneWhiteByrd    = false;
  doneHartVks      = false;
  doneRhoReal      = false;
  doneRHart        = false;
    
//============================================================================
// Send vks back to the states in real space

   int dataSize    = ngrida*ngridb;
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
      double *dataR    = msg->data;

      rho_rs.uPackShrink(dataR,Vks); // down pack Vks for the send
 
      if(config.useCommlibMulticast){mcastInstance.beginIteration();}
      if(config.useCommlibMulticast){
        realSpaceSectionCProxy.doProduct(msg);
      }else{
        realSpaceSectionProxy.doProduct(msg);
      }//enddif
     if(config.useCommlibMulticast){mcastInstance.endIteration();}

   }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::ResumeFromSync(){

    if(config.useCommlibMulticast)
	ComlibResetSectionProxy(&realSpaceSectionCProxy);
    if(config.useRInsRhoGP)
	ComlibResetProxy(&rhoGProxy_com);
}
//============================================================================

void CP_Rho_RealSpacePlane::RHartReport()
{
  countRHart++;
  if(countRHart==countRHartValue)
    {
      doneRHart=true;
      if(doneRhoReal)
	doMulticast();
      
    }
}


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//! return tru if input is power of 2
//============================================================================
bool is_pow2(int input){
    int y=0;
    for(int x=0;x<32;x++)
    {
	y= 1<<x;
	if(y==input)
	    return true;
    }
    return false;

//---------------------------------------------------------------------------
  }//end routine
//============================================================================



