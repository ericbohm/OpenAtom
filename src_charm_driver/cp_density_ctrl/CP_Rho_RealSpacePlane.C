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
#include "sim_subroutines.h"
#include "CP_State_Plane.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern CProxy_CP_State_RealSpacePlane realSpacePlaneProxy;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern ComlibInstanceHandle commRealInstance;
extern ComlibInstanceHandle commRealIGXInstance;
extern ComlibInstanceHandle commRealIGYInstance;
extern ComlibInstanceHandle commRealIGZInstance;
extern ComlibInstanceHandle mcastInstance;
extern CProxy_FFTcache fftCacheProxy;
extern CkGroupID mCastGrpId;
extern Config config;
extern int nstates;

bool is_pow2(int );

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
    RTH_Suspend(); 
    // 3rd entry point is acceptGradRhoVks(RhoRSFFTMsg *) 
    if(c->cp_grad_corr_on!=0){
       c->GradCorr();
       RTH_Suspend(); 
    }//endif
    // 4th entry point is acceptWhiteByrd(RhoRSFFTMsg *) || acceptHartVks.
    // Whichever arrives LAST calls resume (the last shall be first).
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
CP_Rho_RealSpacePlane::CP_Rho_RealSpacePlane(int xdim, size2d yzdim, 
				     int numRealSpace, int numRhoG, bool _useCommlib)
//============================================================================
   {//begin routine
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
    CkPrintf("[%d %d] RhoR constructs \n",thisIndex.x, thisIndex.y);
#endif

//============================================================================
// set flag

//============================================================================
// Malloc Memory, Intialize counters, set constants

    initRhoRealSlab(&rho_rs, xdim, yzdim[0], yzdim[1], numRealSpace, 
                    numRhoG,thisIndex.x,thisIndex.y);

    count          = 0;
    countWhiteByrd = 0;
    doneGradRhoVks = 0;
    doneHartVks    = true;
    doneWhiteByrd  = true;
    rhoGHelpers    = config.rhoGHelpers;
    for(int i=0;i<4;i++){countGradVks[i]=0;}

    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    FFTscale             = 1.0/((double)config.numFFTPoints);
    volumeFactor         = (sim->vol)*FFTscale;
    probScale            = (1.0/(sim->vol));

    cp_grad_corr_on      = sim->cp_grad_corr_on;

//============================================================================
// make sections in the realSpacePlane array. These will be used when 
// computing real-space densities and multicasting v_ks values 

    numMcastSent = 0;
    CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
    CkArrayIndexMax *elems = new CkArrayIndexMax[nstates];

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
    realSpaceSectionProxy.init(dummyProductMessage);

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

    usesAtSync = CmiTrue;
    //    if(config.lbdensity){
    //      setMigratable(true);
    //    }else{
      setMigratable(false);
      //    }//endif
//============================================================================

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

void CP_Rho_RealSpacePlane::pup(PUP::er &p)
{
  ArrayElement2D::pup(p);
  p|cp_grad_corr_on;
  p|FFTscale;        
  p|volumeFactor;        
  p|probScale;             
  p|count;
  p|countFFTdata;
  p|numMcastSent;
  p|rhoGHelpers;
  p|realSpaceSectionProxy;
  PUParray(p,countGradVks,5);
  p|doneGradRhoVks;
  p|countWhiteByrd;
  p|doneWhiteByrd;
  p|doneHartVks;
  rho_rs.pup(p); 
  p|realSpaceSectionCProxy;
  if(p.isUnpacking())
    {
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

    }

  /* 
  if(p.isUnpacking())
    {
      run_thread = RTH_Runtime_create(RTH_Routine_lookup(CP_Rho_RealSpacePlane,run),this);
      RTH_Runtime_resume(run_thread);
    }
  RTH_Runtime_pup(run_thread,p,this);
  */

}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * Here the density from all the states is added up. The data from all the
 * states is received via an array section reduction.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity(CkReductionMsg *msg) {
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
    CkPrintf("RhoReal accepting Density %d %d %d\n",
              thisIndex.x,thisIndex.y,CkMyPe());
#endif
   if(cp_grad_corr_on!=0){
     doneWhiteByrd  = false;
   }else{
     doneWhiteByrd  = true;
   }
    doneHartVks    = false;
#ifdef _CP_DEBUG_HARTEEXT_OFF_
    doneHartVks    = true;
#endif

//============================================================================

    CkAssert(msg->getSize() == rho_rs.sizeZ * rho_rs.sizeX * sizeof(double));
    double *realValues = (double *) msg->getData(); 
    double *density    = rho_rs.density;
    for(int i = 0; i < rho_rs.sizeZ*rho_rs.sizeX; i++){
      density[i] = realValues[i] * probScale;	
    }//endfor
    delete msg;

//============================================================================

#ifdef _CP_DEBUG_RHOR_RHO_
      char myFileName[MAX_CHAR_ARRAY_LENGTH];
      sprintf(myFileName, "Rho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
      FILE *fp = fopen(myFileName,"w");
        for (int i = 0; i <rho_rs.sizeZ*rho_rs.sizeX; i++){
          fprintf(fp,"%g\n",realValues[i]*probScale);
        }//endfor
      fclose(fp);
#endif

//============================================================================

    RTH_Runtime_resume(run_thread);

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity() {
	energyComputation();
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** 
 * Make a copy of densities and compute one part of the energies using
 * CP_exc_calc
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
   int npts         = rho_rs.trueSize;
   int nf1          = rho_rs.sizeX;
   int nf2          = rho_rs.sizeY;
   int nf3          = rho_rs.sizeZ;
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
 * do FFT on rho_rs.doFFTonThis
 */
//============================================================================
void CP_Rho_RealSpacePlane::fftRhoRtoRhoG(){
//============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("In RhoRealSpacePlane[%d %d] FFT_RSpacetoGSpace %d %d\n",thisIndex.x, 
                   thisIndex.y, CkMyPe(),CmiMemoryUsage());
#endif

//============================================================================

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
 
  double *density     = rho_rs.density;
  double *doFFTonThis = rho_rs.doFFTonThis;
  rho_rs.uPackAndScale(doFFTonThis,density,volumeFactor);
  fftCacheProxy.ckLocalBranch()->doRhoRealtoRhoG(doFFTonThis);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(RhoRtoGFFT_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

//============================================================================
// Send chunk to RhoGDensity

 int iopt = 0;
 sendPartlyFFTtoRhoG(iopt);

//============================================================================
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::sendPartlyFFTtoRhoG(int iopt){
//============================================================================

    CPcharmParaInfo *sim      = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    int nchareRhoG            = sim->nchareRhoG;
    int **tranpack_rho        = sim->index_tran_upack_rho;
    int *nlines_per_chareRhoG = sim->nlines_per_chareRhoG;

    double *FFTresultR;
    switch(iopt){   
      case 0: FFTresultR = rho_rs.doFFTonThis; break;
      case 1: FFTresultR = rho_rs.rhoIRX;      break;
      case 2: FFTresultR = rho_rs.rhoIRY;      break;
      case 3: FFTresultR = rho_rs.rhoIRZ;      break;
    }//end switch
    complex *FFTresult = reinterpret_cast<complex*> (FFTresultR);

//============================================================================
// Launch the communication

	switch(iopt){
	    case 0 : if(config.useRInsRhoGP) commRealInstance.beginIteration(); break;
	    case 1 : if(config.useRInsIGXRhoGP) commRealIGXInstance.beginIteration(); break;
	    case 2 : if(config.useRInsIGYRhoGP) commRealIGYInstance.beginIteration(); break;
	    case 3 : if(config.useRInsIGZRhoGP) commRealIGZInstance.beginIteration(); break;

	}
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
        case 0 : rhoGProxy_com(ic,0).acceptData(msg);   break;
        case 1 : rhoGProxyIGX_com(ic,0).acceptWhiteByrd(msg); break;
        case 2 : rhoGProxyIGY_com(ic,0).acceptWhiteByrd(msg); break;
        case 3 : rhoGProxyIGZ_com(ic,0).acceptWhiteByrd(msg); break;

      }//end switch

    }//end for : chare sending

	switch(iopt){
	    case 0 : if(config.useRInsRhoGP) commRealInstance.endIteration(); break;
	    case 1 : if(config.useRInsIGXRhoGP) commRealIGXInstance.endIteration(); break;
	    case 2 : if(config.useRInsIGYRhoGP) commRealIGYInstance.endIteration(); break;
	    case 3 : if(config.useRInsIGZRhoGP) commRealIGZInstance.endIteration(); break;
	}

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
/**
 *accept tranpose data : receive grad_rho(z,gy,gx) z is parallel
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================

void CP_Rho_RealSpacePlane::acceptGradRhoVks(RhoRSFFTMsg *msg){

//============================================================================

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoG;
  int **tranUnpack       = sim->index_tran_upack_rho;
  int *nlines_per_chareG = sim->nlines_per_chareRhoG;
   
  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;
  int pSize              = (rho_rs.sizeX+2)*(rho_rs.sizeY);
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
               count,nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=nlines_per_chareG[Index]){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude, %d != %d for rho chare %d %d\n",size,nlines_per_chareG[Index],
                  thisIndex.y,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message

  double *data;
  switch(iopt){
      case 0: data = rho_rs.doFFTonThis; break;
      case 1: data = rho_rs.rhoIRX; break;
      case 2: data = rho_rs.rhoIRY; break;
      case 3: data = rho_rs.rhoIRZ; break;
  }//end switch
  complex *planeArr = reinterpret_cast<complex*> (data);
  if(countGradVks[iopt]==1){memset(data,0,sizeof(double)*pSize);}

  for(int i=0;i<size;i++){planeArr[tranUnpack[Index][i]] = partiallyFFTd[i];}
  delete msg;

//============================================================================
// When you have all the data : finish the FFT back to real space

  if (countGradVks[iopt] == nchareG){
    countGradVks[iopt]=0;
    doneGradRhoVks++;

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

    if(iopt!=0){rho_rs.doFwFFTGtoR(iopt,probScale);}

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(fwFFTGtoRnot0_, StartTime, CmiWallTimer());    
#endif
  }

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
   double *rhoIRX             = rho_rs.rhoIRX;
   double *rhoIRY             = rho_rs.rhoIRY;
   double *rhoIRZ             = rho_rs.rhoIRZ;
   double *Vks                = rho_rs.Vks;
   int size                   = rho_rs.size;
   int npts                   = rho_rs.trueSize;
   int nf1                    = rho_rs.sizeX;
   int nf2                    = rho_rs.sizeY;
   int nf3                    = rho_rs.sizeZ;
   double *gradientCorrection = rho_rs.gradientCorrection;
   double *exc_gga_ret        = &(rho_rs.exc_gga_ret);

//============================================================================

#ifdef _CP_DEBUG_RHOR_VKSC_
    char myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "BGradRho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
    FILE *fp = fopen(myFileName,"w");
      for (int i = 0; i <rho_rs.sizeY*rho_rs.sizeX; i++){
        fprintf(fp,"%g %g %g %g\n",rho_rs.rhoIRX[i],rho_rs.rhoIRY[i],rho_rs.rhoIRZ[i],
                                   rho_rs.Vks[i]);
      }//endfor
    fclose(fp);
#endif

//============================================================================
// Compute the gradient corrected functional

    rho_rs.exc_gga_ret = 0.0;
    for(int i=0;i<size;i++){gradientCorrection[i] = 0.0;}
#define GGA_ON
#ifdef GGA_ON

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif

    CPXCFNCTS::CP_getGGAFunctional(npts,nf1,nf2,nf3,density,
               rhoIRX,rhoIRY,rhoIRZ,gradientCorrection,thisIndex.x,exc_gga_ret);
    for(int i=0;i<npts;i++){Vks[i] += gradientCorrection[i];}

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(GradCorrGGA_, StartTime, CmiWallTimer());    
#endif

#endif

#ifdef CMK_VERSION_BLUEGENE
     CmiNetworkProgress();
#endif

//============================================================================
// Reduce the exchange correlation energy

   double exc[2];
   exc[0]=rho_rs.exc_ret;
   exc[1]=rho_rs.exc_gga_ret;
   contribute(2*sizeof(double),exc,CkReduction::sum_double);

//============================================================================

#ifdef _CP_DEBUG_RHOR_VKSD_
    myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "AGradRho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
    fp = fopen(myFileName,"w");
      for (int i = 0; i <rho_rs.sizeY*rho_rs.sizeX; i++){
        fprintf(fp,"%g %g %g %g\n",rho_rs.rhoIRX[i],rho_rs.rhoIRY[i],rho_rs.rhoIRZ[i],
                                   rho_rs.Vks[i]);
      }//endfor
    fclose(fp);
#endif

//============================================================================
// Start the white bird puppy : back fft of rhoirx, rhoiry, rhoirz

    whiteByrdFFT();

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// The white-bird term : First fwfft of redined delrho(r) to delrho(g)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::whiteByrdFFT(){
//============================================================================

   int ioptx           = 1;   
   int iopty           = 2;
   int ioptz           = 3;
   int npts            = rho_rs.trueSize;
   double *rhoIRX      = rho_rs.rhoIRX;
   double *rhoIRY      = rho_rs.rhoIRY;
   double *rhoIRZ      = rho_rs.rhoIRZ;
   double *gradientCorrection = rho_rs.gradientCorrection; //scratch

//============================================================================
// I) rhoIRX : Unpack for real to complex FFT, perform FFT, transpose

  memcpy(gradientCorrection,rhoIRX,npts*sizeof(double));
  rho_rs.uPackAndScale(rhoIRX,gradientCorrection,FFTscale);
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
  fftCacheProxy.ckLocalBranch()->doRhoRealtoRhoG(rhoIRX);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(WhiteByrdFFTX_, StartTime, CmiWallTimer());    
#endif
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

  sendPartlyFFTtoRhoG(ioptx);

//============================================================================
// II) rhoIRY : Unpack for real to complex FFT, perform FFT, transpose

  memcpy(gradientCorrection,rhoIRY,npts*sizeof(double));
  rho_rs.uPackAndScale(rhoIRY,gradientCorrection,FFTscale);
#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif
  fftCacheProxy.ckLocalBranch()->doRhoRealtoRhoG(rhoIRY);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(WhiteByrdFFTY_, StartTime, CmiWallTimer());    
#endif
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

  sendPartlyFFTtoRhoG(iopty);

//============================================================================
// III) rhoIRZ : Unpack for real to complex FFT, perform FFT, transpose


  memcpy(gradientCorrection,rhoIRZ,npts*sizeof(double));
  rho_rs.uPackAndScale(rhoIRZ,gradientCorrection,FFTscale);
#ifndef CMK_OPTIMIZE
  StartTime=CmiWallTimer();
#endif

  fftCacheProxy.ckLocalBranch()->doRhoRealtoRhoG(rhoIRZ);

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(WhiteByrdFFTZ_, StartTime, CmiWallTimer());    
#endif

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

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
  int pSize              = (rho_rs.sizeX+2)*(rho_rs.sizeY);

//============================================================================
// Perform some error checking

  countWhiteByrd++;
  if (countWhiteByrd > nchareG) {
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Mismatch in allowed rho_gspace chare arrays : %d %d %d %d\n",
               count,nchareG,thisIndex.x,thisIndex.y);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  if(size!=nlines_per_chareG[Index]){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude, %d != %d for rho chare %d %d\n",size,nlines_per_chareG[Index],
                  thisIndex.y,Index);
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//============================================================================
// unpack the data and delete the message

  double *data      = rho_rs.rhoIRX; // note this is just a memory reuse trick
                                     // since we are long done with rhoIRX by
                                     // this point
  complex *planeArr = reinterpret_cast<complex*> (data);
  if(countWhiteByrd==1){memset(data,0,sizeof(double)*pSize);}

  for(int i=0;i<size;i++){planeArr[tranUnpack[Index][i]] = partiallyFFTd[i];}
  delete msg;

//============================================================================
// When you have all the messages, do the last fft, and add in the correction
// and resume

  if (countWhiteByrd == nchareG){
    countWhiteByrd=0;
    int iopt = 1;// reusing rhoIRX
    double scale = 1.0;

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif

    rho_rs.doFwFFTGtoR(iopt,scale);

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(PostByrdfwFFTGtoR_, StartTime, CmiWallTimer());    
#endif

    double *Vks       = rho_rs.Vks;
    double *whitebyrd = data;
    int npts          = (rho_rs.sizeX)*(rho_rs.sizeY);
    for(int i=0;i<npts;i++){Vks[i] -= whitebyrd[i];}

#ifdef _CP_DEBUG_RHOR_VKSE_
    char myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "WaRho_Real_%d_%d.out", thisIndex.x,thisIndex.y);
    FILE *fp = fopen(myFileName,"w");
    for (int i = 0; i <rho_rs.sizeY*rho_rs.sizeX; i++){
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
 *accept hartExt tranpose data : receive grad_rho(z,gy,gx) z is parallel
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::acceptHartVks(RhoHartRSFFTMsg *msg){
//============================================================================

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoG;
  int **tranUnpack       = sim->index_tran_upack_rho;
  int *nlines_per_chareG = sim->nlines_per_chareRhoG;
   
  int size               = msg->size; 
  int Index              = msg->senderBigIndex;
  int istrt_lines        = msg->senderStrtLine;
  int iopt               = msg->iopt;
  complex *partiallyFFTd = msg->data;
  int pSize              = (rho_rs.sizeX+2)*(rho_rs.sizeY);
#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("Data from RhoG arriving at RhoR : %d %d %d %d\n",
	   thisIndex.x,thisIndex.y,iopt,countGradVks[iopt]);
#endif
  countGradVks[iopt]++;
  double *data      = rho_rs.doFFTonThis;
  complex *planeArr = reinterpret_cast<complex*> (data);
  if(countGradVks[iopt]==1){memset(data,0,sizeof(double)*pSize);}

  for(int i=0,j=istrt_lines;i<size;i++,j++){
    planeArr[tranUnpack[Index][j]] = partiallyFFTd[i];
  }//endfor

  delete msg;  
  CkAssert(iopt==0);
  if (countGradVks[iopt] == nchareG*rhoGHelpers){
      countGradVks[iopt]=0;
      double scale = 1.0;
#ifdef _CP_DEBUG_RHOR_VKSA_
      char fmyFileName[MAX_CHAR_ARRAY_LENGTH];
      sprintf(fmyFileName, "HartRho_Realb4fft_%d_%d_0.out", thisIndex.x,thisIndex.y);
      FILE *ffp = fopen(fmyFileName,"w");
      for(int i=0;i<pSize;i++){
	fprintf(ffp,"%g\n",data[i]);
      }//endfor
      fclose(ffp);
#endif

#ifndef CMK_OPTIMIZE
      double StartTime=CmiWallTimer();
#endif

      rho_rs.doFwFFTGtoR(4,scale);

#ifndef CMK_OPTIMIZE
      traceUserBracketEvent(fwFFTGtoR0_, StartTime, CmiWallTimer());    
#endif
      double *vksExc     = rho_rs.Vks;
      double *vksHartExt = data;
      int size           = (rho_rs.sizeX)*(rho_rs.sizeY);
#ifdef _CP_DEBUG_RHOR_VKSA_
      char myFileName[MAX_CHAR_ARRAY_LENGTH];
      sprintf(myFileName, "VksRho_Real_%d_%d_0.out", thisIndex.x,thisIndex.y);
      FILE *fp = fopen(myFileName,"w");
      for(int i=0;i<size;i++){
	fprintf(fp,"%g %g\n",vksExc[i],vksHartExt[i]);
      }//endfor
      fclose(fp);
#endif
      for(int i=0;i<size;i++){vksExc[i]+=vksHartExt[i];}

      doneHartVks = true;
      if(doneWhiteByrd){
	  RTH_Runtime_resume(run_thread);
      }//endif : check if whitebyrd is done

  }//endif : communication from rhog 

//============================================================================
  }//end routine 
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::doMulticast(){
//============================================================================
// Check for flow of control errors

  if(!doneWhiteByrd || !doneHartVks){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Flow of Control Error : Attempting to rho multicast\n");
    CkPrintf("without harteext or gradcorr (whitebyrd) \n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// Send vks back to the states in real space
   if ((config.useGMulticast+config.useCommlibMulticast)!=1) {
     CkAbort("No multicast strategy\n");
   }//endif

   if (config.useGMulticast || config.useCommlibMulticast) {

      int dataSize    = rho_rs.sizeX * rho_rs.sizeZ;
      ProductMsg *msg = new (dataSize, 0) ProductMsg;
      msg->idx        = thisIndex.x;
      msg->datalen    = dataSize;
      msg->hops       = 0;
      double *data    = msg->data;
      double *Vks     = rho_rs.Vks;
      for(int i=0;i<dataSize;i++){data[i] = Vks[i];}
 
      if(config.useCommlibMulticast){mcastInstance.beginIteration();}
      if(config.useCommlibMulticast){
        realSpaceSectionCProxy.doProduct(msg);
      }else{
        realSpaceSectionProxy.doProduct(msg);
      }//enddif
     if(config.useCommlibMulticast){mcastInstance.endIteration();}

   }//endif
   bzero(rho_rs.Vks,rho_rs.sizeX * rho_rs.sizeZ*sizeof(double));

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



