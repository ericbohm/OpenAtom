//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_RealSpacePlane.C
 *
 * @defgroup Density Density
 *
 * \brief Computes electron density in real space, (exchange correlation energy) for transforming to
 * CP_Rho_GSpacePlane which will utilize CP_Rho_GHartExt and
 * CP_Rho_RHartExt (Hartree and external energies). Each plane
 * may be further subdivided into subplanes at runtime for additional
 * parallelism.
 *
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
#include <iostream>
#include <fstream>
#include <cmath>

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "cp_state_ctrl/CP_State_Plane.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "fft_charm.h"

//============================================================================
extern CProxy_TimeKeeper                                TimeKeeperProxy;
extern CkVec <CProxy_CP_State_RealSpacePlane>           UrealSpacePlaneProxy;
extern CProxy_HFCalculator HFCalculatorProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane>        UrealParticlePlaneProxy;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>             UrhoRealProxy;
extern CkVec <CProxy_CP_Rho_GSpacePlane>                UrhoGProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>                   UrhoRHartExtProxy;
extern CkVec <CProxy_GSpaceDriver>                      UgSpaceDriverProxy;

//ids for FFTs
extern CkVec <CProxy_fft2d>                     Urho_fft_xProxy, Urho_fft_yProxy,
                                                Urho_fft_zProxy, Urho_fft_hartProxy;
extern CkGroupID                                mCastGrpId;
extern CPcharmParaInfo                          simReadOnly;
extern Config                                   config;

//! return true if input is power of 2
bool is_pow2(int input){
  unsigned x = input;
  return x && !(x & (x - 1));
#if 0 //old code, why so much pain?
  int y = 0;
  for(int x = 0; x < 32;x++){
    y = 1 << x;
    if(y == input){
      return true;
    }
  }//endfor
  return false;
#endif
}

//#define  _CP_DEBUG_RHOR_VERBOSE_ 1

/** @addtogroup Density
  @{
 */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the states
// Performs lots of operations to get exc, eext energies and vks
//
//============================================================================
CP_Rho_RealSpacePlane::CP_Rho_RealSpacePlane(int _rhokeeperid,
    UberCollection _instance) : rhoKeeperId(_rhokeeperid), thisInstance(_instance)
{

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d %d] RhoR constructs \n",thisIndex.x, thisIndex.y);
#endif

  FFTscale     = 1.0/((double)config.numFFTPoints);  // these are based on the full size
  volumeFactor = simReadOnly.vol * FFTscale;
  probScale    = 1.0 / simReadOnly.vol;

  //NULL initialize all memory
  exc_ret=muxc_ret = exc_gga_ret=0.0;
  Vks       = NULL;
  VksC      = NULL;
  density   = NULL;
  densityC  = NULL;
  rhoIRX    = NULL;
  rhoIRXC   = NULL;
  rhoIRY    = NULL;
  rhoIRYC   = NULL;
  rhoIRZ    = NULL;
  rhoIRZC   = NULL;
  Vks       = NULL;
  VksHartC  = NULL;

  //============================================================================
  // Initialize counters, set booleans.myTime
  myTime          = 0;
  doneGradRhoVks  = 0;
  countRHart      = 0;

  countRHartValue = 1;

  doneHartVks     = false;
  doneRHart       = false;
  doneWhiteByrd   = false;

  //============================================================================
  // Migration

  usesAtSync = true;
  setMigratable(false);
}//end routine
//============================================================================

/** post constructor initialization */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RealSpacePlane::init(){

  //query fft library to find my extents
  Charm_getInputIndex(thisIndex.x, thisIndex.y, fft_xoffset,
      Urho_fft_xProxy[thisInstance.proxyOffset]);
  Charm_getInputIndex(thisIndex.x, thisIndex.y, fft_yoffset,
      Urho_fft_yProxy[thisInstance.proxyOffset]);
  Charm_getInputIndex(thisIndex.x, thisIndex.y, fft_zoffset,
      Urho_fft_zProxy[thisInstance.proxyOffset]);
  Charm_getInputIndex(thisIndex.x, thisIndex.y, fft_hartoffset,
      Urho_fft_hartProxy[thisInstance.proxyOffset]);
  Charm_getInputExtents(myGrid_start[MY_X], myGrid_end[MY_X],
                        myGrid_start[MY_Y], myGrid_end[MY_Y],
                        myGrid_start[MY_Z], myGrid_end[MY_Z],
                        Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);
  for(int i = 0; i < 3; i++) {
    myGrid_length[i] = myGrid_end[i] - myGrid_start[i];
  }

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d %d] RhoR init : grid size %d %d %d \n",thisIndex.x, thisIndex.y,
    myGrid_length[0], myGrid_length[1], myGrid_length[2]);
#endif

  /*Don't be confused: get sublime.
    Science grid coordinates were named "a, b, c" where a is the fastest
    growing dimension (thanks Fortran).
    FFT grid named its coordinates as "x, y, z" with z growing fastest, i.e. it
    does z FFT first.
    To match, we had to then map them accordingly a->z, b->y, c->x.
    This mapping is different from the mapping science team assumes; they
    map a->x, b->y, c->z (where x,y,z is used for God knows what)
   */

  //many reductions are targeted at me (atleast for now)
  redCount = new int[myGrid_length[MY_C]];
  RedMsg = new CkReductionMsg*[myGrid_length[MY_C]];
  num_redn_complete = 0;
  for(int i = 0; i < myGrid_length[MY_C]; i++) {
    redCount[i] = 0;
    RedMsg[i] = NULL;
  }

  //============================================================================
  // Set up the data class : mallocs rho,gradrho, etc.
  // The following replaces the Rho slab and related code.

  int data_size = myGrid_length[MY_C] * myGrid_length[MY_B] *
                  (myGrid_length[MY_A]/2 + 1); //FFTW style storage
  myGrid_size = 2 * data_size; //how many real points I have

  complex *dummy;
  dummy    = (complex*) fftw_malloc(data_size * sizeof(complex));
  VksC     = dummy;
  Vks      = reinterpret_cast<double*> (dummy);

  dummy    = (complex*) fftw_malloc(data_size * sizeof(complex));
  densityC = dummy;
  density  = reinterpret_cast<double*> (dummy);

  dummy    = (complex*) fftw_malloc(data_size * sizeof(complex));
  rhoIRXC  = dummy;
  rhoIRX   = reinterpret_cast<double*> (dummy);
  Charm_setInputMemory((void*)rhoIRX, Urho_fft_xProxy[thisInstance.proxyOffset],
      fft_xoffset);
  Charm_createInputPlan(Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);

  dummy    = (complex*) fftw_malloc(data_size * sizeof(complex));
  rhoIRYC  = dummy;
  rhoIRY   = reinterpret_cast<double*> (dummy);
  Charm_setInputMemory((void*)rhoIRY, Urho_fft_yProxy[thisInstance.proxyOffset],
      fft_yoffset);
  Charm_createInputPlan(Urho_fft_yProxy[thisInstance.proxyOffset], fft_yoffset);

  dummy    = (complex*) fftw_malloc(data_size * sizeof(complex));
  rhoIRZC  = dummy;
  rhoIRZ   = reinterpret_cast<double*> (dummy);
  Charm_setInputMemory((void*)rhoIRZ, Urho_fft_zProxy[thisInstance.proxyOffset],
      fft_zoffset);
  Charm_createInputPlan(Urho_fft_zProxy[thisInstance.proxyOffset], fft_zoffset);

  dummy    = (complex*) fftw_malloc(data_size * sizeof(complex));
  VksHartC = dummy;
  VksHart  = reinterpret_cast<double*> (dummy);
  Charm_setInputMemory((void*)VksHart,
      Urho_fft_hartProxy[thisInstance.proxyOffset], fft_hartoffset);
  Charm_createInputPlan(Urho_fft_hartProxy[thisInstance.proxyOffset],
      fft_hartoffset);

  //make sections in the realSpacePlane array. These will be used when
  //computing real-space densities and multicasting vks values
  CkMulticastMgr *mcastGrp = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
  CkArrayIndexMax **elems   = new
    CkArrayIndexMax*[myGrid_length[MY_C]];

  for(int plane = 0; plane < myGrid_length[MY_C]; plane++) {
    elems[plane] = new CkArrayIndexMax[simReadOnly.nstates];
  }
  // section i has all the portions with all
  for(int c_plane = myGrid_start[MY_C]; c_plane < myGrid_end[MY_C]; c_plane++) {
    CkArrayIndex2D idx(0, c_plane);
    if(is_pow2(simReadOnly.nstates)) {
      for (int j = 0; j < simReadOnly.nstates; j++) {
        idx.index[0] = j ^ ((thisIndex.x + thisIndex.y) % simReadOnly.nstates);
        elems[c_plane - myGrid_start[MY_C]][j] = idx;
      }//endfor
    } else {
      for (int j = 0; j < simReadOnly.nstates; j++) {
        idx.index[0] = (j + thisIndex.x + thisIndex.y) % simReadOnly.nstates;
        elems[c_plane - myGrid_start[MY_C]][j] = idx;
      }//endfor
    }//endif
  }

  // we need one RS proxy for each K-point until the cross proxies with
  // reductions work correctly
  UberCollection RhoReductionSource(thisInstance);
  realSpaceSectionProxyA= new CProxySection_CP_State_RealSpacePlane*[config.UberJmax];
#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d %d] RhoR create %d x %d mcast groups \n", thisIndex.x,
      thisIndex.y, config.UberJmax, myGrid_length[MY_C]);
#endif
  for(int kp = 0; kp < config.UberJmax; kp++)
  {
    realSpaceSectionProxyA[kp] = new CProxySection_CP_State_RealSpacePlane[myGrid_length[MY_C]];
    for(int c_plane = 0; c_plane < myGrid_length[MY_C]; c_plane++) {
      RhoReductionSource.idxU.y = kp; // not at the gamma point
      RhoReductionSource.setPO();
      realSpaceSectionProxyA[kp][c_plane] = CProxySection_CP_State_RealSpacePlane::
        ckNew(UrealSpacePlaneProxy[RhoReductionSource.proxyOffset].ckGetArrayID(),
        elems[c_plane], simReadOnly.nstates);

      realSpaceSectionProxyA[kp][c_plane].ckDelegate
        (CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());

      mcastGrp->setSection(realSpaceSectionProxyA[kp][c_plane]);

      InitDensity *indexMsg = new InitDensity;
      // inform realspace element of this section proxy.
      indexMsg->grid_offset_b = myGrid_start[MY_B];
      indexMsg->grid_num_b = myGrid_length[MY_B];
      indexMsg->pencil_offset_x = thisIndex.x;
      indexMsg->pencil_offset_y = thisIndex.y;
      realSpaceSectionProxyA[kp][c_plane].init(indexMsg);
#ifdef _CP_DEBUG_RHOR_VERBOSE_
      //CkPrintf("[%d %d] RhoR init : section %d %d \n",thisIndex.x, thisIndex.y,
      //    myGrid_start[MY_B], myGrid_length[MY_B]);
#endif
    }
  }

  delete [] elems;
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_RealSpacePlane::~CP_Rho_RealSpacePlane(){
  fftw_free((void*)VksC);
  fftw_free((void*)densityC);
  fftw_free((void*)rhoIRXC);
  fftw_free((void*)rhoIRYC);
  fftw_free((void*)rhoIRZC);
  fftw_free((void*)VksHartC);
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Data comes from StateRspacePlane once an algorithm step.
/**
 * Here the density from all the states is added up. The data from all the
 * states is received via possibly multiple array section reductions. Nothing
 * happens in this chare until the density arrives.
 *
 * 1) We will receive one reduction for every "c" plane part of our sub-grid.
 *
 * 2) If we're not at the gamma point, there will be UberJmax of 1).
 * Otherwise there will be only 1.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptDensity(CkReductionMsg *msg) {
  //============================================================================

  int plane_c = (int)msg->getUserFlag();
  plane_c -= myGrid_start[MY_C];

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d] RhoReal accepting Density %d %d for plane %d \n",
      CkMyPe(), thisIndex.x, thisIndex.y, plane_c);
#endif

#ifdef _NAN_CHECK_
  for(int i=0; i < msg->getSize()/sizeof(double); i++){
    CkAssert(isnan(((double*) msg->getData())[i]) == 0);
  }//endif
#endif

  redCount[plane_c]++;
  if(redCount[plane_c] == 1) {
    RedMsg[plane_c] = msg; // This guy is deleted next time we get a redMsg for the same plane
  } else {
    if(redCount[plane_c] <= config.UberJmax) {
      // there is one per k-point we'll sum them into the first
      // message we receive
      double *indata = (double*) msg->getData();
      double *outdata = (double*) RedMsg[plane_c]->getData();
      int size = msg->getSize()/sizeof(double);
      for(int i = 0; i < size ;i++)
        outdata[i] += indata[i];
      delete(msg);
    } else {
      CkAbort("Too many reduction messages received in density\n");
    }
  }//endif

  //TODO: scale plane_c when it has arrived; why wait for all planes?
  if(redCount[plane_c] == config.UberJmax){
    num_redn_complete++;
    redCount[plane_c] = 0;
    if(num_redn_complete == myGrid_length[MY_C]) {
      num_redn_complete = 0;
      handleDensityReduction();
    }
  }//endif
  //--------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
/**
 * Density has arrived. Set flags, scale, and trigger work.
 **/
//============================================================================
void CP_Rho_RealSpacePlane::handleDensityReduction() {
  //============================================================================

#ifdef _CP_SUBSTEP_TIMING_
  if(rhoKeeperId > 0){
    double rhostart = CmiWallTimer();
    CkCallback cb(CkIndex_TimeKeeper::collectStart(NULL), 0, TimeKeeperProxy);
    contribute(sizeof(double), &rhostart, CkReduction::min_double, cb, rhoKeeperId);
  }//endif
#endif

  //============================================================================
  // Set the flags : you are not done unless certain conditions apply.

  myTime++;
  doneHartVks   = false;
  doneWhiteByrd = false;
  doneRHart     = false;
  if(simReadOnly.cp_grad_corr_on == 0) { doneWhiteByrd = true; }
  if(simReadOnly.ees_eext_on == 0)     { doneRHart     = true; }

#ifdef _CP_DEBUG_HARTEEXT_OFF_
  doneHartVks    = true;
  doneRHart      = true;
#endif

  //============================================================================
  // Unpack into spread out form and delete the message

  double *scale_density = density;
  scaleData(scale_density, probScale);

  //============================================================================
  // If debugging, generate output!

#ifdef _CP_DEBUG_RHOR_RHO_
  char myFileName[MAX_CHAR_ARRAY_LENGTH];
  sprintf(myFileName, "Rho_Real_%d_%d.out", thisIndex.x, thisIndex.y);
  FILE *fp = fopen(myFileName,"w");
  scale_density = density;
  for (int i = 0; i < myGrid_size; i++){
    fprintf(fp, "%g\n", scale_density[i]);
  }//endfor
  fclose(fp);
#endif

  //============================================================================
  // Compute the exchange correlation energy (density no-grad part)

  energyComputation();

  //============================================================================
  // 2nd Launch real-space external-hartree and the G-space non-local
  // The energy comp through RhoG is the more expensive critical path.
  launchEextRNlG();

  //============================================================================
}//end routine
//============================================================================

void CP_Rho_RealSpacePlane::scaleData(double *scaledData, double scaleFac) {
  //read comment in inner loop for this end value
  int i, end = myGrid_length[MY_A] + (2 - (myGrid_length[MY_A] & 1));
  for(int plane = 0; plane < myGrid_length[MY_C]; plane++) {
    double *realValues = (double *) RedMsg[plane]->getData();
    for(int line = 0; line < myGrid_length[MY_B]; line++) {
      for(i = 0; i < myGrid_length[MY_A]; i++) {
        scaledData[i] = realValues[i] * scaleFac;
      }
      //for R2C FFT, we allocate (dim_len/2 + 1) complex numbers for each
      //of the grid line. This is required because R2C FFT generates those many
      //complex numbers out of a dimension dim_len. */
      for(; i < end; i++) {
        scaledData[i] = 0.0;
      }//endfor
      scaledData += 2*(myGrid_length[MY_A]/2 + 1);
      realValues += myGrid_length[MY_A];
    }//endfor
    delete RedMsg[plane];
  }
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
/** Launch ees NL and ees Ext routines
 *  Do this once an algorithm step
 */
//============================================================================
void CP_Rho_RealSpacePlane::launchEextRNlG() {
  //============================================================================
  // Launch the external energy computation in r-space :
#ifndef _CP_DEBUG_HARTEEXT_OFF_
  if(simReadOnly.ees_eext_on == 1){
    if(thisIndex.x + thisIndex.y == 0) {
      UrhoRHartExtProxy[thisInstance.proxyOffset].startEextIter();
    }
#if 0 //old code
    int div    = (ngridcEext / ngridc);
    int rem    = (ngridcEext % ngridc);
    int ind    = thisIndex.x+ngridc;
    for(int j=0;j<config.nchareHartAtmT;j++){
      UrhoRHartExtProxy[thisInstance.proxyOffset](thisIndex.x,thisIndex.y,j).startEextIter();
      if(thisIndex.x<rem){
        UrhoRHartExtProxy[thisInstance.proxyOffset](ind,thisIndex.y,j).startEextIter();
      }//endif : the launch
    }//endfor : atmTyp parallism
#endif
  }//endif : Launch is needed
#endif

  //============================================================================
  // Launch nonlocal g space if it wasn't done in RS
  //  Spread the launch over all the rhoRchares you can.

  if(simReadOnly.ees_nloc_on == 1 && config.launchNLeesFromRho == 1 && simReadOnly.natm_nl > 0) {
    /* Using a broadcast */
    if(thisIndex.x + thisIndex.y == 0) {
      for(int kpoint = 0; kpoint < config.UberJmax; kpoint++) {
        UberCollection destKpointInstance = thisInstance;
        destKpointInstance.idxU.y = kpoint;
        int proxyOffset = destKpointInstance.setPO();
        UgSpaceDriverProxy[proxyOffset].startNonLocalEes(myTime);
      }
    }
#if 0 //old code
    CkAssert(sizeZ>=config.nchareG);
    if(thisIndex.x<config.nchareG){
      int nstates = config.nstates;
      int div     = (nstates/rhoRsubplanes);
      int rem     = (nstates % rhoRsubplanes);
      int add     = (thisIndex.y < rem ? 1 : 0);
      int max     = (thisIndex.y < rem ? thisIndex.y : rem);
      int ist     = div*thisIndex.y + max;
      int iend    = ist + div + add;
      for(int kpoint=0;kpoint < config.UberJmax;kpoint++)
      {
        UberCollection destKpointInstance=thisInstance;
        destKpointInstance.idxU.y=kpoint;
        int proxyOffset=destKpointInstance.setPO();
        //TODO Change this to a section multicast (see Init() in CP_Rho_GSpacePlane)
        for(int ns=ist;ns<iend;ns++){
          CkAssert(ns<config.nstates);
          UgSpaceDriverProxy[proxyOffset](ns,thisIndex.x).startNonLocalEes(myTime);
        }//endfor
      }//endfor
    }//endif
#endif //NEW_FFT_DEBUG
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
  CkPrintf("[%d] In RhoRealSpacePlane[%d,%d] energyComp, Memory %.2lf MB\n", CkMyPe(),
      thisIndex.x, thisIndex.y, CmiMemoryUsage()/(1024.0 * 1024));
#endif

  //============================================================================

  int nf1          = myGrid_length[MY_A];
  int nf2          = myGrid_length[MY_B];
  int nf3          = myGrid_length[MY_C];
  int npts         = config.numFFTPoints;  // total number of points

  //============================================================================
  // Perform exchange correlation computation (no grad corr here).

  CPXCFNCTS::CP_exc_calc(npts, nf1, nf2, nf3, density, Vks, &exc_ret, &muxc_ret,
                         config.nfreq_xcfnctl);

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  double sum = 0.0;
  for (int i = 0; i < myGrid_size; i++) {
    sum += Vks[i];
    if(thisIndex.x + thisIndex.y == 0) {
      CkPrintf("%lf\n", Vks[i]);
    }
  }//endfor
  CkPrintf("[%d] RhoR[%d,%d] Sum_Vks %lf\n", CkMyPe(), thisIndex.x, thisIndex.y,
      sum);
#endif

  if(simReadOnly.cp_grad_corr_on == 0) {
    double exc[2];
    exc[0] = exc_ret;
    exc[1] = 0.0;
    contribute(2 * sizeof(double), exc, CkReduction::sum_double);
  }//endif

  //============================================================================
  // Invoke FFT to take rho(r) to rho(g) : do not over-write the density!!

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d] About to call FFT RhoRealSpacePlane[%d %d] FFT_RSpacetoGSpace, Memory %.2lf MB\n",
      CkMyPe(), thisIndex.x, thisIndex.y, CkMyPe(), CmiMemoryUsage()/(1024.0 * 1024));
#endif

  double  *dataR     = rhoIRX;   // rhoirx is around doing nothing now
  complex *dataC     = rhoIRXC;  // so we can use it to store the FFT

#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif

  //copy scaled density
  for(int i = 0; i < myGrid_size; i++) {
    dataR[i] = density[i] * volumeFactor;
  }

  Charm_doForwardFFT(CkCallback(CkIndex_CP_Rho_GSpacePlane::acceptRhoData(),
        UrhoGProxy[thisInstance.proxyOffset]),
        Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);

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

  if(simReadOnly.ees_nloc_on == 1) {
    //converting to broadcast
    if(thisIndex.x + thisIndex.y == 0) {
      for(int kpoint=0;kpoint < config.UberJmax;kpoint++)
      {
        UberCollection destKpointInstance = thisInstance;
        destKpointInstance.idxU.y = kpoint;
        int proxyOffset = destKpointInstance.setPO();
        UrealParticlePlaneProxy[proxyOffset].launchFFTControl(myTime);
      }
    }

#if 0  //old code
    if(thisIndex.x < simReadOnly.ngrid_nloc_c){
      int nstates =  config.nstates;
      int div     = (nstates/rhoRsubplanes);
      int rem     = (nstates % rhoRsubplanes);
      int add     = (thisIndex.y < rem ? 1 : 0);
      int max     = (thisIndex.y < rem ? thisIndex.y : rem);
      int ist     = div*thisIndex.y + max;
      int iend    = ist + div + add;
      for(int kpoint=0;kpoint < config.UberJmax;kpoint++)
      {
        UberCollection destKpointInstance=thisInstance;
        destKpointInstance.idxU.y=kpoint;
        int proxyOffset=destKpointInstance.setPO();
        for(int ns=ist;ns<iend;ns++){
          CkAssert(ns<config.nstates);
          UrealParticlePlaneProxy[proxyOffset](ns,thisIndex.x).launchFFTControl(myTime);
        }//endfor
      }//endfor
    }//endif
#endif
  }//endif
  //----------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *
 * Called when FFT of divRho is finished.
 *
 * Invoked 3 times per algorithm step : once for each grad_rho
 *
 * Memory required is : rho_igx,rho_igy,rho_igz so stuff can come in any order
 *                    : density and vks are needed later so no reusing for you.
 *                    : VksHart can also arrive at any time and cannot be used
 *                      here.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptGradRhoVks() {
  //============================================================================

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d] In RhoRealSpacePlane[%d,%d] acceptGradRhoVks_%d, Memory %.2lf MB\n",
      CkMyPe(), thisIndex.x, thisIndex.y, doneGradRhoVks, CmiMemoryUsage()/(1024.0 * 1024));
#endif
  doneGradRhoVks++;

  //============================================================================
  // When you have rhoiRX,rhoiRY,rhoiRZ and Vks invoke gradient correction

  if(doneGradRhoVks == 3) {
    doneGradRhoVks = 0;

#ifdef _CP_DEBUG_RHOG_VERBOSE_
    double sumVks = 0, sumX = 0, sumY = 0, sumZ = 0;
    for(int i = 0; i < myGrid_size; i++) {
      sumX += rhoIRX[i];
      sumY += rhoIRY[i];
      sumZ += rhoIRZ[i];
      sumVks += Vks[i];
    }
    CkPrintf("[%d] RhoR[%d,%d] %d div_sums %lf %lf %lf %lf\n", CkMyPe(),
        thisIndex.x, thisIndex.y, doneHartVks, sumX, sumY, sumZ, sumVks);
#endif

   //The gradient of the density is now completed. You can compute the
   //GGA-DFT functional now. Density is now available to be used as scratch.

    if(simReadOnly.cp_grad_corr_on == 0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Don't come in the grad corr routines when\n");
      CkPrintf("gradient corrections are off\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  double sumVks = 0, sumX = 0, sumY = 0, sumZ = 0;
  for(int i = 0; i < ngridb*(ngrida+2); i++) {
    sumX += rhoIRX[i];
    sumY += rhoIRY[i];
    sumZ += rhoIRZ[i];
    sumVks += Vks[i];
  }
  CkPrintf("[%d] RhoR[%d,%d] %d div_sums %lf %lf %lf %lf\n", CkMyPe(), thisIndex.x,
      thisIndex.y, doneHartVks, sumX, sumY, sumZ, sumVks);
#endif

#ifdef _CP_DEBUG_RHOR_VKSC_
    char myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "BGradRho_Real_%d_%d.out", thisIndex.x, thisIndex.y);
    FILE *fp = fopen(myFileName,"w");
    for (int i = 0; i < myGrid_size; i++){
      fprintf(fp, "%g %g %g %g\n", rhoIRX[i], rhoIRY[i], rhoIRZ[i], Vks[i]);
    }//endfor
    fclose(fp);
#endif

    //============================================================================
    // Compute the gradient corrected functional : Density is toast after this.
    int nf1                    = myGrid_length[MY_A];
    int nf2                    = myGrid_length[MY_B];
    int nf3                    = myGrid_length[MY_C];
    int npts                   = config.numFFTPoints;  // total number of points

#define GGA_ON
#ifdef GGA_ON
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif
    CPXCFNCTS::CP_getGGAFunctional(npts, nf1, nf2, nf3, density, rhoIRX, rhoIRY,
        rhoIRZ, Vks, thisIndex.x, &exc_gga_ret, config.nfreq_xcfnctl);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(GradCorrGGA_, StartTime, CmiWallTimer());
#endif
#endif // GGA ON

#ifdef _CP_DEBUG_RHOG_VERBOSE_
    double sumVks2 = 0.0, sumX2 = 0, sumY2 = 0, sumZ2 = 0;
    for(int i = 0; i < myGrid_size; i++) {
      sumX2 += rhoIRX[i];
      sumY2 += rhoIRY[i];
      sumZ2 += rhoIRZ[i];
      sumVks2 += Vks[i];
    }
    CkPrintf("[%d] RhoR[%d,%d] %d after_gga %lf %lf %lf %lf\n", CkMyPe(),
        thisIndex.x, thisIndex.y, doneHartVks, sumX2, sumY2, sumZ2, sumVks2);
#endif


    //============================================================================
    // Reduce the exchange correlation energy

    double exc[2];
    exc[0] = exc_ret;
    exc[1] = exc_gga_ret;
    contribute(2*sizeof(double), exc, CkReduction::sum_double);

    //============================================================================
    // output

#ifdef _CP_DEBUG_RHOR_VKSD_
    myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "AGradRho_Real_%d_%d.out", thisIndex.x, thisIndex.y);
    fp = fopen(myFileName,"w");
    for (int i = 0; i < myGrid_size; i++){
      fprintf(fp,"%g %g %g %g\n", rhoIRX[i], rhoIRY[i], rhoIRZ[i], Vks[i]);
    }//endfor
    fclose(fp);
#endif

    //============================================================================
    // Start the white bird : back fft of rhoirx, rhoiry, rhoirz

    whiteByrdFFT();
  }
  //---------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
  Send to RhoGspacePlane. RhoGspacePlane sends you back back another term.
  After this routine, rhoIRX, rhoIRY and rhoIRZ are `free'.

  Invoked once per algorithm step

 **/
//============================================================================
void CP_Rho_RealSpacePlane::whiteByrdFFT(){
  //============================================================================
  // Constants and pointers

  //============================================================================
  // I) scale, real to complex FFT, perform FFT

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d] In RhoRealSpacePlane[%d,%d] whiteByrdFFT, Memory %.2lf MB\n",
      CkMyPe(), thisIndex.x, thisIndex.y, CmiMemoryUsage()/(1024.0 * 1024));
#endif
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif

  for(int i = 0; i < myGrid_size; i++) {
    rhoIRX[i] *= FFTscale;
    rhoIRY[i] *= FFTscale;
    rhoIRZ[i] *= FFTscale;
  }

#ifdef  _CP_DEBUG_RHOG_VERBOSE_
    double sumVks2 = 0.0, sumX2 = 0, sumY2 = 0, sumZ2 = 0;
    for(int i = 0; i < myGrid_size; i++) {
      sumX2 += rhoIRX[i];
      sumY2 += rhoIRY[i];
      sumZ2 += rhoIRZ[i];
      sumVks2 += Vks[i];
      if(thisIndex.x + thisIndex.y == 0) {
        CkPrintf("%lf\n", Vks[i]);
      }
    }
    CkPrintf("[%d] RhoR[%d,%d] %d before white_fft %lf %lf %lf %lf\n", CkMyPe(),
        thisIndex.x, thisIndex.y, doneHartVks, sumX2, sumY2, sumZ2, sumVks2);
#endif

  Charm_doForwardFFT(CkCallback(CkIndex_CP_Rho_GSpacePlane::acceptWhiteByrd(),
        UrhoGProxy[thisInstance.proxyOffset]),
        Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);
  Charm_doForwardFFT(CkCallback(CkIndex_CP_Rho_GSpacePlane::acceptWhiteByrd(),
        UrhoGProxy[thisInstance.proxyOffset]),
        Urho_fft_yProxy[thisInstance.proxyOffset], fft_yoffset);
  Charm_doForwardFFT(CkCallback(CkIndex_CP_Rho_GSpacePlane::acceptWhiteByrd(),
        UrhoGProxy[thisInstance.proxyOffset]),
        Urho_fft_zProxy[thisInstance.proxyOffset], fft_zoffset);

#ifdef  _CP_DEBUG_RHOR_VERBOSE_
  double sumVks2 = 0, sumX2 = 0, sumY2 = 0, sumZ2 = 0;
  double *Vks                = rho_rs.Vks;
  for(int i = 0; i < ngridb*(ngrida+2); i++) {
    sumX2 += rhoIRX[i];
    sumY2 += rhoIRY[i];
    sumZ2 += rhoIRZ[i];
    sumVks2 += Vks[i];
    if(thisIndex.x + thisIndex.y == 0) {
      CkPrintf("%lf\n", Vks[i]);
    }
  }
  CkPrintf("[%d] RhoR[%d,%d] %d before white_fft %lf %lf %lf %lf\n", CkMyPe(), thisIndex.x,
      thisIndex.y, doneHartVks, sumX2, sumY2, sumZ2, sumVks2);
#endif


  //============================================================================
}
//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
  The white bird vks correction has returned from RhoG
  Invoked once per algorithm step
 **/
//============================================================================
void CP_Rho_RealSpacePlane::acceptWhiteByrd() {

  //============================================================================
  // Add the whitebyrd contrib to vks

#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d] In RhoRealSpacePlane[%d,%d] acceptWhiteByrd, Memory %.2lf MB\n",
      CkMyPe(), thisIndex.x, thisIndex.y, CmiMemoryUsage()/(1024.0 * 1024));
#endif
  double *dataR  = rhoIRX;  // whitebyrd correction stored here

#ifdef  _CP_DEBUG_RHOR_VERBOSE_
  double sum = 0.0;
  for (int i = 0; i < myGrid_size; i++) {
    sum += dataR[i];
    if(thisIndex.x + thisIndex.y == 0) {
      CkPrintf("%lf\n", dataR[i]);
    }
  }//endfor
  CkPrintf("[%d] RhoR[%d,%d] Sum_Vks acceptWhiteByrd %lf\n", CkMyPe(), thisIndex.x,
      thisIndex.y, sum);
#endif
  for(int i = 0; i < myGrid_size;i++) {
    Vks[i] -= dataR[i];
  }

  //============================================================================
  // Are we done yet?
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
 * Accept hartExt - FFT back from GHart
 * Invoked once per algorithm step.
 */
//============================================================================
void CP_Rho_RealSpacePlane::acceptHartVks(){
  //============================================================================
#ifdef _CP_DEBUG_RHOR_VERBOSE_
  CkPrintf("[%d] In RhoRealSpacePlane[%d,%d] acceptHartVks, Memory %.2lf MB\n",
      CkMyPe(), thisIndex.x, thisIndex.y, CmiMemoryUsage()/(1024.0 * 1024));
#endif
#ifdef  _CP_DEBUG_RHOR_VERBOSE_
    double sum = 0.0;
    for (int i = 0; i < myGrid_size; i++) {
      sum += VksHart[i];
      if(thisIndex.x + thisIndex.y == 0) {
        CkPrintf("%lf\n", VksHart[i]);
      }
    }//endfor
    CkPrintf("[%d] RhoR[%d,%d] Sum_Vks acceptHartVks %lf\n", CkMyPe(), thisIndex.x,
        thisIndex.y, sum);
#endif
  for(int i = 0; i < myGrid_size; i++)
  {
    Vks[i] += VksHart[i];
  }

  //============================================================================
  // Are we done yet?
  doneHartVks = true;
  doMulticastCheck();

  //============================================================================
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
  Under ees-eext Rhart chare reports its completion :  Set the done flag.
 **/
//============================================================================
void CP_Rho_RealSpacePlane::RHartReport(){
  doneRHart = true;
  doMulticastCheck();
  countRHart = 0;
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
  If all the parts of exc-eext-hart are done, invoke blast of vks to states
 **/
//============================================================================
void CP_Rho_RealSpacePlane::doMulticastCheck(){
  //============================================================================

  if(doneWhiteByrd && doneRHart && doneHartVks){
#ifdef  _CP_DEBUG_RHOR_VERBOSE_
    char myFileName[MAX_CHAR_ARRAY_LENGTH];
    sprintf(myFileName, "vks_final_%d_%d.out", thisIndex.x, thisIndex.y);
    FILE *fp = fopen(myFileName,"w");
    double sum = 0.0;
    for (int i = 0; i < myGrid_size; i++) {
      sum += Vks[i];
      fprintf(fp, "%lf\n", Vks[i]);
    }//endfor
    CkPrintf("[%d] RhoR[%d,%d] Sum_Vks final %lf\n", CkMyPe(), thisIndex.x,
        thisIndex.y, sum);
    fclose(fp);
#endif
    doMulticast();
  }

  //============================================================================
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
  Send vks back to the states
 **/
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
  countRHart = 0;
  doneWhiteByrd    = false;
  doneHartVks      = false;
  doneRHart        = false;

  //============================================================================
  // Send vks back to the states in real space

  int dataSize    = myGrid_length[MY_B] * myGrid_length[MY_A];

  if (config.useGMulticast) {

    for(int c_plane = 0; c_plane < myGrid_length[MY_C]; c_plane++) {

      VksMsg *msg = new (dataSize, 0) VksMsg;
      msg->pencil_offset_y       = thisIndex.y;
      double *dataR   = msg->data;

      int offset = c_plane * (myGrid_length[MY_B] * 2 * (myGrid_length[MY_A]/2 + 1));
      int jump = 2 - (myGrid_length[MY_A] & 1);
      int dest_off = 0;

      for(int b = 0; b < myGrid_length[MY_B]; b++) {
        for(int a = 0; a < myGrid_length[MY_A]; a++) {
          dataR[dest_off] = Vks[offset];
          offset++;
          dest_off++;
        }//endfor
        offset += jump;
      }//endfor

#ifdef _CP_DEBUG_RHOR_VKSE_
      char myFileName[MAX_CHAR_ARRAY_LENGTH];
      sprintf(myFileName, "vks_rho_Real_%d_%d.out", thisIndex.x, thisIndex.y);
      FILE *fp = fopen(myFileName,"w");
      //TODO : fix this to suit the new layout
      for (int j = 0, iii = 0; j < myNgridb; j++){
        for (int i = 0; i <ngrida; i++,iii++){
          fprintf(fp,"%d %d %g\n",i,j+myBoff,dataR[iii]);
        }}//endfor
      fclose(fp);
      contribute(CkCallback(CkIndex_CP_Rho_RealSpacePlane::exitForDebugging(),
            UrhoRealProxy[thisInstance.proxyOffset](0,0)));
#else
      for(int kp = 0; kp < config.UberJmax; kp++) {
        VksMsg *loopm;
        if(kp+1 < config.UberJmax) {
          // copying a message is perilous
          loopm = new (dataSize, 0) VksMsg;
          loopm->pencil_offset_y       = thisIndex.y;
          memcpy(loopm->data, msg->data, dataSize * sizeof(double));
        } else {
          loopm = msg;
        }
        realSpaceSectionProxyA[kp][c_plane].acceptVks(loopm);
      } //endfor kp
    }
#ifdef _CP_SUBSTEP_TIMING_
    if(rhoKeeperId>0)
    {
      double rhoend=CmiWallTimer();
      contribute(sizeof(double), &rhoend, CkReduction::max_double,
          CkCallback(CkIndex_TimeKeeper::collectEnd(NULL), 0, TimeKeeperProxy),
          rhoKeeperId);
    }
#endif
#endif //else of _CP_DEBUG_RHOR_VKSE_
  }//endif
  //============================================================================
}//end routine
//============================================================================


//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
  RhoReal exit for debugging
 **/
//============================================================================
void CP_Rho_RealSpacePlane::exitForDebugging(){
  CkPrintf("I am in the exitfordebuging rhoreal. Bye-bye\n");
  CkExit();
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *  Pup my variables for migration
 **/
//============================================================================
void CP_Rho_RealSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);
  PUParray(p, myGrid_length, 3);
  PUParray(p, myGrid_start, 3);
  PUParray(p, myGrid_end, 3);
  p|myGrid_size;
  p|myTime;

  p|rhoKeeperId;
  p|doneGradRhoVks;
  p|doneWhiteByrd;
  p|doneHartVks;
  p|doneRHart;
  p|countRHart;
  p|countRHartValue;
  p|FFTscale;
  p|volumeFactor;
  p|probScale;

  p|exc_ret;
  p|muxc_ret;
  p|exc_gga_ret;

  //TODO: pup for rest to be written
  //redCount, num_redn_complete, RedMsg
  //realSpaceSectionProxyA
  //Vks, density, rhoIR{X,Y,Z}, VksHart
  //---------------------------------------------------------------------------
}//end routine
//============================================================================

/*@}*/

