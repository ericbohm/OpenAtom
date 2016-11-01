//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_GSpacePlane.C
 ** @addtogroup Density
 *  @{
 *
 *  This is a description of the "life" of a CP_Rho_GSpacePlane  object
 *
 *  At the start of the program, the constructor CP_Rho_GSpacePlane is called.
 *  The RealSpaceDensity objects send data to CP_Rho_GSpacePlane using the
 *  acceptData() method. Inverse ffts in z and x directions are performed
 *  before the data is received, so here inverse fft in y direction is
 *  performed. This data is processed using the CP_hart_eext_calc. Then forward
 *  fft in the y direction is performed and data is send back to
 *  RealSpaceDensity objects.
 *
 *  The CP_Rho_GSpacePlaneHelper objects essentially help to split the work involved
 *  in g-space density computation. They receive their share of the work
 *  through the method recvCP_Rho_GSpacePlanePart() and send the processed
 *  data back to CP_Rho_GSpacePlane objects using the recvProcessedPart() method.
 *
 */
//============================================================================

#include <iostream>
#include <fstream>
#include <cmath>

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "cp_state_ctrl/CP_State_Plane.h"

#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "fft_charm.h"

//============================================================================

extern Config                                   config;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>     UrhoRealProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>           UrhoGHartExtProxy;
extern CkVec <CProxy_CP_State_GSpacePlane>      UgSpacePlaneProxy;
extern CkVec <CProxy_GSpaceDriver>              UgSpaceDriverProxy;
extern CkVec <CProxy_CP_Rho_GSpacePlane>        UrhoGProxy;
extern CkVec <CProxy_fft2d>                     Urho_fft_xProxy, Urho_fft_yProxy,
                                                Urho_fft_zProxy, Urho_fft_hartProxy;
extern CPcharmParaInfo                          simReadOnly;

//#define _CP_DEBUG_RHOG_VERBOSE_ 1
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::CP_Rho_GSpacePlane(UberCollection _instance) :
    thisInstance(_instance)
{
#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("{%d} Rho GS [%d] constructor\n", thisInstance.proxyOffset, thisIndex);
#endif

  myTime        = 1;
  divRhoX  = NULL;
  divRhoY  = NULL;
  divRhoZ  = NULL;

  //RAZ:  Added spin flags;
  cp_lsda              = simReadOnly.cp_lsda;
  mySpinIndex          = thisInstance.idxU.s;

  if(mySpinIndex==1 && cp_lsda!=1){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Error in spin index and LSDA flag in Hartree call.\n");
    CkPrintf("Good-bye.\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  doneWhiteByrd = 0;
  usesAtSync = true;
  if(config.lbdensity){
    setMigratable(true);
  }else{
    setMigratable(false);
  }
  //---------------------------------------------------------------------------
}//end routine
//============================================================================

//post constructor initialization
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::init()
{
  //query fft library to find my extents
  Charm_getOutputIndex(thisIndex, fft_xoffset,
      Urho_fft_xProxy[thisInstance.proxyOffset]);
  Charm_getOutputIndex(thisIndex, fft_yoffset,
      Urho_fft_yProxy[thisInstance.proxyOffset]);
  Charm_getOutputIndex(thisIndex, fft_zoffset,
      Urho_fft_zProxy[thisInstance.proxyOffset]);
  Charm_getOutputExtents(myGrid_start[MY_X], myGrid_end[MY_X],
                        myGrid_start[MY_Y], myGrid_end[MY_Y],
                        myGrid_start[MY_Z], myGrid_end[MY_Z],
                        Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);
  Charm_getAllPoints(myPoints, Urho_fft_xProxy[thisInstance.proxyOffset],
      fft_xoffset);

  for(int i = 0; i < 3; i++) {
    myGrid_length[i] = myGrid_end[i] - myGrid_start[i];
  }

  numPoints = (*myPoints).size();

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("{%d} Rho GS [%d] init: grid size %d %d %d, numPoints %d\n",
      thisInstance.proxyOffset, thisIndex, myGrid_length[0], myGrid_length[1], myGrid_length[2],
      numPoints);
#endif

  int data_size = myGrid_length[MY_X] * myGrid_length[MY_Y];
  myGrid_size = data_size;

  divRhoX   = (complex *) fftw_malloc(data_size * sizeof(complex));
  Charm_setOutputMemory((void*)divRhoX,
      Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);
  Charm_createOutputPlan(Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);

  divRhoY   = (complex *) fftw_malloc(data_size * sizeof(complex));
  Charm_setOutputMemory((void*)divRhoY,
      Urho_fft_yProxy[thisInstance.proxyOffset], fft_yoffset);
  Charm_createOutputPlan(Urho_fft_yProxy[thisInstance.proxyOffset], fft_yoffset);

  divRhoZ   = (complex *) fftw_malloc(data_size * sizeof(complex));
  Charm_setOutputMemory((void*)divRhoZ,
      Urho_fft_zProxy[thisInstance.proxyOffset], fft_zoffset);
  Charm_createOutputPlan(Urho_fft_zProxy[thisInstance.proxyOffset], fft_zoffset);
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::~CP_Rho_GSpacePlane(){
  fftw_free((void*)divRhoX);
  fftw_free((void*)divRhoY);
  fftw_free((void*)divRhoZ);
}
//============================================================================

//============================================================================
// This function when the density FFT from Rho_R is complete
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptRhoData() {
#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("{%d} Rho GS [%d] acceptRhoData\n", thisInstance.proxyOffset, thisIndex);
#endif
#ifdef _CP_DEBUG_RHOG_RHOG_
  char myFileName[100];
  sprintf(myFileName, "Rho_Gspace_%d.out", thisIndex);
  std::vector< gridPoint > & dpoints = (*myPoints);
  FILE *fp = fopen(myFileName,"w");
  for (int i = 0; i < numPoints; i++){
    fprintf(fp," %d %d %d : %g %g\n", dpoints[i].d3, dpoints[i].d2, dpoints[i].d1,
      divRhoX[dpoints[i].offset].re, divRhoX[dpoints[i].offset].im);
  }//endfor
  fclose(fp);
#endif
  //============================================================================
  // I) Communicate rho(g) to RHoGHartExt to compute eext and hart part of vks
  //     or print out that you are taking the day off

  //------------------------------------------------
  // Hartree is on, send rhoG to the harteext-G chare
#ifndef _CP_DEBUG_HARTEEXT_OFF_  // hartree is cooking
  //create first message to be able to copy
  RhoGHartMsg *msg, *base_msg = new (numPoints) RhoGHartMsg;
  complex *dest = base_msg->data;
  std::vector< gridPoint > & points = (*myPoints);
  for(int p = 0; p < numPoints; p++) {
    int off = points[p].offset;
    dest[p] = divRhoX[off];
  }

  int nchareHartAtmT = config.nchareHartAtmT;
  for(int j = 0; j < nchareHartAtmT; j++) {
    if(j != (nchareHartAtmT - 1)) {
      msg = new (numPoints) RhoGHartMsg;
      CmiMemcpy(msg->data, base_msg->data, numPoints * sizeof(complex));
    } else {
      msg = base_msg;
    }
    msg->size        = numPoints;
    //------------------------------------------------
    //RAZ: This is important:  
    // This is spin modified logic to send to up 
    // instance GHartExt Calc.
    // No way spin down should go there.
    // We are sending Grho from data_out
    //------------------------------------------------
    if (mySpinIndex==0){
      UrhoGHartExtProxy[thisInstance.proxyOffset](thisIndex,j).acceptData(msg);
    }else if(mySpinIndex==1){
      UberCollection upinstance=thisInstance;
      // flip the spin
      upinstance.idxU.s= !upinstance.idxU.s;
      // get new offset
      int Offset2RhoUp=upinstance.setPO();
      UrhoGHartExtProxy[Offset2RhoUp](thisIndex,j).acceptData(msg);
    }//endif spin logic

  }//endfor : atmType parallelization
#else // Hartree is off, chill
  if(thisIndex == 0) {
    CkPrintf("EHART       = OFF FOR DEBUGGING\n");
    CkPrintf("EExt        = OFF FOR DEBUGGING\n");
    CkPrintf("EWALD_recip = OFF FOR DEBUGGING\n");
  }//endif
#endif

  //============================================================================
  // III) Start grad corr computations if necessary

#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif
  if(simReadOnly.cp_grad_corr_on != 0) {
    divRhoVksGspace();
  }
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(divRhoVksGspace_, StartTime, CmiWallTimer());
#endif

  //============================================================================
  //kick off NL if its our job

  launchNlG();

  //---------------------------------------------------------------------------
}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// The density is here : Launch ees NL
// Do this once an algorithm step
//
//============================================================================
void CP_Rho_GSpacePlane::launchNlG() {
  //============================================================================
  // Launch the nonlocal energy computation

  if(simReadOnly.ees_nloc_on == 1 && config.launchNLeesFromRho == 2){
    CkAssert(config.UberJmax <= 1);
    //TODO: converted to broadcast, ensure correctness
    if(thisIndex == 0) {
     UgSpaceDriverProxy[thisInstance.proxyOffset].startNonLocalEes(myTime);
    }//endif
  }//endif

  //============================================================================
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Generate divX, divY, divZ, and begin the three FFTs
//============================================================================

void CP_Rho_GSpacePlane::divRhoVksGspace() {

  double tpi,*hmati;

  CPXCFNCTS::CP_fetch_hmati(&hmati,&tpi);

  memset(divRhoY, 0, sizeof(complex) * myGrid_size);
  memset(divRhoZ, 0, sizeof(complex) * myGrid_size);
  double gx,gy,gz;

  std::vector< gridPoint > & points = (*myPoints);
  double sumX = 0, sumY = 0, sumZ = 0;
  for(int p = 0; p < numPoints; p++) {
    int offset = points[p].offset;
    gx = tpi * (points[p].d3 * hmati[1] + points[p].d2 * hmati[2] +
        points[p].d1 * hmati[3]);
    gy = tpi * (points[p].d3 * hmati[4] + points[p].d2 * hmati[5] +
        points[p].d1 * hmati[6]);
    gz = tpi * (points[p].d3 * hmati[7] + points[p].d2 * hmati[8] +
        points[p].d1 * hmati[9]);
    complex tmp = (divRhoX[offset].multiplyByi())*(-1.0);
    divRhoX[offset] = tmp * gx;
    divRhoY[offset] = tmp * gy;
    divRhoZ[offset] = tmp * gz;
#if _CP_DEBUG_RHOG_VERBOSE_
    sumX += divRhoX[offset].re + divRhoX[offset].im;
    sumY += divRhoY[offset].re + divRhoY[offset].im;
    sumZ += divRhoZ[offset].re + divRhoZ[offset].im;
#endif
  }//endfor

#if _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("{%d} Rho GS [%d] divSums %lf %lf %lf\n", thisInstance.proxyOffset, thisIndex,
    sumX, sumY, sumZ);
#endif

  Charm_doBackwardFFT(CkCallback(CkIndex_CP_Rho_RealSpacePlane::acceptGradRhoVks(),
        UrhoRealProxy[thisInstance.proxyOffset]),
        Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset,
        1 / simReadOnly.vol);
  Charm_doBackwardFFT(CkCallback(CkIndex_CP_Rho_RealSpacePlane::acceptGradRhoVks(),
        UrhoRealProxy[thisInstance.proxyOffset]),
        Urho_fft_yProxy[thisInstance.proxyOffset], fft_yoffset,
        1 / simReadOnly.vol);
  Charm_doBackwardFFT(CkCallback(CkIndex_CP_Rho_RealSpacePlane::acceptGradRhoVks(),
        UrhoRealProxy[thisInstance.proxyOffset]),
        Urho_fft_zProxy[thisInstance.proxyOffset], fft_zoffset,
        1 / simReadOnly.vol);

  //---------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 *
 *  Callback for all 3 divRho FFT
 *
 *  We cannot receive whitebyrds until every divRho sent above has been received
 *  and processed by RhoReal. Therefore, our divRho memory is safe and warm
 *  while it it processing in the routines before the send. It is also safe
 *  during the send
*/
//============================================================================
void CP_Rho_GSpacePlane::acceptWhiteByrd() {
  //============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("{%d} Rho GS [%d] acceptWhiteByrd_%d\n", thisInstance.proxyOffset, thisIndex,
      doneWhiteByrd);
#endif

  doneWhiteByrd++;

  // When all 3 gradients are in g-space, then you will ready for the next step.
  if(doneWhiteByrd == 3){
    doneWhiteByrd = 0;
    /** The partially FFT'ed white byrd correction to VKS arrives to RhoG
      and ffts invoked. Only happens if gradient corrections are on.  */
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif

    //============================================================================
    // Compute my whiteByrd : store it in divrhox
    double tpi,*hmati;
    CPXCFNCTS::CP_fetch_hmati(&hmati, &tpi);

    double gx, gy, gz;
    complex *whitebyrd = divRhoX; // zeroing done carefully inside loop

    complex zero;
    zero.re = 0.0; zero.im = 0.0;
    std::vector< gridPoint > & points = (*myPoints);
    int last_offset = -1;
    for(int p = 0; p < numPoints; p++) {
      int offset = points[p].offset;
      if(offset != (last_offset + 1)) {
        for(int cur_off = last_offset + 1; cur_off < offset; cur_off++) {
          whitebyrd[cur_off] = zero;
        }
      }
      gx = tpi * (points[p].d3 * hmati[1] + points[p].d2 * hmati[2] +
          points[p].d1 * hmati[3]);
      gy = tpi * (points[p].d3 * hmati[4] + points[p].d2 * hmati[5] +
          points[p].d1 * hmati[6]);
      gz = tpi * (points[p].d3 * hmati[7] + points[p].d2 * hmati[8] +
          points[p].d1 * hmati[9]);
      complex tmp = divRhoX[offset]*gx + divRhoY[offset]*gy + divRhoZ[offset]*gz;
      whitebyrd[offset] = tmp.multiplyByi()*(-1.0);
      last_offset = offset;
    }

    for(int cur_off = last_offset + 1; cur_off < myGrid_size; cur_off++) {
      whitebyrd[cur_off] = zero;
    }

    Charm_doBackwardFFT(CkCallback(CkIndex_CP_Rho_RealSpacePlane::acceptWhiteByrd(),
          UrhoRealProxy[thisInstance.proxyOffset]),
          Urho_fft_xProxy[thisInstance.proxyOffset], fft_xoffset);
    myTime++;
  }
    //---------------------------------------------------------------------------
}//end routine
//============================================================================

//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Glenn's RhoG exit
 */
//============================================================================
void CP_Rho_GSpacePlane::exitForDebugging(){
  CkPrintf("I am in the exitfordebuging rhoG puppy. Bye-bye\n");
  CkExit();
}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::pup(PUP::er &p){
  ArrayElement1D::pup(p);
  PUParray(p, myGrid_length, 3);
  PUParray(p, myGrid_start, 3);
  PUParray(p, myGrid_end, 3);
  p|myGrid_size;
  p|numPoints;
  p|myTime;
  p|doneWhiteByrd;
  p|ehart_ret;
  p|eext_ret;
  p|ewd_ret;

  //TODO: pup needs to be done for following
  //myPoints, nlsectproxy
  //divRho{X,Y,Z}
  //--------------------------------------------------------------------------
}//end routine
//============================================================================
/*@}*/
