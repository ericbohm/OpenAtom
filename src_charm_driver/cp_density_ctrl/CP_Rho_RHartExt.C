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
#include "main/eesCache.h"
#include "cp_state_ctrl/CP_State_Plane.h"

#include "fft_charm.h"

#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern CkVec <CProxy_CP_Rho_RealSpacePlane> UrhoRealProxy;
extern CkVec <CProxy_AtomsCache>            UatomsCacheProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>       UrhoRHartExtProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>       UrhoGHartExtProxy;
extern CkVec <CProxy_eesCache>              UeesCacheProxy;
extern Config                               config;
extern CPcharmParaInfo                      simReadOnly;
extern CkVec <CProxy_fft2d>                 Urho_fft_hartProxy;
extern CkVec < CkVec<CProxy_fft2d> >             Urho_fft_atmSFProxy;
extern CkVec <CkVec<CProxySection_CP_Rho_RHartExt> >        Urhart_sectionProxy;
extern CkGroupID                                mCastGrpId;

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
CP_Rho_RHartExt::CP_Rho_RHartExt(UberCollection _instance)
: thisInstance(_instance)
{
  if(simReadOnly.ees_eext_on == 0) {
    CkAbort("RHartExt created in non-EES method run\n");
  }

  int natmTypTot  = simReadOnly.natm_typ;
  int nchareHartAtmT = config.nchareHartAtmT;

  int div        = (natmTypTot / nchareHartAtmT);
  int rem        = (natmTypTot % nchareHartAtmT);
  int max        = (thisIndex.z < rem ? thisIndex.z : rem);

  natmTyp        = (thisIndex.z < rem ? div+1 : div);
  atmTypoff      = div*thisIndex.z + max;

  //============================================================================
  // Compute messages sizes and zero message counters

  nAtmTypRecv      = 0;
  registrationFlag = 0;
  launchFlag       = 0;
  iteration        = 0;
  iterAtmTyp       = 0;
  countDebug       = 0;
}
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
  Charm_getInputIndex(thisIndex.x, thisIndex.y, fft_atmSFOffset,
      Urho_fft_atmSFProxy[thisIndex.z][thisInstance.proxyOffset]);

  Charm_getInputExtents(myGrid_start[MY_X], myGrid_end[MY_X],
                        myGrid_start[MY_Y], myGrid_end[MY_Y],
                        myGrid_start[MY_Z], myGrid_end[MY_Z],
                        Urho_fft_atmSFProxy[thisIndex.z][thisInstance.proxyOffset],
                        fft_atmSFOffset);

  for(int i = 0; i < 3; i++) {
    myGrid_length[i] = myGrid_end[i] - myGrid_start[i];
  }

  eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  eesRHart_index = eesData->registerCacheRHart(thisIndex.x, thisIndex.y,
      myGrid_start[MY_C], myGrid_end[MY_C], myGrid_start[MY_B],
      myGrid_end[MY_B]);

  int data_size = myGrid_length[MY_C] * myGrid_length[MY_B] *
    (myGrid_length[MY_A]/2 + 1); //FFTW style storage
  myGrid_size = 2 * data_size; //how many real points I have

  atmSFC      = (complex*) fftw_malloc(myGrid_size * sizeof(complex));
  atmSFR      = reinterpret_cast<double*> (atmSFC);
  atmForcC    = atmSFC;
  atmForcR    = atmSFR;
  Charm_setInputMemory((void*)atmSFC,
      Urho_fft_atmSFProxy[thisIndex.z][thisInstance.proxyOffset],
      fft_atmSFOffset);
  Charm_createInputPlan(
      Urho_fft_atmSFProxy[thisIndex.z][thisInstance.proxyOffset],
      fft_atmSFOffset);

  if(thisIndex.z == 0) {
    atmEwdSFC   = (complex*) fftw_malloc(myGrid_size * sizeof(complex));
    atmEwdSFR   = reinterpret_cast<double*> (atmEwdSFC);
    atmEwdForcC = atmEwdSFC;
    atmEwdForcR = atmEwdSFR;
    Charm_getInputIndex(thisIndex.x, thisIndex.y, fft_atmSFTotOffset,
        Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset]);
    Charm_setInputMemory((void*)atmEwdSFC,
        Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset],
        fft_atmSFTotOffset);
    Charm_createInputPlan(
        Urho_fft_atmSFProxy[config.nchareHartAtmT][thisInstance.proxyOffset],
        fft_atmSFTotOffset);
  }
  
  if(config.nchareHartAtmT > 1) {
    if(thisIndex.x + thisIndex.y == 0) {
      CkMulticastMgr *mg = CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch();
      mg->resetSection(Urhart_sectionProxy[thisIndex.z][thisInstance.proxyOffset]);
    }
  }

  CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(), thisProxy);
  contribute(cb);

}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_RHartExt::~CP_Rho_RHartExt(){
  fftw_free((void*)atmSFC);
  if(thisIndex.z == 0) {
    fftw_free((void*)atmEwdSFC);
  }
}
//============================================================================

//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
/**
 * Make sure everyone is registered on the 1st time step
 */
//==========================================================================
void CP_Rho_RHartExt::registrationDone() {
  //==========================================================================

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("[%d] I am rhart chare %d %d in reg\n", CkMyPe(), thisIndex.x,
    thisIndex.y);
#endif

  registrationFlag = 1;
  if(launchFlag == 1) {
    computeAtmSF();
  }

}//end routine
//==========================================================================


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
  int cp_min_opt = simReadOnly.cp_min_opt;

  iteration ++;
  // the atoms haven't moved yet
  if(cp_min_opt == 0) {
    if(UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration
        != iteration-1) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Flow of Control Error in GHartExtVks : atoms slow %d %d\n",
          UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration,
          iteration);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

  //============================================================================
  // This is a new time step : Increment time step counter, zero atm typ counter
  //                           Launch if we are ready

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("[%d] I am rhart chare %d %d in start : reg flag %d\n", CkMyPe(),
      thisIndex.x, thisIndex.y, registrationFlag);
#endif

  if(iterAtmTyp == natmTyp) {
    CkPrintf("%d %d signing off\n", thisIndex.x, thisIndex.y);
    CkExit();
  }

  iterAtmTyp  = 0;
  nAtmTypRecv = 0;

  launchFlag = 1;
  if(registrationFlag == 1) {
    computeAtmSF();
  }

  //============================================================================
}//end routine
//============================================================================


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
  if(iterAtmTyp > natmTyp) {
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Too many iterations RHartExtVks\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("[%d] I am rhart chare %d %d %d in computeatmsf at iter %d\n",
      CkMyPe(), thisIndex.x, thisIndex.y, thisIndex.z, iterAtmTyp, CkMyPe());
#endif

  //============================================================================
  // Look up some constants and fill the r-space grid

  int itime          = iteration;
  int iterAtmTypFull = iterAtmTyp + atmTypoff;

  eesCache *eesData= UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  if(iterAtmTyp == 1){
    eesData->queryCacheRHart(eesRHart_index, itime, iterAtmTyp);
  }

  int **plane_index = eesData->RhoRHartData[eesRHart_index].plane_index;
  int **nSub       = eesData->RhoRHartData[eesRHart_index].nSub;
  int ***igrid     = eesData->RhoRHartData[eesRHart_index].igrid;
  double ***mn     = eesData->RhoRHartData[eesRHart_index].mn;

  // find me the local copy
  AtomsCache *ag   = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
  int natm         = ag->natm;

  CPLOCAL::eesPackGridRchare(natm, iterAtmTypFull, atmSFR, myGrid_start[MY_C],
      myGrid_end[MY_C], myGrid_start[MY_B], myGrid_end[MY_B], myGrid_size,
      igrid, mn, plane_index, nSub);

#ifdef _CP_RHART_VERBOSE_DUMP
  int off = 0;
  for(int c = myGrid_start[MY_C]; c < myGrid_end[MY_C]; c++) {
    for(int b = myGrid_start[MY_B]; b < myGrid_end[MY_B]; b++) {
      for(int a = 0; a < 2 * (myGrid_length[MY_A]/2 + 1); a++) {
        CkPrintf("%d %d %d %lf\n", c, b, a, atmSFR[off]);
        off++;
      }
    }
  }
#endif

  //============================================================================
  // FFT the result to G-space
  if(config.nchareHartAtmT > 1) {
    //use section proxy
    Charm_doForwardFFT(CkCallback(CkIndex_CP_Rho_GHartExt::doneAtmSF_FFT(),
          UrhoGHartExtProxy[thisInstance.proxyOffset](0,thisIndex.z)),
        Urho_fft_atmSFProxy[thisIndex.z][thisInstance.proxyOffset],
        fft_atmSFOffset);
  } else {
    Charm_doForwardFFT(CkCallback(CkIndex_CP_Rho_GHartExt::recvAtmSFFromRhoRHart(),
          UrhoGHartExtProxy[thisInstance.proxyOffset]),
        Urho_fft_atmSFProxy[thisIndex.z][thisInstance.proxyOffset],
        fft_atmSFOffset);
  }

  //============================================================================
}// end routine
//============================================================================

//============================================================================
/**
 *  Auxilary functions to mark the ending of FFT and informing the section
 */
//============================================================================

void CP_Rho_RHartExt::doneAtmSF_FFT() {
  if(thisIndex.x + thisIndex.y != 0) {
    CkAbort("Rho_RHartExt::doneAtmSF_FFT called on non-zero indexed chare\n");
  }
  FFT_Done_Msg* msg = new FFT_Done_Msg;
  Urhart_sectionProxy[thisIndex.z][thisInstance.proxyOffset].doneAtmSF_Multicast(msg);
}

void CP_Rho_RHartExt::doneAtmSF_Multicast(FFT_Done_Msg* msg) {
  delete msg;
  recvAtmForcFromRhoGHart();
}

void CP_Rho_RHartExt::doneAtmSFTot_FFT() {
  if(thisIndex.x + thisIndex.y + thisIndex.z != 0) {
    CkAbort("Rho_RHartExt::doneAtmSFTot_FFT called on non-zero indexed chare\n");
  }
  FFT_Done_Msg* msg = new FFT_Done_Msg;
  Urhart_sectionProxy[thisIndex.z][thisInstance.proxyOffset].doneAtmSFTot_Multicast(msg);
}

void CP_Rho_RHartExt::doneAtmSFTot_Multicast(FFT_Done_Msg *msg) {
  delete msg;
  recvAtmForcTotFromRhoGHart();
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * Hartree sends back atom forces from ext interation
 * Depending on the flag, it is Ewald or e-atm interation
 */
//============================================================================
void CP_Rho_RHartExt::recvAtmForcFromRhoGHart() {
  //============================================================================
  nAtmTypRecv++;
  computeAtmForc(0);
}//end routine
//============================================================================
void CP_Rho_RHartExt::recvAtmForcTotFromRhoGHart() {
  //============================================================================
  nAtmTypRecv++;
  computeAtmForc(1);
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

  // Ewald total SF  arrives on the last iteration only .
  // In parallel it can arrive before the last atmtyp SF
  if( (flagEwd == 1) && (iterAtmTyp != natmTyp)){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("You can't have flagEwd==1 unless you are on the last step\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  // If we are not on the last iteration when ewald can show up
  // the iteration # must match the number of SF received
  if( (flagEwd==0) && (iterAtmTyp!=nAtmTypRecv) && (iterAtmTyp!=natmTyp)){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Atm SF Recv and atm SF calc out of sync %d - %d %d %d\n",
        thisIndex.x, iterAtmTyp, nAtmTypRecv, natmTyp);
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

  int iterAtmTypFull = iterAtmTyp + atmTypoff;

  double *data;
  if(flagEwd == 0) {
    data = atmSFR;
  } else {
    data = atmEwdSFR;
  }

  // you have already queried for this step:
  eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
  AtomsCache *ag     = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();

  int **plane_index = eesData->RhoRHartData[eesRHart_index].plane_index;
  int **nSub       = eesData->RhoRHartData[eesRHart_index].nSub;
  int ***igrid     = eesData->RhoRHartData[eesRHart_index].igrid;
  double ***dmn_x    = eesData->RhoRHartData[eesRHart_index].dmn_x;
  double ***dmn_y    = eesData->RhoRHartData[eesRHart_index].dmn_y;
  double ***dmn_z    = eesData->RhoRHartData[eesRHart_index].dmn_z;

  FastAtoms *fastAtoms = &(ag->fastAtoms);
  int natm             = ag->natm;

#if CMK_TRACE_ENABLED
  double  StartTime=CmiWallTimer();
#endif

#ifdef _CP_RHART_VERBOSE_DUMP
  char myFileName[MAX_CHAR_ARRAY_LENGTH];
  sprintf(myFileName, "atmsf_%d_%d_%d.out", thisIndex.x, thisIndex.y, nAtmTypRecv);
  FILE *fp = fopen(myFileName,"w");
  int off = 0;
  for(int c = myGrid_start[MY_C]; c < myGrid_end[MY_C]; c++) {
    for(int b = myGrid_start[MY_B]; b < myGrid_end[MY_B]; b++) {
      for(int a = 0; a < 2 * (myGrid_length[MY_A]/2 + 1); a++) {
        fprintf(fp,"%d %d %d %lf\n", c, b, a, data[off]);
        off++;
      }
    }
  }
  fclose(fp);
#endif

  CPLOCAL::eesAtmForceRchare(natm, fastAtoms, iterAtmTypFull, igrid,
      dmn_x, dmn_y, dmn_z, plane_index, nSub, data, myGrid_start[MY_C],
      myGrid_end[MY_C], myGrid_start[MY_B], myGrid_end[MY_B], flagEwd);
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesAtmForcR_, StartTime, CmiWallTimer());
#endif

  //============================================================================
  // If you aren't done, start another round.

  if(iterAtmTyp < natmTyp){
    computeAtmSF();
  }

  //============================================================================
  // If your are done, Tell rho real and reset yourself for the next time step

  int nrecvMax = natmTyp;
  if(thisIndex.z == 0){
    nrecvMax++;
  }

#ifndef _CP_DEBUG_STOP_RHART_
  if(iterAtmTyp == natmTyp && nAtmTypRecv == nrecvMax) {
    contribute(CkCallback(CkIndex_CP_Rho_RealSpacePlane::RHartReport(),
        UrhoRealProxy[thisInstance.proxyOffset]));
    iterAtmTyp  = 0;
    nAtmTypRecv = 0;
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
  CkExit();
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::pup(PUP::er &p){
  //============================================================================
  ArrayElement3D::pup(p);
  p|countDebug;
  p|registrationFlag;
  p|launchFlag;
  p|atmTypoff;
  p|natmTyp;
  p|iteration;
  p|iterAtmTyp;
  p|nAtmTypRecv;

  //TODO fill in remaining pup
  //============================================================================
}//end routine
//============================================================================

/*@}*/
