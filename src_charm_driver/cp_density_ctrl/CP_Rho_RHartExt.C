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

  ngrida      = _ngrida; 
  ngridb      = _ngridb; 
  ngridc      = _ngridc; 
  ees_eext_on = _ees_eext_on;
  natmTyp     = _natmTyp;

  countFFT[0]      = 0;
  countFFT[1]      = 0;
  nAtmTypRecv      = 0;
  registrationFlag = 0;
  launchFlag       = 0;
  iteration        = 0;
  iterAtmTyp       = 0;

  csize = (ngrida/2+1)*ngridb;  // complex variable size
  if(ees_eext_on==1){
    eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
    eesData->registerCacheRHart(thisIndex.x);

    atmSFC      = (complex*) fftw_malloc(csize*sizeof(complex));
    atmSFR      = reinterpret_cast<double*> (atmSFC);
    atmForcC    = atmSFC;
    atmForcR    = atmSFR;

    atmEwdSFC   = (complex*) fftw_malloc(csize*sizeof(complex));
    atmEwdSFR   = reinterpret_cast<double*> (atmEwdSFC);
    atmEwdForcC = atmEwdSFC;
    atmEwdForcR = atmEwdSFR;
    int i=1;
    CkCallback cb(CkIndex_CP_Rho_RHartExt::registrationDone(NULL),rhoRHartExtProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

  rhoGHartProxy_com    = rhoGHartExtProxy;;
  if (config.useRHartInsGHart){
    ComlibAssociateProxy(&commRHartGHartIns,rhoGHartProxy_com);          
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

   p|registrationFlag;
   p|launchFlag;
   p|ngrida;
   p|ngridb;
   p|ngridc;
   p|ees_eext_on;
   p|natmTyp;
   p|nAtmTypRecv;
   p(countFFT,2);
   p|iteration;
   p|iterAtmTyp;
   p|csize;

   if (p.isUnpacking()) {
     atmSFC      = (complex*) fftw_malloc(csize*sizeof(complex));
     atmSFR      = reinterpret_cast<double*> (atmSFC);
     atmForcC    = atmSFC;
     atmForcR    = atmSFR;
     atmEwdSFC   = (complex*) fftw_malloc(csize*sizeof(complex));
     atmEwdSFR   = reinterpret_cast<double*> (atmEwdSFC);
     atmEwdForcC = atmEwdSFC;
     atmEwdForcR = atmEwdSFR;
   }//endif

   p((char*)atmSFC,csize*sizeof(complex));
   p((char*)atmEwdSFC,csize*sizeof(complex));

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
  int itime        = iteration;

  eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
  if(iterAtmTyp==1){eesData->queryCacheRHart(myPlane,itime,iterAtmTyp);}

  int *plane_index = eesData->RhoRHartData[myPlane].plane_index;
  int **igrid      = eesData->RhoRHartData[myPlane].igrid;
  double **mn      = eesData->RhoRHartData[myPlane].mn;

  AtomsGrp *ag     = atomsGrpProxy.ckLocalBranch(); // find me the local copy
  int natm         = ag->natm;

  CPLOCAL::eesPackGridRchare(natm,iterAtmTyp,atmSFR,myPlane,igrid,mn,plane_index);

//============================================================================
// FFT the result to G-space

  FFTSFBck();

//============================================================================
  }//end routine 
//============================================================================


//============================================================================
// FFT of SF(x,y,z,iatmTyp) -> SF(gx,gy,z,iatmTyp)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::FFTSFBck(){

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in fftsback at %d\n",thisIndex.x,iterAtmTyp);
#endif

  int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_rho_x;
#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

  fftCacheProxy.ckLocalBranch()->doEextFFTRtoG_Rchare(atmSFC,atmSFR,nplane_x,
                                                      ngrida,ngridb);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(doEextFFTRtoG_, StartTime, CmiWallTimer());    
#endif

  SendAtmSFRhoGHart();
}
//============================================================================


//============================================================================
// Send SF(gx,gy,z,iatmTYP) to g-space whence the FFT will be completed
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::SendAtmSFRhoGHart(){

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in sendsf at %d\n",thisIndex.x,iterAtmTyp);
#endif

  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 
  int nchareG            = sim->nchareRhoGEext;
  int **tranpack         = sim->index_tran_upack_eext;
  int *nlines_per_chareG = sim->nlines_per_chareRhoGEext;

//===================================================================
// Perform the transpose and then the blast off the final 1D-FFT

  if(config.useRHartInsGHart){commRHartGHartIns.beginIteration();}

    for (int ic=0;ic<nchareG;ic++) { // chare arrays to which we will send
      int sendFFTDataSize = nlines_per_chareG[ic];
      RhoGHartMsg *msg    = new (sendFFTDataSize, 8 * sizeof(int)) RhoGHartMsg; 
      msg->iter           = iterAtmTyp;
      msg->size           = sendFFTDataSize;
      msg->offset         = thisIndex.x;    // c-plane-index
      complex *data       = msg->data;
      if(config.prioEextFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.rhogHartpriority+thisIndex.x*10;
      }//endif
      for(int i=0;i<sendFFTDataSize;i++){data[i] = atmSFC[tranpack[ic][i]];}
      rhoGHartProxy_com(ic,0).recvAtmSFFromRhoRHart(msg); // send the message
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
    }//end for : chare sending

  if(config.useRHartInsGHart){commRHartGHartIns.endIteration();}

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


  CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG            = sim->nchareRhoGEext;
  int **tranUnpack       = sim->index_tran_upack_eext;
  int *nlines_per_chareG = sim->nlines_per_chareRhoGEext;
   
  int size               = msg->size; 
  int Index              = msg->senderIndex;
  int iopt               = msg->iopt;
  int iter               = msg->iter;
  complex *partiallyFFTd = msg->data;
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
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

  if(size!=nlines_per_chareG[Index]){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Dude.1, %d != %d for rho rhart chare %d %d\n",size,nlines_per_chareG[Index],
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
  switch(iopt){
    case 0: data = atmSFC; break;
    case 1: data = atmEwdSFC; break;
  }//end switch

  if(countFFT[iopt]==1){memset(data,0,sizeof(complex)*csize);}
  for(int i=0;i<size;i++){data[tranUnpack[Index][i]] = partiallyFFTd[i];}

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
    FFTAtmForcFwd(iopt);
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Complete the FFT : sf(gx,gy,z) ->sf(x,y,z)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_RHartExt::FFTAtmForcFwd(int flagEwd){

#ifdef _CP_RHART_VERBOSE_
  CkPrintf("HI, I am rhart chare %d in FFTAtmForc w %d at %d \n",
        thisIndex.x,flagEwd,iterAtmTyp);
#endif

  double *dataR;   if(flagEwd==0){dataR=atmSFR;}else{dataR=atmEwdSFR;}
  complex *dataC;  if(flagEwd==0){dataC=atmSFC;}else{dataC=atmEwdSFC;}
  int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_rho_x;

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    
  
  fftCacheProxy.ckLocalBranch()->doEextFFTGtoR_Rchare(dataC,dataR,nplane_x,
                                                      ngrida,ngridb);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(doEextFFTGtoR_, StartTime, CmiWallTimer());    
#endif

  computeAtmForc(flagEwd);

}
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

  int myPlane      = thisIndex.x;
  int itime        = iteration;
  double *data;    if(flagEwd==0){data=atmSFR;}else{data=atmEwdSFR;}

  eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
  // you have already queried for this step:
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

  CPLOCAL::eesAtmForceRchare(natm,fastAtoms,nAtmTypRecv,igrid,dmn_x,dmn_y,dmn_z,
                             plane_index,data,myPlane,flagEwd);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(eesAtmForcR_, StartTime, CmiWallTimer());    
#endif

//============================================================================
// If you aren't done, start another round otherwise reset yourself

 // another round of SFs
 if(iterAtmTyp<natmTyp){computeAtmSF();}

#define _CP_DEBUG_STOP_RHART_OFF

#ifndef _CP_DEBUG_STOP_RHART_

 // we're done here
 if(iterAtmTyp==natmTyp && nAtmTypRecv==(natmTyp+1)){
#ifdef _CP_RHART_VERBOSE_
   CkPrintf("HI, I am rhart chare %d in computeAtmForc w %d at %d done\n",
             thisIndex.x,flagEwd,iterAtmTyp);
#endif
   iterAtmTyp  = 0;
   nAtmTypRecv = 0;

 }//endif

#else

 // we're done here and we exit
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
