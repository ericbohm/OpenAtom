//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_GHartExc.C
 *
 *  This is a description of the "life" of a CP_Rho_GHartExc  object
 * 
 *  At the start of the program, the constructor CP_Rho_GHartExc is
 *  called.  We create our own rho_gs slab because we need the
 *  kvectors and the fft code.
 *
 *  Each iteration, the CP_Rho_GpacePlane object sends us the same rho
 *  data it will use in its own FFTs.  We then call the very intensive
 *  HartExtVksG function on the data.  We contribute the energy to the
 *  reduction, fft the vks, and ship the vks to rhoReal.
 *
 *  Then we're done until the next iteration.
 * 
 *  There is no RthThread control loop here because there is no
 *  meaningful flow of entry methods.  We get a message, we calculate,
 *  we send, we're done.  This object exists solely as a way to
 *  parallelize the uncomfortably long HartExtVks computation.
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

extern int sizeX;
extern Config config;

extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CPcharmParaInfoGrp    scProxy;
extern CProxy_AtomsGrp              atomsGrpProxy;
extern CProxy_CP_Rho_GHartExt       rhoGHartExtProxy;
extern CProxy_CP_Rho_RHartExt       rhoRHartExtProxy;
extern CProxy_eesCache              eesCacheProxy;
extern CProxy_FFTcache              fftCacheProxy;

extern ComlibInstanceHandle         commGHartInstance;
extern ComlibInstanceHandle         commGHartRHartIns0;
extern ComlibInstanceHandle         commGHartRHartIns1;

#define _CP_GHART_VERBOSE_OFF_

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
/**
 *  This object just gets a rho message, computes GHartExt, and sends
 *  vks.  
 */
//============================================================================
CP_Rho_GHartExt::CP_Rho_GHartExt(size2d sizeYZ,
          int _ngridaEext, int _ngridbEext, int _ngridcEext, int _ees_eext_on,
          int _natmTyp)
//============================================================================
   {//begin routine
//============================================================================

  CkAssert(sizeX>0); //check for startup wackiness
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
  CkVec <RunDescriptor> *sortedRunDescriptors = sim->RhosortedRunDescriptors;

//============================================================================

  ngridaEext  = _ngridaEext; 
  ngridbEext  = _ngridbEext; 
  ngridcEext  = _ngridcEext; 
  ees_eext_on = _ees_eext_on;
  natmTyp     = _natmTyp;

  iopt            = 0;
  iteration       = 0;
  iterAtmTyp      = 0;
  atmSFHere       = 0;
  densityHere     = 0;
  countEextFFT    = 0;
  ehart_ret       = 0.0;
  eext_ret        = 0.0;
  ewd_ret         = 0.0;
  registrationFlag= 0;
  launchFlag      = 0;
  nsendAtmTyp     = 0;


  rhoGHelpers     = config.rhoGHelpers;
  rho_gs.sizeX    = sizeX;
  rho_gs.sizeY    = sizeYZ[0];
  rho_gs.sizeZ    = sizeYZ[1];
  rho_gs.xdim     = rho_gs.sizeX;
  rho_gs.ydim     = rho_gs.sizeY;
  rho_gs.zdim     = 1;

//==================================================================================
// Decomposition rhoG lines into slices of size rhoGHelper

  ind_x            = thisIndex.x;
  ind_xdiv         = (ind_x/rhoGHelpers);
  ind_xrem         = (ind_x%rhoGHelpers);
  int numLines_tot = sortedRunDescriptors[ind_xdiv].size()/2;

  getSplitDecomp(&istrt_lines,&iend_lines,&numLines,numLines_tot,
                 rhoGHelpers,ind_xrem);

//==================================================================================
// Carve out your rundescriptor, make the k-vectors , malloc the memory

  rho_gs.numLines  = numLines;
  rho_gs.numRuns   = (numLines*2);
  rho_gs.numFull   = (numLines*rho_gs.sizeZ);
  rho_gs.size      = rho_gs.numFull;
  rho_gs.runs      = new RunDescriptor[(rho_gs.numRuns)];
  rho_gs.numPoints = 0;
  for (int r = (2*istrt_lines),s=0; r < (2*iend_lines); r++,s++) {
    rho_gs.numPoints += sortedRunDescriptors[ind_xdiv][r].length;
    rho_gs.runs[s]    = sortedRunDescriptors[ind_xdiv][r];
  }//endfor

  int nPacked;
  rho_gs.setKVectors(&nPacked);
  rho_gs.nPacked=nPacked;

  CkAssert(nPacked==rho_gs.numPoints);

  int numFull = rho_gs.numFull;
  numFullEext = numLines*ngridcEext;

  rho_gs.packedRho = (complex *)fftw_malloc(nPacked*sizeof(complex));
  rho_gs.Vks       = (complex *)fftw_malloc(nPacked*sizeof(complex));

  if(ees_eext_on==1){
    atmSF          = (complex *)fftw_malloc(numFullEext*sizeof(complex));
    atmSFtot       = (complex *)fftw_malloc(numFullEext*sizeof(complex));
  }//endif

  rho_gs.divRhoX   = NULL;
  rho_gs.divRhoY   = NULL;
  rho_gs.divRhoZ   = NULL;
  rho_gs.Rho       = NULL;

//==================================================================================

  if(ees_eext_on==1){
    eesCache *eesData  = eesCacheProxy.ckLocalBranch ();
    eesData->registerCacheGHart(thisIndex.x,nPacked,rho_gs.k_x,rho_gs.k_y,rho_gs.k_z);
    int i=1;
    CkCallback cb(CkIndex_CP_Rho_GHartExt::registrationDone(NULL),rhoGHartExtProxy);
    contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

//==================================================================================
// Set some proxies, set the migratable flag

  setMigratable(false);

  rhoRealProxy_com = rhoRealProxy;
  if(config.useGHartInsRhoRP){
     ComlibAssociateProxy(&commGHartInstance,rhoRealProxy_com);
  }//endif

  rhoRHartProxy_com0 = rhoRHartExtProxy;
  rhoRHartProxy_com1 = rhoRHartExtProxy;
  if(config.useGHartInsRHart){
     ComlibAssociateProxy(&commGHartRHartIns0,rhoRHartProxy_com0);
     ComlibAssociateProxy(&commGHartRHartIns1,rhoRHartProxy_com1);
  }//endif

  usesAtSync = CmiTrue;
  if(config.lbdensity){
    setMigratable(true);
  }else{
    setMigratable(false);
  }//endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
// Destructor
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GHartExt::~CP_Rho_GHartExt(){

  if(ees_eext_on==1){
    fftw_free(atmSF);
    fftw_free(atmSFtot);
  }//endif
  atmSF    = NULL;
  atmSFtot = NULL;

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::pup(PUP::er &p){
//============================================================================

  ArrayElement2D::pup(p);
  rho_gs.pup(p);

  p|ngridaEext;
  p|ngridbEext;
  p|ngridcEext;
  p|ees_eext_on;
  p|numFullEext;
  p|registrationFlag;
  p|launchFlag;

  p|natmTyp;
  p|iterAtmTyp;
  p|nsendAtmTyp;
  p|atmSFHere;
  p|densityHere;
  p|countEextFFT;
  p|ehart_ret;
  p|eext_ret;
  p|ewd_ret;

  p|iopt;
  p|iteration;
  p|ind_x;
  p|ind_xdiv;
  p|ind_xrem;
  p|rhoGHelpers; 
  p|istrt_lines; 
  p|iend_lines;
  p|numLines;
  p|rhoRealProxy_com;
  p|rhoRHartProxy_com0;
  p|rhoRHartProxy_com1;

  //  if(config.useCommlib)
  //    ComlibResetProxy(&rhoRealProxy_com);

//---------------------------------------------------------------------------
  }//endif
//============================================================================

//============================================================================
// The density arrives from RhoGspace
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::acceptData(RhoGHartMsg *msg){
//============================================================================
// Check the flow of control to see if we can use the data.

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int cp_min_opt = sim->cp_min_opt;

  if(cp_min_opt==0){
     if(atomsGrpProxy.ckLocalBranch()->iteration != iteration){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Flow of Control Error in HartExtVks : atoms slow %d %d\n",
              atomsGrpProxy.ckLocalBranch()->iteration,iteration);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
     }//endif
  }//endif

//============================================================================
// Copy out the data and flip arrival flag

  int ncoef  = rho_gs.numPoints;
  CkAssert(ncoef==msg->size);
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
    }
#endif

  memcpy(rho_gs.packedRho,msg->data,sizeof(complex)*ncoef);
  delete msg;  

  densityHere = 1;

//============================================================================
// If ees is off, go ahead with the N^2 method.
// If ees is on, either go ahead or chill depending on whether atmSF is here

  if(ees_eext_on==0){HartExtVksG();}
  if(ees_eext_on==1 && atmSFHere==1){getHartEextEes();}

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Compute hartree eext and vks using the N^2 method
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::HartExtVksG() { 
//============================================================================
// Get the variables

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in hartextvskg at %d\n",thisIndex.x);
#endif

   AtomsGrp *ag         = atomsGrpProxy.ckLocalBranch(); // find me the local copy
   int natm             = ag->natm;
   FastAtoms *fastAtoms = &(ag->fastAtoms);

   int numPoints = rho_gs.numPoints;
   int numFull   = rho_gs.numFull;
   complex *rho  = rho_gs.packedRho;
   complex *vks  = rho_gs.Vks;
   int *k_x      = rho_gs.k_x;
   int *k_y      = rho_gs.k_y;
   int *k_z      = rho_gs.k_z;

   ehart_ret  = 0.0;
   eext_ret  = 0.0;
   ewd_ret   = 0.0;

//============================================================================
// compute vks(g) from hart eext and reduce eext and ehart

#ifndef CMK_OPTIMIZE
   double  StartTime=CmiWallTimer();
#endif    

   CPLOCAL::CP_hart_eext_calc(numPoints,rho,natm,fastAtoms,vks,
                              &ehart_ret,&eext_ret,&ewd_ret,k_x,k_y,k_z,
                              thisIndex.x);
   iteration ++;
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(HartExcVksG_, StartTime, CmiWallTimer());    
#endif

#ifdef _CP_DEBUG_RHOG_VKSA_
     char myFileName[100];
     sprintf(myFileName, "Vks_Gspace_%d%d.out", thisIndex.x,thisIndex.y);
     FILE *fp = fopen(myFileName,"w");
       for (int i = 0; i < numPoints; i++){ 
              fprintf(fp," %d %d %d : %g %g\n",
                 k_x[i],k_y[i],k_z[i],
                 vks[i].re,vks[i].im);
       }//endfor
     fclose(fp);
#endif

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in hartextvskg completed at %d\n",thisIndex.x);
#endif

//============================================================================
// Reduce the energies computed then FFT Vks into Real space (part way)

   double e[3];
   e[0] = ehart_ret;
   e[1] = eext_ret;
   e[2] = ewd_ret;
   contribute(3 * sizeof(double),e,CkReduction::sum_double);

   FFTVks();

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// Partly fft vks(gx,gy,gz) -> vks(gx,gy,z)
//============================================================================
void CP_Rho_GHartExt::FFTVks() { 
//============================================================================
// Perform the FFT(gx,gy,gz) to FFT(gx,gy,z)

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in fftvks at %d\n",thisIndex.x,iterAtmTyp);
#endif

   FFTcache *fftcache = fftCacheProxy.ckLocalBranch();
   complex *vks       = rho_gs.Vks;
   complex *vksScr    = fftcache->tmpData; // scratch from the cache has FFT output

   fftcache->doHartFFTGtoR_Gchare(vks,vksScr,rho_gs.numFull,rho_gs.numPoints,numLines,
   			          rho_gs.numRuns,rho_gs.runs,rho_gs.sizeZ);

#ifdef _CP_DEBUG_RHOG_VKSA_
   int numFull        = rho_gs.numFull;
   sprintf(myFileName, "Vks_GspaceAFFT_%d%d.out", thisIndex.x,thisIndex.y);
   fp = fopen(myFileName,"w");
   for(int i = 0; i < numFull; i++){ 
     fprintf(fp," %g %g\n",VksExpd[i].re, VksExpd[i].im);
   }//endfor
   fclose(fp);
#endif

//============================================================================
//  Send partly ffted vks off to real space where FFT will be completed

   sendVks();

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Send vks_hart_ext back to rho_real where fft(gx,gy) will be performed
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::sendVks() { 
//============================================================================

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in sendvks at %d\n",thisIndex.x,iterAtmTyp);
#endif

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("Communicating data from RhoGHart to RhoR : %d %d\n",
	   thisIndex.x,thisIndex.y);
#endif

   FFTcache *fftcache = fftCacheProxy.ckLocalBranch();
   complex *vksScr    = fftcache->tmpData; // scratch from the cache has FFT output

//============================================================================
// Do a Comlib Dance

  if(config.useGHartInsRhoRP){commGHartInstance.beginIteration();}
  
//============================================================================

  int sizeZ=rho_gs.sizeZ;
  for(int z=0; z < sizeZ; z++) {

    RhoHartRSFFTMsg *msg = new (numLines,8*sizeof(int)) RhoHartRSFFTMsg;
    msg->size           = numLines;   // number of z-lines in this batch
    msg->senderBigIndex = ind_xdiv;   // big line batch index
    msg->senderStrtLine = istrt_lines;// where my lines start in big batch
    msg->iopt           = 0;          // iopt always 0 for us
    if(config.prioFFTMsg){
       CkSetQueueing(msg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(msg) = config.rhorpriority + thisIndex.x + thisIndex.y;
    }//endif

    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = vksScr[j];}
    rhoRealProxy_com(z,0).acceptHartVks(msg);
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
  }//endfor

//============================================================================
// Complete the commlib dance and hang out.
    
  if (config.useGHartInsRhoRP){commGHartInstance.endIteration();}

//---------------------------------------------------------------------------
  }// end routine
//============================================================================

//==========================================================================
// Make sure everyone is registered on the 1st time step
//==========================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==========================================================================
  void CP_Rho_GHartExt::registrationDone(CkReductionMsg *msg) {
//==========================================================================

  int sum = ((int *)msg->getData())[0];
  delete msg;

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("HI, I am Ghart chare %d in reg : %d %d\n",thisIndex.x,sum,natmTyp);
#endif

  registrationFlag=1;
  if(launchFlag==1){launchFlag=0;FFTEesBck();}

}
//==========================================================================


//============================================================================
// Recv Atm SF from RhoRhart : Euler Exponential spline based method
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::recvAtmSFFromRhoRHart(RhoGHartMsg *msg){
//============================================================================
// Unpack the message

  if(ees_eext_on==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees eext is off. You can't recvFromRhoRhart\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// Unpack the message and then delete it
//  atmSFT: numLines collections of lines of lth ngridcEext each with constant (gx,gy) 
//          Each message contains 1 pt on each line

  int size             = msg->size;
  int offset           = msg->offset;
  int iter             = msg->iter;
  complex *partlyIFFTd = msg->data;
#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++)
    {
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
    }
#endif

  if(iter!= (iterAtmTyp+1) ){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Atm Type Iteration out of whack. %d %d\n",iter,iterAtmTyp);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  CkAssert(numLines == size);

  if(countEextFFT==0){memset(atmSF, 0, sizeof(complex)*numFullEext);}
  for(int i=0,j=offset; i< numLines; i++,j+=ngridcEext){atmSF[j] = partlyIFFTd[i];}

  delete msg;

//============================================================================
// You must receive 1 message from each R-chare before continuing

  countEextFFT++;
  if (countEextFFT == ngridcEext) {
    countEextFFT = 0;
    launchFlag   = 1;
    if(registrationFlag==1){launchFlag=0; FFTEesBck();} 
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Finish FFting to G-space atmSF(gx,gy,z) -> atmSF(gx,gy,gz)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::FFTEesBck(){
//============================================================================
// Increment Counters set faux booleans

  atmSFHere = 1;
  iterAtmTyp ++;

  if(iterAtmTyp>1 && densityHere==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, the density has not arrived. You can't keep going. GhartEext\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// FFT yourself to heaven

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("GHart %d starting FFTBck at %d\n",thisIndex.x,iterAtmTyp);
#endif

  int numPoints = rho_gs.numPoints;
  fftCacheProxy.ckLocalBranch()->doEextFFTRtoG_Gchare(atmSF,numFullEext,numPoints,
                                 numLines,rho_gs.numRuns,rho_gs.runs,ngridcEext);

//============================================================================

  if(densityHere==1){getHartEextEes();}

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// compute HartreeEextEes
//============================================================================
void CP_Rho_GHartExt::getHartEextEes(){
//============================================================================

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("GHart %d starting getHartEextEes at %d\n",thisIndex.x,iterAtmTyp);
#endif

  if(ees_eext_on==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees eext is off. You can't getHartEextEes\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// Compute eext energy, hartree, total SF and VKS

  eesCache *eesData  = eesCacheProxy.ckLocalBranch ();

  int myChareG = thisIndex.x;
  int ncoef    = rho_gs.nPacked;
  complex *vks = rho_gs.Vks;
  complex *rho = rho_gs.packedRho;
  int *k_x     = rho_gs.k_x;
  int *k_y     = rho_gs.k_y;
  int *k_z     = rho_gs.k_z;
  double *b_re = eesData->RhoGHartData[myChareG].b_re;
  double *b_im = eesData->RhoGHartData[myChareG].b_im;

  if(iterAtmTyp==1){  
     bzero(vks,ncoef*sizeof(complex));  // no getting around these zeros
     bzero(atmSFtot,ncoef*sizeof(complex));
     ehart_ret = 0.0;
     eext_ret  = 0.0;
     ewd_ret   = 0.0;
  }//endif
#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif
  CPLOCAL::eesHartEextGchare(ncoef,iterAtmTyp,rho,vks,atmSF,atmSFtot,
                             b_re,b_im,&ehart_ret,&eext_ret,k_x,k_y,k_z,myChareG);
#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(eesHartExcG_, StartTime, CmiWallTimer());    
#endif

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in hartees: EHart %.10g Eext %.10g at %d\n",
	   thisIndex.x,ehart_ret,eext_ret,iterAtmTyp);
#endif

//============================================================================
// If you have SFtot get the ewald energy and FFT SFtot to real space

  if(iterAtmTyp==natmTyp){
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    CPLOCAL::eesEwaldGchare(ncoef,atmSFtot,b_re,b_im,&ewd_ret,k_x,k_y,k_z,myChareG);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(eesEwaldG_, StartTime, CmiWallTimer());    
#endif
#ifdef _CP_GHART_VERBOSE_
    CkPrintf("Ghart %d iter %d Here in hartees: EwaldRecip : %.10g\n",
	   thisIndex.x,iterAtmTyp,ewd_ret);
#endif
  }//endif

//============================================================================
// Blast out the energies when done

  if(iterAtmTyp==natmTyp){
#ifdef _CP_GHART_VERBOSE_
    CkPrintf("Ghart %d reduces energies at %d\n",iterAtmTyp);
#endif
    double e[3];
    e[0] = ehart_ret;
    e[1] = eext_ret;
    e[2] = ewd_ret;
    contribute(3*sizeof(double),e,CkReduction::sum_double);
  }//endif

//============================================================================
// Perform the back FFT SF and SFtot to get atm forces and VKS to get e-forces

  FFTEesFwd(0);                          // we always geneterate an SF
  if(iterAtmTyp==natmTyp){FFTVks();}     // Vks needs contribs from all SF
  if(iterAtmTyp==natmTyp){FFTEesFwd(1);} // Ewd needs the total SF

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Statr FFting to R-space atmSF(gx,gy,gz) -> atmSF(gx,gy,z)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::FFTEesFwd(int flag){
//============================================================================
#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in FFT to r: %.10g %.10g : with %d at %d %d\n",
            thisIndex.x,ehart_ret,eext_ret,flag,iterAtmTyp,natmTyp);
#endif

 //--------------------------------------------
 // Get the data and some scratch
  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
  complex *data_out  = fftcache->tmpData;
  complex *data;     if(flag==0){data=atmSF;}else{data=atmSFtot;}

 //--------------------------------------------
 // FFT the data : result goes into scratch
  fftcache->doEextFFTGtoR_Gchare(data,data_out,numFullEext,rho_gs.numPoints,numLines,
 			         rho_gs.numRuns,rho_gs.runs,ngridcEext);

 //--------------------------------------------
 // send yourself back to real space
  sendAtmSF(flag);

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// Send the SF data to back to Rhart to get atm forces
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::sendAtmSF(int flag){
//============================================================================
#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d Here in send to Rhart with opt %d at %d %d\n",
         thisIndex.x,flag,iterAtmTyp,natmTyp);
#endif

  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
  complex *senddata  = fftcache->tmpData;

//============================================================================
// start commlib

   switch(flag){
     case 0 : if(config.useGHartInsRHart) commGHartRHartIns0.beginIteration();break;
     case 1 : if(config.useGHartInsRHart) commGHartRHartIns1.beginIteration();break;
   }//endif

//============================================================================
// Send the message : 1 pt from each line to each chareR

  for(int z=0; z < ngridcEext; z++) {

    RhoRHartMsg *msg = new (numLines,8*sizeof(int)) RhoRHartMsg;
    msg->size        = numLines;     // number of z-lines in this batch
    msg->senderIndex = thisIndex.x;  // line batch index
    msg->iopt        = flag;         // iopt
    msg->iter        = iterAtmTyp;
    
    if(config.prioEextFFTMsg){
       CkSetQueueing(msg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(msg) = config.rhorHartpriority + thisIndex.x*10;
    }//endif

    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    for (int i=0,j=z; i<numLines; i++,j+=ngridcEext){data[i] = senddata[j];}

    switch(iopt){
      case 0 : rhoRHartProxy_com0(z,0).recvAtmForcFromRhoGHart(msg); break;
      case 1 : rhoRHartProxy_com1(z,0).recvAtmForcFromRhoGHart(msg); break;
    }//end switch
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif
  }//endfor

//============================================================================
// end commlib

   switch(flag){
     case 0 : if(config.useGHartInsRHart) commGHartRHartIns0.endIteration();break;
     case 1 : if(config.useGHartInsRHart) commGHartRHartIns1.endIteration();break;
   }//endif

//============================================================================
// We are done when when have sent out all SFs and the Ewald total SF

  nsendAtmTyp ++;
  if(nsendAtmTyp==(natmTyp+1)){
    nsendAtmTyp = 0;
    iterAtmTyp  = 0;
    atmSFHere   = 0;
    densityHere = 0; 
  }//endif

//============================================================================
   }//end routine
//============================================================================

