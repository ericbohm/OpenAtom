//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_GHartExt.C
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

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "main/AtomsCache.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "main/eesCache.h"
#include "cp_state_ctrl/CP_State_Plane.h"

#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"
#include "src_piny_physics_v1.0/include/proto_defs/proto_cp_ewald_corrs.h"
#include "src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "main/PhysScratchCache.h"
#include <iostream>
#include <fstream>
#include <cmath>

//============================================================================

extern int sizeX;
extern Config config;

extern CkVec <CProxy_CP_Rho_RealSpacePlane> UrhoRealProxy;
extern CProxy_CPcharmParaInfoGrp    scProxy;
extern CProxy_PhysScratchCache     pScratchProxy;
extern CkVec <CProxy_AtomsCache>              UatomsCacheProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>       UrhoGHartExtProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>       UrhoRHartExtProxy;
extern CkVec <CProxy_eesCache>              UeesCacheProxy;
extern CkVec <CProxy_FFTcache>              UfftCacheProxy;

extern ComlibInstanceHandle         commGHartInstance;
extern ComlibInstanceHandle         commGHartRHartIns0;
extern ComlibInstanceHandle         commGHartRHartIns1;

//#define _DEBUG_INT_TRANS_FWD_
//#define _CP_GHART_VERBOSE_

/** @addtogroup Density
    @{
*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
/**
 *  This object just gets a rho message, computes GHartExt, and sends
 *  vks.  
 */
//============================================================================
CP_Rho_GHartExt::CP_Rho_GHartExt(
          int _ngridaEext, int _ngridbEext, int _ngridcEext, int _ees_eext_on,
          int _natmTyp, UberCollection _instance) :thisInstance(_instance)
//============================================================================
   {//begin routine
//============================================================================

  CkAssert(sizeX>0); //check for startup wackiness
  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  CkVec <RunDescriptor> *sortedRunDescriptors = sim->RhosortedRunDescriptors;

//============================================================================
// Fetch and compute the useful variables : zero counters

  rhoRsubplanes  = config.rhoRsubplanes;
  rhoGHelpers    = config.rhoGHelpers;
  nchareHartAtmT = config.nchareHartAtmT;
  iperd          = sim->iperd;

  ngridaEext     = _ngridaEext; 
  ngridbEext     = _ngridbEext; 
  ngridcEext     = _ngridcEext; 
  ees_eext_on    = _ees_eext_on;
  natmTypTot     = _natmTyp;
  CkAssert(natmTypTot>=nchareHartAtmT);

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
  countVksTot     = 0;
  countAtmSFtot   = 0;
  CountDebug      = 0;

  if(rhoRsubplanes>1){
     recvCountFromRHartExt = 0;
     for(int i=0;i<rhoRsubplanes;i++){
	if(sim->nline_send_eext_y[thisIndex.x][i]>0){recvCountFromRHartExt++;}
     }//endfor
     recvCountFromRHartExt*=ngridcEext;
  }else{
     recvCountFromRHartExt=ngridcEext;
  }//endif

//==================================================================================
// AtmTyp parallelization

  int div        = (natmTypTot/nchareHartAtmT); 
  int rem        = (natmTypTot % nchareHartAtmT);
  int max        = (thisIndex.y < rem ? thisIndex.y : rem);
  natmTyp        = (thisIndex.y<rem ? div+1 : div);  
  atmTypoff      = div*thisIndex.y + max;         

  if(ees_eext_on==0 && nchareHartAtmT>1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Parallel atom type not supported without ees Eext\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//==================================================================================
// Decomposition rhoG lines into slices of size rhoGHelper

  ind_x            = thisIndex.x;
  ind_xdiv         = (ind_x/rhoGHelpers);
  ind_xrem         = (ind_x%rhoGHelpers);
  int numLines_tot = sortedRunDescriptors[ind_xdiv].size()/2;

  getSplitDecomp(&istrt_lines,&iend_lines,&numLines,numLines_tot,
                 rhoGHelpers,ind_xrem);

//==================================================================================
// Set up rho_gs : Carve out your rundescriptor : Make the k-vectors  (donuts)

  rho_gs.sizeX     = sizeX;
  rho_gs.sizeY     = sim->sizeY;
  rho_gs.sizeZ     = sim->sizeZ;
  rho_gs.xdim      = rho_gs.sizeX;
  rho_gs.ydim      = rho_gs.sizeY;
  rho_gs.zdim      = 1;

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

//==================================================================================
// Malloc your memory : more juice for ees_ext and yet to support atmtyp parallel

  atmSF          = NULL;
  atmSFtot       = NULL;
  VksRecv        = NULL;
  atmSFtotRecv   = NULL;
  rho_gs.divRhoX = NULL;
  rho_gs.divRhoY = NULL;
  rho_gs.divRhoZ = NULL;
  rho_gs.Rho     = NULL;

  rho_gs.packedRho = (complex *)fftw_malloc(nPacked*sizeof(complex));
  rho_gs.Vks       = (complex *)fftw_malloc(nPacked*sizeof(complex));
  if(ees_eext_on==1 && nchareHartAtmT>1 && thisIndex.y==1){
    VksRecv        = (complex *)fftw_malloc(nPacked*sizeof(complex));
  }//endif

  if(ees_eext_on==1){
    atmSF          = (complex *)fftw_malloc(numFullEext*sizeof(complex));
    atmSFtot       = (complex *)fftw_malloc(numFullEext*sizeof(complex));
    if(nchareHartAtmT>1 && thisIndex.y==0){
      atmSFtotRecv = (complex *)fftw_malloc(numFullEext*sizeof(complex));
    }//endif
  }//endif


//==================================================================================
// Set some proxies, set the migratable flag

  setMigratable(false);

  usesAtSync = true;
  if(config.lbdensity){
    setMigratable(true);
  }else{
    setMigratable(false);
  }//endif

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// Post constructor initialization
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::init(){
//==================================================================================
// Register in the cache : contribute to a reduction to be sure everyone is done

  if(ees_eext_on==1){
     eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
     eesData->registerCacheGHart(thisIndex.x,rho_gs.nPacked,rho_gs.k_x,rho_gs.k_y,rho_gs.k_z);
     int i=1;
     CkCallback cb(CkIndex_CP_Rho_GHartExt::registrationDone(NULL),UrhoGHartExtProxy[thisInstance.proxyOffset]);
     contribute(sizeof(int),&i,CkReduction::sum_int,cb);
  }//endif

//==================================================================================
// Set up periodicity corrections if necessary

  rho_gs.iperd = iperd;
  if(iperd!=3){
    int nPacked     = rho_gs.nPacked;
    rho_gs.perdCorr = (double *)fftw_malloc(nPacked*sizeof(double));
    setput_nd_eext_corrs(nPacked,rho_gs.k_x,rho_gs.k_y,rho_gs.k_z,rho_gs.perdCorr);
  }//endif

//==================================================================================
// Set some proxies, set the migratable flag

  setMigratable(false);

  rhoRealProxy_com = UrhoRealProxy[thisInstance.proxyOffset];

#ifdef USE_COMLIB
  if(config.useGHartInsRhoRP){
     ComlibAssociateProxy(commGHartInstance,rhoRealProxy_com);
  }//endif
#endif

  if (ees_eext_on == 1)
  {
      rhoRHartProxy_com0 = UrhoRHartExtProxy[thisInstance.proxyOffset];
      rhoRHartProxy_com1 = UrhoRHartExtProxy[thisInstance.proxyOffset];
      #ifdef USE_COMLIB
        if(config.useGHartInsRHart)
        {
            ComlibAssociateProxy(commGHartRHartIns0,rhoRHartProxy_com0);
            ComlibAssociateProxy(commGHartRHartIns1,rhoRHartProxy_com1);
        }//endif
      #endif
  }
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
    if(nchareHartAtmT>1 && thisIndex.y==1){fftw_free(VksRecv);}
    if(nchareHartAtmT>1 && thisIndex.y==0){fftw_free(atmSFtotRecv);}
  }//endif
  atmSF        = NULL;
  atmSFtot     = NULL;
  VksRecv      = NULL;
  atmSFtotRecv = NULL;

}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::pup(PUP::er &p){
//============================================================================

  ArrayElement2D::pup(p);
  rho_gs.pup(p);

  p|iperd;
  p|rhoRsubplanes;
  p|ngridaEext;
  p|ngridbEext;
  p|ngridcEext;
  p|ees_eext_on;
  p|numFullEext;
  p|registrationFlag;
  p|launchFlag;
  p|natmTypTot;
  p|atmTypoff;
  p|nchareHartAtmT;
  p|natmTyp;
  p|iterAtmTyp;
  p|nsendAtmTyp;
  p|atmSFHere;
  p|densityHere;
  p|countEextFFT;
  p|ehart_ret;
  p|eext_ret;
  p|ewd_ret;
  p|countAtmSFtot;
  p|countVksTot;
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
  //    ComlibResetProxy(&UrhoRealProxy[thisInstance.proxyOffset]_com);

//---------------------------------------------------------------------------
  }//endif
//============================================================================

//============================================================================
// The density arrives from RhoGspace ONCE a time step
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::acceptData(RhoGHartMsg *msg){
//============================================================================
// Check the flow of control to see if we can use the data.

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int cp_min_opt = sim->cp_min_opt;

  iteration ++;

 // the atoms haven't moved yet
  if(cp_min_opt==0){
     if(UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration != iteration-1){
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkPrintf("Flow of Control Error in GHartExtVks : atoms slow %d %d\n",
              UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->iteration,iteration);
        CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
        CkExit();
     }//endif
  }//endif

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
  }
#endif

//============================================================================
// Copy out the data and flip arrival flag

  int ncoef  = rho_gs.numPoints;
  CkAssert(ncoef==msg->size);

  CmiMemcpy(rho_gs.packedRho,msg->data,sizeof(complex)*ncoef);
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
  CkPrintf("Ghart %d Here in hartextvskg at %d on %d\n",thisIndex.x,CkMyPe());
#endif

   AtomsCache *ag         = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch(); // find me the local copy
   int natm             = ag->natm;
   FastAtoms *fastAtoms = &(ag->fastAtoms);

   int numPoints = rho_gs.numPoints;
   int numFull   = rho_gs.numFull;
   complex *rho  = rho_gs.packedRho;
   complex *vks  = rho_gs.Vks;
   int *k_x      = rho_gs.k_x;
   int *k_y      = rho_gs.k_y;
   int *k_z      = rho_gs.k_z;

   ehart_ret     = 0.0;
   eext_ret      = 0.0;
   ewd_ret       = 0.0;

//============================================================================
// compute vks(g) from hart eext and reduce eext and ehart

#if CMK_TRACE_ENABLED
   double  StartTime=CmiWallTimer();
#endif    

   double *perdCorr = rho_gs.perdCorr;
   CPLOCAL::CP_hart_eext_calc(numPoints,rho,natm,fastAtoms,vks,
                              &ehart_ret,&eext_ret,&ewd_ret,k_x,k_y,k_z,perdCorr,
                              thisIndex.x, pScratchProxy.ckLocalBranch()->psscratch, config.nfreq_cplocal_hartext);
#if CMK_TRACE_ENABLED
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
  CkPrintf("Ghart %d Here in hartextvskg completed at %d on %d\n",thisIndex.x,CkMyPe());
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
  CkPrintf("Ghart %d %d Here in fftvks at %d on %d\n",
             thisIndex.x,thisIndex.y,iterAtmTyp,CkMyPe());
#endif

   FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
   fftcache->getCacheMem("CP_Rho_GHartExt::FFTVks");
   complex *vks       = rho_gs.Vks;
   complex *vksScr    = fftcache->tmpData; // scratch from the cache has FFT output

   fftcache->doHartFFTGtoR_Gchare(vks,vksScr,rho_gs.numFull,rho_gs.numPoints,numLines,
   			          rho_gs.numRuns,rho_gs.runs,rho_gs.sizeZ,ind_x);

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
  CkPrintf("Ghart %d %d Here in sendvks at %d on %d\n",
            thisIndex.x,thisIndex.y,iterAtmTyp,CkMyPe());
#endif

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("Communicating data from RhoGHart to RhoR : %d %d\n",
	   thisIndex.x,thisIndex.y);
#endif

//============================================================================

   CPcharmParaInfo *sim = CPcharmParaInfo::get();
   FFTcache *fftcache   = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
   complex *vksScr      = fftcache->tmpData; // scratch from the cache has FFT output

   int ix               = thisIndex.x;
   int sendLines        = numLines;

   int **index_pack;
   int *nline_send;
   if(rhoRsubplanes>1){
     nline_send  = sim->nline_send_eext_y[ix];
     index_pack  = sim->index_tran_pack_eext_ys[ix];
   }//endif

//============================================================================
// Do a Comlib Dance
#ifdef USE_COMLIB
   if(rhoRsubplanes==1){
#ifdef OLD_COMMLIB
    if(config.useGHartInsRhoRP){commGHartInstance.beginIteration();}
#else
    if(config.useGHartInsRhoRP){ComlibBegin(rhoRealProxy_com,0);}
#endif
   }//endif
#endif

//============================================================================

  int sizeZ=rho_gs.sizeZ;
  for(int z=0; z < sizeZ; z++) {
  for(int s=0;s<rhoRsubplanes;s++){

   if(rhoRsubplanes>1){sendLines = nline_send[s];} 
   if(sendLines >0)
      {

	RhoHartRSFFTMsg *msg = new (sendLines,8*sizeof(int)) RhoHartRSFFTMsg;
	msg->size           = sendLines;   // number of z-lines in this batch
	msg->index          = thisIndex.x;  // regular old index
	msg->senderBigIndex = ind_xdiv;     // big line batch index
	msg->senderStrtLine = istrt_lines;  // where my lines start in big batch
	msg->iopt           = 0;            // iopt always 0 for us
	if(config.prioFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = config.rhorpriority + thisIndex.x + thisIndex.y;
	}//endif

	// beam out all points with same z to chare array index z
	complex *data = msg->data;
	if(rhoRsubplanes==1){
	  for (int i=0,j=z; i<sendLines; i++,j+=sizeZ){data[i] = vksScr[j];}
	}else{
	  for(int i=0; i< sendLines; i++){
	    data[i] = vksScr[(z+index_pack[s][i])];
	  }//endif
	}//endif

	if(rhoRsubplanes==1){
	  rhoRealProxy_com(z,0).acceptHartVks(msg);
	}else{
	  UrhoRealProxy[thisInstance.proxyOffset](z,s).acceptHartVks(msg);
	}//endif

      } //endif
  }//endfor : subplanes
#ifdef CMK_BLUEGENEL
  CmiNetworkProgress();
#endif
  }//endfor : z plane parallel index

//============================================================================
// Complete the commlib dance and hang out.

#ifdef USE_COMLIB    
  if(rhoRsubplanes==1){
#ifdef OLD_COMMLIB
    if(config.useGHartInsRhoRP){commGHartInstance.endIteration();}
#else
    if(config.useGHartInsRhoRP){ComlibEnd(rhoRealProxy_com,0);}
#endif
  }//endif
#endif

  fftcache->freeCacheMem("CP_Rho_GHartExt::sendVks");

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
  CkPrintf("HI, I am Ghart chare %d in reg : %d %d on %d\n",
            thisIndex.x,sum,natmTyp,CkMyPe());
#endif

  registrationFlag=1;
  if(launchFlag==1){launchFlag=0;FFTEesBck();}

}//end routine
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

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  int size             = msg->size;
  int offset           = msg->offset;   // z index
  int isub             = msg->offsetGx; // subplane index
  int iter             = msg->iter;
  complex *partlyIFFTd = msg->data;
 
  int ix              = thisIndex.x;
  int numLinesNow     = numLines;
  int **index_pack;
  if(rhoRsubplanes>1){
    int *nline_send = sim->nline_send_eext_y[ix];
    index_pack      = sim->index_tran_pack_eext_y[ix];
    numLinesNow     = nline_send[isub];
  }//endif

//============================================================================
// Check for errors

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
  }//endfor
#endif

  if(iter!= (iterAtmTyp+1) ){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Atm Type Iteration out of whack. %d %d : chare %d %d : %d %d\n",
              iter,iterAtmTyp,thisIndex.x,thisIndex.y,natmTyp,atmTypoff);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  CkAssert(numLinesNow == size);

//============================================================================
// No need to zero : Every value is set.

  if(rhoRsubplanes==1){
    for(int i=0,j=offset; i< numLines; i++,j+=ngridcEext){atmSF[j] = partlyIFFTd[i];}
  }else{
    for(int i=0; i< numLinesNow; i++){
      atmSF[(offset+index_pack[isub][i])] = partlyIFFTd[i];
    }//endif
  }//endif

  delete msg;

//============================================================================
// You must receive 1 message from each R-chare before continuing

  countEextFFT++;
  if (countEextFFT == recvCountFromRHartExt) {
#ifdef _DEBUG_GLENN_KPT_
    char name[100];
    sprintf(name,"SfAtmGxGyZ_inG.p%d.t%d.out",thisIndex.x,iter);
    FILE *fp = fopen(name,"w");
    for(int ix =0;ix<numLines;ix++){
     for(int iy =0;iy<ngridcEext;iy++){
       int i = ix*ngridcEext + iy;
       fprintf(fp,"%d %d : %g %g\n",iy,ix,atmSF[i].re,atmSF[i].im);
     }//endfor
    }//endof
    fclose(fp);
#endif
    countEextFFT = 0;
#ifndef _DEBUG_INT_TRANS_FWD_
    launchFlag   = 1;
    if(registrationFlag==1){launchFlag=0; thisProxy(thisIndex.x, thisIndex.y).FFTEesBck();} 
#else   
    char name[100];
    sprintf(name,"partFFTGxGyZT%d.out.%d",rhoRsubplanes,thisIndex.x);
    FILE *fp = fopen(name,"w");
    for(int ix =0;ix<numLines;ix++){
     for(int iy =0;iy<ngridcEext;iy++){
       int i = ix*ngridcEext + iy;
       fprintf(fp,"%d %d : %g %g\n",iy,ix,atmSF[i].re,atmSF[i].im);
     }//endfor
    }//endof
    fclose(fp);
    UrhoGHartExtProxy[thisInstance.proxyOffset](0,0).exitForDebugging();
#endif
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Finish FFting to G-space  ::
//         2D)  atmSF(gx,gy,z) -> atmSF(gx,gy,gz)
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
  CkPrintf("GHart %d %d starting FFTBck at %d on %d\n",thisIndex.x, thisIndex.y,
	   iterAtmTyp,CkMyPe());
#endif

  int numPoints = rho_gs.numPoints;
  UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch()->doEextFFTRtoG_Gchare(atmSF,numFullEext,numPoints,
                                 numLines,rho_gs.numRuns,rho_gs.runs,ngridcEext,ind_x);

//============================================================================
// If the density has arrived, you can get the energy otherwise chill

  if(densityHere==1){getHartEextEes();}

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
// compute HartreeEextEes
//============================================================================
void CP_Rho_GHartExt::getHartEextEes(){
//============================================================================
// Output and Error check

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("GHart %d %d starting getHartEextEes at %d on %d\n",
           thisIndex.x,thisIndex.y,iterAtmTyp,CkMyPe());
#endif

  if(ees_eext_on==0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Yo dawg, ees eext is off. You can't getHartEextEes\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// Compute eext energy, hartree, total SF and VKS

 //----------------------------------------------------------
 // Local Variables
  eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();

  int myChareG       = thisIndex.x;
  int iterAtmTypFull = iterAtmTyp+atmTypoff;

  int ncoef    = rho_gs.nPacked;
  complex *vks = rho_gs.Vks;
  complex *rho = rho_gs.packedRho;
  int *k_x     = rho_gs.k_x;
  int *k_y     = rho_gs.k_y;
  int *k_z     = rho_gs.k_z;
  double *b_re = eesData->RhoGHartData[myChareG]->b_re;
  double *b_im = eesData->RhoGHartData[myChareG]->b_im;

 //----------------------------------------------------------
 // Initialize
  if(iterAtmTyp==1){  
     bzero(vks,ncoef*sizeof(complex));  // no getting around these zeros
     bzero(atmSFtot,ncoef*sizeof(complex));
     ehart_ret = 0.0;
     eext_ret  = 0.0;
     ewd_ret   = 0.0;
  }//endif

 //----------------------------------------------------------
 // Get the energy, vks, modifiy atmSF, contribute to total SF
#if CMK_TRACE_ENABLED
  double StartTime=CmiWallTimer();
#endif
  double *perdCorr = rho_gs.perdCorr;
  CPLOCAL::eesHartEextGchare(ncoef,iterAtmTypFull,rho,vks,atmSF,atmSFtot,
                             b_re,b_im,&ehart_ret,&eext_ret,k_x,k_y,k_z,perdCorr,myChareG,config.nfreq_cplocal_eeshart);
#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesHartExcG_, StartTime, CmiWallTimer());    
#endif

#ifdef _CP_GHART_VERBOSE_
  CkPrintf("Ghart %d %d Here in hartees: EHart %.10g Eext %.10g at %d on %d\n",
	   thisIndex.x,thisIndex.y,ehart_ret,eext_ret,iterAtmTyp,CkMyPe());
#endif

//============================================================================
// If you have SFtot get the ewald energy  : A reduction is required when atmTyp parallel

  if(iterAtmTyp==natmTyp && nchareHartAtmT==1){
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif
    double *perdCorr = rho_gs.perdCorr;
    CPLOCAL::eesEwaldGchare(ncoef,atmSFtot,b_re,b_im,&ewd_ret,k_x,k_y,k_z,perdCorr,myChareG,config.nfreq_cplocal_eesewald);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(eesEwaldG_, StartTime, CmiWallTimer());    
#endif
#ifdef _CP_GHART_VERBOSE_
    CkPrintf("Ghart %d %d iter %d Here in hartees: EwaldRecip : %.10g on\n",
              thisIndex.x,thisIndex.y,iterAtmTyp,ewd_ret,CkMyPe());
#endif
  }//endif
     

//============================================================================
// Blast out the energies when done with atom type stuff : index=0 may have to wait

  if(iterAtmTyp==natmTyp){
#ifdef _CP_GHART_VERBOSE_
    CkPrintf("Ghart %d reduces energies at %d on %d\n",iterAtmTyp,CkMyPe());
#endif
    if(thisIndex.y!=0 || nchareHartAtmT==1){
        double e[3];
	e[0] = ehart_ret;
	e[1] = eext_ret;
	e[2] = ewd_ret;
	contribute(3*sizeof(double),e,CkReduction::sum_double);
    }//endif
  }//endif

//============================================================================
// Perform the back FFT SF and SFtot to get atm forces and VKS to get e-forces

 //-----------------------------------------------------------------
 // Ewald needs contribs from all SF  : chare index 0 is large and in charge
  if(iterAtmTyp==natmTyp){
    if(nchareHartAtmT==1){
      FFTEesFwd(1);
    }else{
      UrhoGHartExtProxy[thisInstance.proxyOffset](thisIndex.x,0).acceptAtmSFTot(ncoef,atmSFtot);  
    }//endif
  }//endif

 //-----------------------------------------------------------------
 // Vks needs contribs from all SF : chare index 1 is large and in charge
  if(iterAtmTyp==natmTyp){
    if(nchareHartAtmT==1)  { 
      FFTVks(); 
    }else{
      UrhoGHartExtProxy[thisInstance.proxyOffset](thisIndex.x,1).acceptVks(ncoef,vks);  
    }//endif
  }//endif

 //-----------------------------------------------------------------
 // we always geneterate an SF : we have everthing for this atmtype
 //                            : flow of control says ``this guy goes last''
 //                            : This guy controls exit condition.
  FFTEesFwd(0);                // DON'T MOVE HIM

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
  CkPrintf("Ghart %d %d Here in FFT to r: %.10g %.10g : with %d at %d %d on %d\n",
            thisIndex.x,thisIndex.y,ehart_ret,eext_ret,flag,iterAtmTyp,natmTyp,CkMyPe());
#endif

 //--------------------------------------------
 // Get the data and some scratch
  FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
  fftcache->getCacheMem("CP_Rho_GHartExt::FFTEesFwd");
  complex *data_out  = fftcache->tmpData;
  complex *data;     if(flag==0){data=atmSF;}else{data=atmSFtot;}

 //--------------------------------------------
 // FFT the data : result goes into scratch
  fftcache->doEextFFTGtoR_Gchare(data,data_out,numFullEext,rho_gs.numPoints,numLines,
 			         rho_gs.numRuns,rho_gs.runs,ngridcEext,ind_x);

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
  CkPrintf("Ghart %d %d Here in send to Rhart with opt %d at %d %d on %d\n",
         thisIndex.x,thisIndex.y,flag,iterAtmTyp,natmTyp,CkMyPe());
#endif

  CPcharmParaInfo *sim = CPcharmParaInfo::get();
  FFTcache *fftcache = UfftCacheProxy[thisInstance.proxyOffset].ckLocalBranch();  
  complex *senddata  = fftcache->tmpData;
  int ix             = thisIndex.x;
  int numLinesNow    = numLines;
  int **index_pack;
  int *nline_send;
  if(rhoRsubplanes>1){
    nline_send  = sim->nline_send_eext_y[ix];
    index_pack  = sim->index_tran_pack_eext_y[ix];
  }//endif

//============================================================================
// start commlib
#ifdef USE_COMLIB
  if(rhoRsubplanes==1){
    switch(flag){
#ifdef OLD_COMMLIB
    case 0 : if(config.useGHartInsRHart) commGHartRHartIns0.beginIteration();break;
    case 1 : if(config.useGHartInsRHart) commGHartRHartIns1.beginIteration();break;
#else
    case 0 : if(config.useGHartInsRHart) ComlibBegin(rhoRHartProxy_com0,1); break;
    case 1 : if(config.useGHartInsRHart) ComlibBegin(rhoRHartProxy_com1,1);break;
#endif
    }//endif
  }//endif
#endif

//============================================================================
// Send the message : 1 pt from each line to each chareR

  int priority =config.rhorHartpriority + thisIndex.x*10;
  int prioritybits = 8*sizeof(int);

  for(int z=0; z < ngridcEext; z++) {
  for(int s=0;s<rhoRsubplanes;s++){

    if(rhoRsubplanes>1){numLinesNow = nline_send[s];}
    if(numLinesNow >0){
	RhoRHartMsg *msg = new (numLinesNow, prioritybits) RhoRHartMsg;
	msg->size        = numLinesNow;     // number of z-lines in this batch
	msg->senderIndex = thisIndex.x;  // line batch index
	msg->iopt        = flag;         // iopt
	msg->iter        = iterAtmTyp;

	if(config.prioEextFFTMsg){
	  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
	  *(int*)CkPriorityPtr(msg) = priority;
	}//endif

	// beam out all points with same z to chare array index z
	complex *data = msg->data;
	if(rhoRsubplanes==1){
	  for (int i=0,j=z; i<numLines; i++,j+=ngridcEext){data[i] = senddata[j];}
	}else{
	  int *index_packs=index_pack[s];    
	  for(int i=0; i< numLinesNow; i++){
	    data[i] = senddata[(z+index_packs[i])];
	  }//endif
	}//endif

	if(rhoRsubplanes==1){
	  switch(iopt){
 	   case 0 : rhoRHartProxy_com0(z,0,thisIndex.y).recvAtmForcFromRhoGHart(msg);break;
	   case 1 : rhoRHartProxy_com1(z,0,thisIndex.y).recvAtmForcFromRhoGHart(msg);break;
	  }//end switch
	}else{
	  switch(iopt){
	   case 0 : UrhoRHartExtProxy[thisInstance.proxyOffset](z,s,thisIndex.y).recvAtmForcFromRhoGHart(msg);break;
	   case 1 : UrhoRHartExtProxy[thisInstance.proxyOffset](z,s,thisIndex.y).recvAtmForcFromRhoGHart(msg);break;
	  }//end switch
	}//endif

    }// endif : we need to send to this subplane
  }// endfor : subplane loop
#ifdef CMK_BLUEGENEL
  if(z%8==0){CmiNetworkProgress();}
#endif
  }//endfor : z-plane parallel loop

//============================================================================
// end commlib

#ifdef USE_COMLIB
  if(rhoRsubplanes==1){
   switch(flag){
#ifdef OLD_COMMLIB
     case 0 : if(config.useGHartInsRHart) commGHartRHartIns0.endIteration();break;
     case 1 : if(config.useGHartInsRHart) commGHartRHartIns1.endIteration();break;
#else
    case 0 : if(config.useGHartInsRHart) ComlibEnd(rhoRHartProxy_com0,1); break;
    case 1 : if(config.useGHartInsRHart) ComlibEnd(rhoRHartProxy_com1,1);break;
#endif

   }//end switch
  }//endif
#endif

//============================================================================
// We are done when when have sent out all SFs and the Ewald total SF (index=0)

  int nsendExpect=natmTyp;
  if(thisIndex.y==0){nsendExpect++;} // chare 0 sends out SFtot, too.
    
  nsendAtmTyp ++;
  if(nsendAtmTyp==nsendExpect){
    nsendAtmTyp = 0;
    iterAtmTyp  = 0;
    atmSFHere   = 0;
    densityHere = 0; 
  }//endif

  fftcache->freeCacheMem("CP_Rho_GHartExt::FFTEesFwd");

//============================================================================
   }//end routine
//============================================================================

//============================================================================
// Collect the SF from all the atm type chares on chare 0
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::acceptAtmSFTot(int size, complex *inSF){
//============================================================================
// Recv the contribs.

  CkAssert(thisIndex.y==0);

  if(countAtmSFtot==0){bzero(atmSFtotRecv,size*sizeof(complex));}
  countAtmSFtot++;

  for(int i=0;i<size;i++){atmSFtotRecv[i]+=inSF[i];}

//============================================================================
// Once you have it all, compute the energy, contribute, fft back, send to Rhart

  if(countAtmSFtot==nchareHartAtmT){
    countAtmSFtot=0;

   //---------------------------------------------------------
   // Compute ewald energy, modify the atmSFtot appropriately
    int myChareG = thisIndex.x;
    int ncoef    = rho_gs.nPacked;
    int *k_x     = rho_gs.k_x;
    int *k_y     = rho_gs.k_y;
    int *k_z     = rho_gs.k_z;
    eesCache *eesData  = UeesCacheProxy[thisInstance.proxyOffset].ckLocalBranch ();
    double *b_re = eesData->RhoGHartData[myChareG]->b_re;
    double *b_im = eesData->RhoGHartData[myChareG]->b_im;
#if CMK_TRACE_ENABLED
    double StartTime=CmiWallTimer();
#endif
    double *perdCorr = rho_gs.perdCorr;
    CPLOCAL::eesEwaldGchare(ncoef,atmSFtotRecv,b_re,b_im,&ewd_ret,k_x,k_y,k_z,perdCorr,myChareG,config.nfreq_cplocal_eesewald);
#if CMK_TRACE_ENABLED
    traceUserBracketEvent(eesEwaldG_, StartTime, CmiWallTimer());    
#endif

   //---------------------------------------------------------
   // Contribute your energies now that you have them all
    double e[3];
    e[0] = ehart_ret;
    e[1] = eext_ret;
    e[2] = ewd_ret;
    contribute(3*sizeof(double),e,CkReduction::sum_double);

   //---------------------------------------------------------
   // FFT back, which generates a send back to RHart
    complex *junk=atmSFtot;
    atmSFtot=atmSFtotRecv;  // avoid a CmiMemcpy : tricky, tricky
    FFTEesFwd(1);
    atmSFtot=junk;

  }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
// Collect the VKS from all the atm type chares on chare 1
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GHartExt::acceptVks(int size, complex * inVks){
//============================================================================
// Recv the contributions

  CkAssert(thisIndex.y==1);

  if(countVksTot==0){bzero(VksRecv,size*sizeof(complex));}
  countVksTot++;

  for(int i=0;i<size;i++){VksRecv[i] += inVks[i];}

//============================================================================
// When all the guys have reported, you can do the fft and then send vks to RhoReal

  if(countVksTot==nchareHartAtmT){
     countVksTot=0;
     complex *junk=rho_gs.Vks;
     rho_gs.Vks= VksRecv;  // save a CmiMemcpy : tricky, tricky
     FFTVks(); 
     rho_gs.Vks=junk;
  }//endif

//============================================================================
   }//end routine
//============================================================================


//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Glenn's special exit 
//============================================================================
void CP_Rho_GHartExt::exitForDebugging(){
//============================================================================

  CPcharmParaInfo *sim      = CPcharmParaInfo::get();
  int nchareG               = sim->nchareRhoGEext;

  CountDebug++;  
  if(CountDebug==nchareG*nchareHartAtmT){
    CountDebug=0;
    CkPrintf("I am in the exitfordebuging rhoghartext puppy. Bye-bye\n");
    CkExit();
  }//endif

//============================================================================
  }//end routine
//============================================================================
/*@}*/
