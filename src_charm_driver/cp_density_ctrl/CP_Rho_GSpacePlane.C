//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_Rho_GSpacePlane.C
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

#include "charm++.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "../../include/debug_flags.h"
#include "util.h"
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"

#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern Config config;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern int nstates;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern ComlibInstanceHandle commInstance;
extern CProxy_ComlibManager mgrProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_FFTcache fftCacheProxy;
extern ComlibInstanceHandle mssInstance;

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::CP_Rho_GSpacePlane(int xdim, size2d sizeYZ, 
	       int numRealSpace, int numRealSpaceDensity, bool _useCommlib,
	       ComlibInstanceHandle _fftcommInstance) 
//============================================================================
    {//begin routine
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
    CkPrintf("[%d %d] Rho GS constructor\n",thisIndex.x,thisIndex.y);
#endif

//============================================================================

    CkAssert(numRealSpaceDensity == 1);
    CkAssert(config.rhoGPPC == 1);
    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
    CkVec <RunDescriptor> *sortedRunDescriptors = sim->RhosortedRunDescriptors;

    rho_gs.sizeX    = xdim;
    rho_gs.sizeY    = sizeYZ[0];
    rho_gs.sizeZ    = sizeYZ[1];
    rho_gs.xdim     = rho_gs.sizeX;
    rho_gs.ydim     = rho_gs.sizeY;
    rho_gs.zdim     = 1;

    count = 0;
    for(int i=1;i<=3;i++){countWhiteByrd[i]=0;}
    doneWhiteByrd = 0;

    int x = thisIndex.x;
    rho_gs.numRuns  = sortedRunDescriptors[x].size();
    rho_gs.numLines = sortedRunDescriptors[x].size()/2;
    rho_gs.numFull  = (rho_gs.numLines)*rho_gs.sizeZ;
    rho_gs.size     = rho_gs.numFull;
    rho_gs.runs     = new RunDescriptor[rho_gs.numRuns];
    rho_gs.numPoints = 0;
    for (int r = 0; r < rho_gs.numRuns; r++) {
      rho_gs.numPoints += sortedRunDescriptors[x][r].length;
      rho_gs.runs[r]    = sortedRunDescriptors[x][r];
    }//endfor
    rho_gs.setKVectors(&nPacked);
    
    rho_gs.Rho       = (complex *)fftw_malloc(rho_gs.numFull*sizeof(complex));
    rho_gs.divRhoX   = (complex *)fftw_malloc(rho_gs.numFull*sizeof(complex));
    rho_gs.divRhoY   = (complex *)fftw_malloc(rho_gs.numFull*sizeof(complex));
    rho_gs.divRhoZ   = (complex *)fftw_malloc(rho_gs.numFull*sizeof(complex));
    rho_gs.packedRho = (complex *)fftw_malloc(nPacked*sizeof(complex));
    rho_gs.packedVks = (complex *)fftw_malloc(nPacked*sizeof(complex));

    setMigratable(false);
    rhoRealProxy_com = rhoRealProxy;
    if(config.useCommlib){
       ComlibAssociateProxy(&commInstance,rhoRealProxy_com);
    }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::~CP_Rho_GSpacePlane()
{
}
//============================================================================


//============================================================================
//  RhoReal sends rho(gx,gy,z) here such that it is now decomposed 
//  with lines of constant gx,gy in funny order to load balance.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptData(RhoGSFFTMsg *msg) {
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("RGS [%d %d] receives message size %d offset %d numlines %d\n",
             thisIndex.x,thisIndex.y,msg->size,msg->offset,rho_gs.numLines);
#endif

//============================================================================

  int size             = msg->size;
  int offset           = msg->offset;
  complex *partlyIFFTd = msg->data;

  int numLines         = rho_gs.numLines;
  int sizeZ            = rho_gs.sizeZ;
  int expandedDataSize = numLines*sizeZ;

  if(size!=numLines){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("size %d != numLines %d\n",size,numLines);
    CkPrintf("RGS [%d %d] receives message size %d offset %d numlines %d\n",
             thisIndex.x,thisIndex.y,msg->size,msg->offset,rho_gs.numLines);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================

  if(rho_gs.Rho==NULL) {
    rho_gs.Rho = (complex *)fftw_malloc(expandedDataSize*sizeof(complex));
    memset(rho_gs.Rho,0,sizeof(complex)*expandedDataSize);
  }//endif

  complex *chunk = rho_gs.Rho;
  for(int i=0,j=offset; i< size; i++,j+=sizeZ){chunk[j] = partlyIFFTd[i];}
  delete msg;

  count++;
  if(count==sizeZ){
    count=0;
    acceptData();
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptData() { 
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
    CkPrintf("Data has arrive in rhogs : peforming FFT %d %d\n",
              thisIndex.x,thisIndex.y);
#endif

//============================================================================
// I) Finish ffting  the density into g-space :
//      a)FFT b) Compress to sphere c) store in packedRho

    int ioptRho=0;
    rho_gs.doBwFFTRtoG(ioptRho);

#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif

//============================================================================
// II) Debug output

#ifdef _CP_DEBUG_RHOG_RHOG_
     char myFileName[100];
     sprintf(myFileName, "Rho_Gspace_%d%d.out", thisIndex.x,thisIndex.y);
     FILE *fp = fopen(myFileName,"w");
       for (int i = 0; i < rho_gs.numPoints; i++){ 
              fprintf(fp," %d %d %d : %g %g\n",
                 rho_gs.k_x[i],rho_gs.k_y[i],rho_gs.k_z[i],
                 rho_gs.packedRho[i].re,rho_gs.packedRho[i].im);
       }//endfor
     fclose(fp);
#endif

//============================================================================
// II) First compute hart+Exc and vks(g) using rho(g) then start gradient comp
    
    HartExcVksG();

#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();
#endif

    divRhoVksGspace();

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::HartExcVksG() { 
//============================================================================

   AtomsGrp *ag = atomsGrpProxy.ckLocalBranch(); // find me the local copy
   complex *rho      = rho_gs.packedRho;
   complex *vks      = rho_gs.packedVks;
   int nsize         = rho_gs.numPoints;
   int natm          = ag->natm;
   Atom *atoms       = ag->atoms;
   double *ehart_ret = &(rho_gs.ehart_ret);
   double *eext_ret  = &(rho_gs.eext_ret);
   double *ewd_ret   = &(rho_gs.ewd_ret);
   int *k_x          = rho_gs.k_x;
   int *k_y          = rho_gs.k_y;
   int *k_z          = rho_gs.k_z;

   CPLOCAL::CP_hart_eext_calc(nsize,rho,natm,atoms,vks,
                              ehart_ret,eext_ret,ewd_ret,k_x,k_y,k_z,
                              thisIndex.x);
   double e[3];
   e[0] = rho_gs.ehart_ret;
   e[1] = rho_gs.eext_ret;
   e[2] = rho_gs.ewd_ret;
   contribute(3 * sizeof(double),e,CkReduction::sum_double);

#ifdef _CP_DEBUG_RHOG_VKSA_
     char myFileName[100];
     sprintf(myFileName, "Vks_Gspace_%d%d.out", thisIndex.x,thisIndex.y);
     FILE *fp = fopen(myFileName,"w");
       for (int i = 0; i < rho_gs.numPoints; i++){ 
              fprintf(fp," %d %d %d : %g %g\n",
                 rho_gs.k_x[i],rho_gs.k_y[i],rho_gs.k_z[i],
                 rho_gs.packedVks[i].re,rho_gs.packedVks[i].im);
       }//endfor
     fclose(fp);
#endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// invoke fft subroutine : doFwFFT() unpack do an FFT 1D along z rho(gx,gy,gz)
// parallelized in full z-lines with constant (gx,gy)
// invoke a transpose back to r-space
//============================================================================

void CP_Rho_GSpacePlane::divRhoVksGspace() { 

//============================================================================

  int ioptx  = 1; 
  int iopty  = 2; 
  int ioptz  = 3; 
  int ioptvks= 0; 
  double tpi,*hmati;
  CPXCFNCTS::CP_fetch_hmati(&hmati,&tpi);

//--------------------------------------------------------------------------
// 0) Fill divRho(i) = g(i)*packedrho for i = x,y,z

  rho_gs.divRhoGdot(hmati,tpi); 
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//--------------------------------------------------------------------------
// I)   fft_gz(divRhoX), launch transpose to R-space

  rho_gs.doFwFFTGtoR(ioptx); 
  RhoGSendRhoR(ioptx);
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//--------------------------------------------------------------------------
// II)  fft_gz(divRhoY), launch transpose to R-space

  rho_gs.doFwFFTGtoR(iopty); 
  RhoGSendRhoR(iopty);
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//--------------------------------------------------------------------------
// III)  fft_gz(divRhoZ), launch transpose  to R-space

  rho_gs.doFwFFTGtoR(ioptz); 
  RhoGSendRhoR(ioptz);
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//--------------------------------------------------------------------------
// IV)  fft_gz(vks), launch transpose  to R-space

  rho_gs.expandRhoGSpace(rho_gs.Rho,rho_gs.packedVks);
  rho_gs.doFwFFTGtoR(ioptvks); 
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif
  RhoGSendRhoR(ioptvks);

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::RhoGSendRhoR(int iopt) { 
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("Communicating data from RhoG to RhoR : %d %d %d\n",
	   thisIndex.x,thisIndex.y,iopt);
#endif

//============================================================================

  int numLines         = rho_gs.numLines; // same amount of data to each Rspace chare
  int sizeZ            = rho_gs.sizeZ;

  complex *ffttempdata;
  switch(iopt){
    case 0 : ffttempdata = rho_gs.Rho;     break; //its vks
    case 1 : ffttempdata = rho_gs.divRhoX; break;
    case 2 : ffttempdata = rho_gs.divRhoY; break;
    case 3 : ffttempdata = rho_gs.divRhoZ; break;
    case 4 : ffttempdata = rho_gs.Rho;     break;  //its whitebyrd
  }//end switch

//============================================================================
// Do a Comlib Dance

//  if (config.useCommlib){commInstance.beginIteration();}

//============================================================================
  for(int z=0; z < sizeZ; z++) {

    RhoRSFFTMsg *msg = new (numLines,8*sizeof(int)) RhoRSFFTMsg;
    msg->size        = numLines;     // number of z-lines in this batch
    msg->senderIndex = thisIndex.x;  // line batch index
    msg->iopt        = iopt;         // iopt
    
    if(config.prioFFTMsg){
       CkSetQueueing(msg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(msg) = config.rhorpriority + thisIndex.x*numLines
                                                       + thisIndex.y;
    }//endif

    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = ffttempdata[j];}

    switch(iopt){
      /*
      case 0 : rhoRealProxy_com(z,0).acceptGradRhoVks(msg); break;
      case 1 : rhoRealProxy_com(z,0).acceptGradRhoVks(msg); break;
      case 2 : rhoRealProxy_com(z,0).acceptGradRhoVks(msg); break;
      case 3 : rhoRealProxy_com(z,0).acceptGradRhoVks(msg); break;
      case 4 : rhoRealProxy_com(z,0).acceptWhiteByrd(msg); break;
      */
      case 0 : rhoRealProxy(z,0).acceptGradRhoVks(msg); break;
      case 1 : rhoRealProxy(z,0).acceptGradRhoVks(msg); break;
      case 2 : rhoRealProxy(z,0).acceptGradRhoVks(msg); break;
      case 3 : rhoRealProxy(z,0).acceptGradRhoVks(msg); break;
      case 4 : rhoRealProxy(z,0).acceptWhiteByrd(msg); break;

    }//end switch

  }//endfor

//============================================================================
// Complete the commlib dance
    
//  if (config.useCommlib){commInstance.endIteration();}

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//  RhoReal sends rho(gx,gy,z) here such that it is now decomposed 
//  with lines of constant gx,gy in funny order to load balance.
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptWhiteByrd(RhoGSFFTMsg *msg) {
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("RGS [%d %d] receives white-byrd message size %d offset %d numlines %d\n",
             thisIndex.x,thisIndex.y,msg->size,msg->offset,rho_gs.numLines);
#endif

//============================================================================

  int size             = msg->size;
  int offset           = msg->offset;
  int iopt             = msg->iopt;
  complex *partlyIFFTd = msg->data;

  int numLines         = rho_gs.numLines;
  int sizeZ            = rho_gs.sizeZ;
  int expandedDataSize = numLines*sizeZ;

//============================================================================
// Error checking

  countWhiteByrd[iopt]++;

  if(size!=numLines){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("size %d != numLines %d\n",size,numLines);
    CkPrintf("RGS [%d %d] receives message size %d offset %d numlines %d opt %d\n",
             thisIndex.x,thisIndex.y,size,offset,numLines,iopt);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// unpack the message

  complex *chunk;
  switch(iopt){
    case 1 : chunk = rho_gs.divRhoX; break;
    case 2 : chunk = rho_gs.divRhoY; break;
    case 3 : chunk = rho_gs.divRhoZ; break;
  }//end switch

  for(int i=0,j=offset; i< size; i++,j+=sizeZ){chunk[j]=partlyIFFTd[i];}
  delete msg;

//============================================================================
// If all chares for a gradient report, perform final FFT. 
// If all gradient are complete, construct whitebyrd.

  if(countWhiteByrd[iopt]==sizeZ){
    countWhiteByrd[iopt]=0;
    doneWhiteByrd++;
    rho_gs.doBwFFTRtoG(iopt);
  }//endif

  if(doneWhiteByrd==3){
    doneWhiteByrd=0;
#ifdef CMK_VERSION_BLUEGENE
    CmiNetworkProgress();    
#endif
    acceptWhiteByrd();
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptWhiteByrd() {
//============================================================================

  int ioptFFTWhite  = 0; 
  int ioptSendWhite = 4; 
  double tpi,*hmati;
  CPXCFNCTS::CP_fetch_hmati(&hmati,&tpi);

  rho_gs.createWhiteByrd(hmati,tpi);
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();    
#endif

  rho_gs.doFwFFTGtoR(ioptFFTWhite); // stored in Rho
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();    
#endif

  RhoGSendRhoR(ioptSendWhite);

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::ResumeFromSync(){

  if (config.useCommlib) {
    ComlibResetProxy(&rhoRealProxy_com);
  }//endif

}
//============================================================================



