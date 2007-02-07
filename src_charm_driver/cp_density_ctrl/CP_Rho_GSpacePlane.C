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
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"

#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern Config config;
extern CProxy_CP_Rho_RealSpacePlane rhoRealProxy;
extern CProxy_CP_Rho_GHartExt rhoGHartExtProxy;
extern int nstates;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CProxy_AtomsGrp atomsGrpProxy;
extern ComlibInstanceHandle commGInstance0;
extern ComlibInstanceHandle commGInstance1;
extern ComlibInstanceHandle commGInstance2;
extern ComlibInstanceHandle commGInstance3;
extern ComlibInstanceHandle commGByrdInstance;

extern CProxy_ComlibManager mgrProxy;
extern CProxy_CP_Rho_GSpacePlane rhoGProxy;
extern CProxy_FFTcache fftCacheProxy;
//extern ComlibInstanceHandle mssInstance;

//#define _CP_DEBUG_RHOG_VERBOSE_
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::CP_Rho_GSpacePlane(int sizeX, size2d sizeYZ, 
	       int numRealSpace, int numRealSpaceDensity, bool _useCommlib)
//============================================================================
    {//begin routine
//============================================================================
#ifdef _CP_DEBUG_RHOG_VERBOSE_
    CkPrintf("[%d %d] Rho GS constructor\n",thisIndex.x,thisIndex.y);
#endif
//============================================================================
// Set counters local variables

    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
    cp_grad_corr_on = sim->cp_grad_corr_on;

    iplane_ind    = thisIndex.x;
    rhoRsubplanes = config.rhoRsubplanes;
    count = 0;
    doneWhiteByrd = 0;
    for(int i=1;i<=3;i++){countWhiteByrd[i]=0;}
    countDebug=0;

//============================================================================
// Deal with the run descriptors then malloc

   //------------------------------------------------------------------------
   // Run descriptors
    CkVec <RunDescriptor> *sortedRunDescriptors = sim->RhosortedRunDescriptors;

    rho_gs.sizeX    = sizeX;
    rho_gs.sizeY    = sizeYZ[0];
    rho_gs.sizeZ    = sizeYZ[1];

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

   //------------------------------------------------------------------------
   // Malloc and set your kvectors
    rho_gs.setKVectors(&nPacked); // consistency with numPoints checked within
    rho_gs.nPacked   = nPacked;

   //------------------------------------------------------------------------
   //  Malloc density (packed Rho) and a master scratch array (Rho)
    int numFull      = rho_gs.numFull;
    // Full size Needed to receive from realSpace whenever it feels like sending
    rho_gs.divRhoX   = (complex *)fftw_malloc(numFull*sizeof(complex));
    rho_gs.divRhoY   = (complex *)fftw_malloc(numFull*sizeof(complex));
    rho_gs.divRhoZ   = (complex *)fftw_malloc(numFull*sizeof(complex));

    rho_gs.Rho       = NULL;  // clean out my memory
    rho_gs.packedVks = NULL;
    rho_gs.packedRho = NULL;

//============================================================================
// Decompose your g-space further to create/define GHartEext Chare

    rhoGHelpers       = config.rhoGHelpers;
    int num_line_tot  = rho_gs.numLines;
    numSplit          = new int [rhoGHelpers];
    istrtSplit        = new int [rhoGHelpers];
    iendSplit         = new int [rhoGHelpers];

    x = thisIndex.x;
    int istrt_pts = 0;
    for(int i=0;i<rhoGHelpers;i++){

      int istrt_line,iend_line,num_line;
      getSplitDecomp(&istrt_line,&iend_line,&num_line,num_line_tot,rhoGHelpers,i);

      int num_pts=0;
      for (int r=(2*istrt_line);r<(2*iend_line);r++){
        num_pts += sortedRunDescriptors[x][r].length;
      }//endif

      istrtSplit[i] = istrt_pts;
      numSplit[i]   = num_pts;     
      iendSplit[i]  = istrt_pts+num_pts;

      istrt_pts    += num_pts;
    }//endfor

//============================================================================
// Set up proxies and set migration options

    rhoRealProxy0_com = rhoRealProxy;
    rhoRealProxy1_com = rhoRealProxy;
    rhoRealProxy2_com = rhoRealProxy;
    rhoRealProxy3_com = rhoRealProxy;
    rhoRealProxyByrd_com = rhoRealProxy;
    if(config.useGIns0RhoRP)
       ComlibAssociateProxy(&commGInstance0,rhoRealProxy0_com);
    if(config.useGIns1RhoRP)
       ComlibAssociateProxy(&commGInstance1,rhoRealProxy1_com);
    if(config.useGIns2RhoRP)
       ComlibAssociateProxy(&commGInstance2,rhoRealProxy2_com);
    if(config.useGIns3RhoRP)
       ComlibAssociateProxy(&commGInstance3,rhoRealProxy3_com);
    if(config.useGByrdInsRhoRBP)
       ComlibAssociateProxy(&commGByrdInstance,rhoRealProxyByrd_com);

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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_Rho_GSpacePlane::~CP_Rho_GSpacePlane(){
    if(numSplit!=NULL)
	delete [] numSplit;
    if(istrtSplit!=NULL)
	delete [] istrtSplit;
    if(iendSplit != NULL)
	delete [] iendSplit;
}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);
  p|iplane_ind;
  p|nPacked;
  p|count;
  PUParray(p,countWhiteByrd,4);
  p|doneWhiteByrd;
  p|rhoGHelpers;
  p|rhoRsubplanes;
  rho_gs.pup(p);  // I pup my data class
  if(p.isUnpacking()){
      numSplit   = new int[rhoGHelpers];
      istrtSplit = new int[rhoGHelpers];
      iendSplit  = new int[rhoGHelpers];
  }//endif
  PUParray(p,numSplit,rhoGHelpers);
  PUParray(p,istrtSplit,rhoGHelpers);
  PUParray(p,iendSplit,rhoGHelpers);
  if(p.isUnpacking()){ 
      //pupping of comlib proxies unreliable
      rhoRealProxy0_com = rhoRealProxy;
      rhoRealProxy1_com = rhoRealProxy;
      rhoRealProxy2_com = rhoRealProxy;
      rhoRealProxy3_com = rhoRealProxy;
      rhoRealProxyByrd_com = rhoRealProxy;
      if(config.useGIns0RhoRP)
	ComlibAssociateProxy(&commGInstance0,rhoRealProxy0_com);
      if(config.useGIns1RhoRP)
	ComlibAssociateProxy(&commGInstance1,rhoRealProxy1_com);
      if(config.useGIns2RhoRP)
	ComlibAssociateProxy(&commGInstance2,rhoRealProxy2_com);
      if(config.useGIns3RhoRP)
	ComlibAssociateProxy(&commGInstance3,rhoRealProxy3_com);
      if(config.useGByrdInsRhoRBP)
	ComlibAssociateProxy(&commGByrdInstance,rhoRealProxyByrd_com);
   }//endif

}//end routine
//============================================================================


//============================================================================
//  RhoReal sends rho(gx,gy,z) here such that it is now decomposed 
//  with lines of constant gx,gy in funny order to load balance.
//  Nothing can return to this chare UNTIL this message is fully received
//  and processed. Thus, we can receive into divRhox and use cache to post 
//  process
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptRhoData(RhoGSFFTMsg *msg) {
//============================================================================
#ifdef _CP_DEBUG_RHOG_VERBOSE_
  CkPrintf("RGS [%d %d] receives message size %d offset %d numlines %d\n",
             thisIndex.x,thisIndex.y,msg->size,msg->offset,rho_gs.numLines);
#endif
//============================================================================
// Set local pointers, constants,

  int size             = msg->size;
  int offset           = msg->offset;   // z index
  int isub             = msg->offsetGx; // subplane index
  complex *partlyIFFTd = msg->data;

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int ix               = thisIndex.x;   // chare array index
  int sizeZ            = rho_gs.sizeZ;
  int numLines         = rho_gs.numLines;
  int ***index_pack;
  if(rhoRsubplanes>1){
    int **nline_send = sim->nline_send_rho_y;
    index_pack       = sim->index_tran_pack_rho_y;
    numLines         = nline_send[ix][isub];
  }//endif

//============================================================================
// Check for errors

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
  }
#endif
  if(size!=numLines){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("size %d != numLines %d\n",size,numLines);
    CkPrintf("RGS [%d %d] receives message size %d offset %d numlines %d\n",
             thisIndex.x,thisIndex.y,msg->size,msg->offset,rho_gs.numLines);
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

//============================================================================
// upack the message : No need to zero, every value is set

  complex *chunk = rho_gs.divRhoX;

  if(rhoRsubplanes==1){
    for(int i=0,j=offset; i< size; i++,j+=sizeZ){chunk[j] = partlyIFFTd[i];}
  }else{
    for(int i=0; i< size; i++){
      chunk[(offset+index_pack[ix][isub][i])] = partlyIFFTd[i];
    }//endif
  }//endif

  delete msg;

  if(count==0){count_stuff=0;}
  count_stuff += size;

  count++;
  if(count==sizeZ*rhoRsubplanes){
    if(count_stuff!=rho_gs.numLines*sizeZ){
      CkPrintf("Not enough stuff %d %d on %d\n",count_stuff,rho_gs.numLines,thisIndex.x);
      CkExit();
    }//endif
    count=0;
#ifndef _DEBUG_INT_TRANS_FWD
    acceptRhoData();
#else   
    char name[100];
    sprintf(name,"partFFTGxGyZT%d.out.%d",rhoRsubplanes,thisIndex.x);
    FILE *fp = fopen(name,"w");
    numLines = rho_gs.numLines;
    for(int ix =0;ix<numLines;ix++){
     for(int iy =0;iy<sizeZ;iy++){
       int i = ix*sizeZ + iy;
       fprintf(fp,"%d %d : %g %g\n",iy,ix,chunk[i].re,chunk[i].im);
     }//endfor
    }//endof
    fclose(fp);
    rhoGProxy(0,0).exitForDebugging();
#endif
  }//endif

//============================================================================
  }//end routine
//============================================================================


//============================================================================
// Complete the FFT to give use Rho(gx,gy,gz)
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::acceptRhoData() { 
//============================================================================

#ifdef _CP_DEBUG_RHOG_VERBOSE_
    CkPrintf("Data has arrive in rhogs : peforming FFT %d %d\n",
              thisIndex.x,thisIndex.y);
#endif

//============================================================================
// I) Finish ffting  the density into g-space :
//      a)FFT b) Compress to sphere and store in data_out

    int ioptRho=0;
#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif
    FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
    fftcache->getCacheMem("CP_Rho_GSpacePlane::acceptRhoData");
    complex *data_out  = fftcache->tmpData;
    complex *data_in   = rho_gs.divRhoX;
    fftcache->doRhoFFTRtoG_Gchare(data_in,data_out,
                                  rho_gs.numFull,rho_gs.numPoints,
                                  rho_gs.numLines,rho_gs.numRuns,rho_gs.runs, 
                                  rho_gs.sizeZ,1,iplane_ind);
#ifndef CMK_OPTIMIZE    
    traceUserBracketEvent(BwFFTRtoG_, StartTime, CmiWallTimer());    
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
// II) Communicate rho(g) to RHoGHartExt to compute eext and hart part of vks
//     or print out that you are taking the day off

#ifndef _CP_DEBUG_HARTEEXT_OFF_  // hartree is cooking

     for(int i=0;i<rhoGHelpers;i++){
  
       RhoGHartMsg *msg = new (numSplit[i],8*sizeof(int)) RhoGHartMsg;
       msg->size        = numSplit[i];     
       msg->senderIndex = thisIndex.x;
       memcpy(msg->data,&data_out[istrtSplit[i]],numSplit[i]*sizeof(complex));
    
      if(config.prioFFTMsg){
         CkSetQueueing(msg, CK_QUEUEING_IFIFO);
         *(int*)CkPriorityPtr(msg) = config.rhorpriority + thisIndex.x+ thisIndex.y;
       }//endif
    
       int index = thisIndex.x*rhoGHelpers + i;
       rhoGHartExtProxy(index,thisIndex.y).acceptData(msg);

#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress(); 
#endif
     }//endfor

#else //hartree is sitting this one out

    if(thisIndex.x==0 && thisIndex.y==0){
       CkPrintf("EHART       = OFF FOR DEBUGGING\n");
       CkPrintf("EExt        = OFF FOR DEBUGGING\n");
       CkPrintf("EWALD_recip = OFF FOR DEBUGGING\n");
     }//endif

#endif

//============================================================================
// III) Start grad corr computations

#ifndef CMK_OPTIMIZE
     StartTime=CmiWallTimer();
#endif    

   if(cp_grad_corr_on!=0){
      divRhoVksGspace();
   }else{
    fftcache->freeCacheMem("CP_Rho_GSpacePlane::acceptRhoData");
   }//endif

#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(divRhoVksGspace_, StartTime, CmiWallTimer());    
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
// Setup up options, get box and -i*rho(g)

  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
  complex *packedRho = fftcache->tmpData;
  complex *divRhoX   = rho_gs.divRhoX;
  complex *divRhoY   = rho_gs.divRhoY;
  complex *divRhoZ   = rho_gs.divRhoZ;
  int ioptx          = 1;
  int iopty          = 2; 
  int ioptz          = 3; 
  double tpi,*hmati;

  CPXCFNCTS::CP_fetch_hmati(&hmati,&tpi);
  rho_gs.divRhoGdot(hmati,tpi,packedRho); //divRhoAlpha's are class variable
  fftcache->freeCacheMem("CP_Rho_GSpacePlane::divRhoVksGspace");

//--------------------------------------------------------------------------
// I)   fft_gz(divRhoX), launch transpose to R-space

   fftcache->doRhoFFTGtoR_Gchare(divRhoX,divRhoX,rho_gs.numFull,rho_gs.numPoints,
                                 rho_gs.numLines,rho_gs.numRuns,rho_gs.runs,
                                 rho_gs.sizeZ,0,iplane_ind);
   RhoGSendRhoR(ioptx);

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//--------------------------------------------------------------------------
// II)  fft_gz(divRhoY), launch transpose to R-space

   fftcache->doRhoFFTGtoR_Gchare(divRhoY,divRhoY,rho_gs.numFull,rho_gs.numPoints,
                                 rho_gs.numLines,rho_gs.numRuns,rho_gs.runs,
                                 rho_gs.sizeZ,0,iplane_ind);
   RhoGSendRhoR(iopty);


#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//--------------------------------------------------------------------------
// III)  fft_gz(divRhoZ),  launch transpose to R-space

   fftcache->doRhoFFTGtoR_Gchare(divRhoZ,divRhoZ,rho_gs.numFull,rho_gs.numPoints,
                                rho_gs.numLines,rho_gs.numRuns,rho_gs.runs,
                                rho_gs.sizeZ,0,iplane_ind);
   RhoGSendRhoR(ioptz);

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

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
// Local Pointers and Variables

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int sizeZ    = rho_gs.sizeZ;
  int ix       = thisIndex.x;
  int numLines = rho_gs.numLines;
  int ***index_pack;
  int **nline_send;
  if(rhoRsubplanes>1){
    nline_send  = sim->nline_send_rho_y;
    index_pack  = sim->index_tran_pack_rho_y;
  }//endif

  complex *ffttempdata;
  switch(iopt){
    case 1 : ffttempdata = rho_gs.divRhoX; break;
    case 2 : ffttempdata = rho_gs.divRhoY; break;
    case 3 : ffttempdata = rho_gs.divRhoZ; break;
    case 4 : ffttempdata = rho_gs.divRhoX; break;  //its whitebyrd
    default: CkAbort("impossible iopt");   break;
  }//end switch

//============================================================================
// Do a Comlib Dance

  if(rhoRsubplanes==1){
   switch(iopt){
    case 1 : if(config.useGIns1RhoRP)     commGInstance1.beginIteration();   break;
    case 2 : if(config.useGIns2RhoRP)     commGInstance2.beginIteration();   break;
    case 3 : if(config.useGIns3RhoRP)     commGInstance3.beginIteration();   break;
    case 4 : if(config.useGByrdInsRhoRBP) commGByrdInstance.beginIteration();break;
    default: CkAbort("impossible iopt"); break;
   }//end switc
  }//endif

//============================================================================
// Send the message

  for(int z=0;z<sizeZ;z++) {
  for(int s=0;s<rhoRsubplanes;s++){

    if(rhoRsubplanes>1){numLines = nline_send[ix][s];}
    RhoRSFFTMsg *msg = new (numLines,8*sizeof(int)) RhoRSFFTMsg;
    msg->size        = numLines;     // number of z-lines in this batch
    msg->senderIndex = thisIndex.x;  // line batch index
    msg->iopt        = iopt;         // iopt
    
    if(config.prioFFTMsg){
       CkSetQueueing(msg, CK_QUEUEING_IFIFO);
       *(int*)CkPriorityPtr(msg) = config.rhorpriority + thisIndex.x+ thisIndex.y;
    }//endif

    // beam out all points with same z to chare array index z
    complex *data = msg->data;
    if(rhoRsubplanes==1){
      for (int i=0,j=z; i<numLines; i++,j+=sizeZ){data[i] = ffttempdata[j];}
    }else{
      for(int i=0; i< numLines; i++){
        data[i] = ffttempdata[(z+index_pack[ix][s][i])];
      }//endif
    }//endif

    if(rhoRsubplanes==1){    
     switch(iopt){
      case 1 : rhoRealProxy1_com(z,0).acceptGradRhoVks(msg); break;
      case 2 : rhoRealProxy2_com(z,0).acceptGradRhoVks(msg); break;
      case 3 : rhoRealProxy3_com(z,0).acceptGradRhoVks(msg); break;
      case 4 : rhoRealProxyByrd_com(z,0).acceptWhiteByrd(msg); break;
      default: CkAbort("impossible iopt"); break;
     }//end switch
    }else{
     switch(iopt){
      case 1 : rhoRealProxy(z,s).acceptGradRhoVks(msg); break;
      case 2 : rhoRealProxy(z,s).acceptGradRhoVks(msg); break;
      case 3 : rhoRealProxy(z,s).acceptGradRhoVks(msg); break;
      case 4 : rhoRealProxy(z,s).acceptWhiteByrd(msg); break;
      default: CkAbort("impossible iopt"); break;
     }//end switch
    }//endif

  }//endfor
#ifdef CMK_VERSION_BLUEGENE
       CmiNetworkProgress();
#endif

  }//endfor

//============================================================================
// Complete the commlib dance

  if(rhoRsubplanes==1){    
   switch(iopt){
    case 1 : if(config.useGIns1RhoRP)     commGInstance1.endIteration();break;
    case 2 : if(config.useGIns2RhoRP)     commGInstance2.endIteration();break;
    case 3 : if(config.useGIns3RhoRP)     commGInstance3.endIteration();break;
    case 4 : if(config.useGByrdInsRhoRBP) commGByrdInstance.endIteration();break;
    default: CkAbort("impossible iopt"); break;
   }//end switc
  }//endif

//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//  RhoReal sends rho(gx,gy,z) here such that it is now decomposed 
//  with lines of constant gx,gy in funny order to load balance.
//
//  We cannot receive whitebyrds until every divRho sent above has been received
//  and processed by RhoReal. Therefore, our divRho memory is safe and warm 
//  while it it processing in the routines before the send. It is also safe
//  during the send
//
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

  int iopt             = msg->iopt;
  int size             = msg->size;
  int offset           = msg->offset;   // z index
  int isub             = msg->offsetGx; // subplane index
  complex *partlyIFFTd = msg->data;

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int ix               = thisIndex.x;   // chare array index
  int sizeZ            = rho_gs.sizeZ;
  int numLines         = rho_gs.numLines;
  int ***index_pack;
  if(rhoRsubplanes>1){
    int **nline_send = sim->nline_send_rho_y;
    index_pack       = sim->index_tran_pack_rho_y;
    numLines         = nline_send[ix][isub];
  }//endif

//============================================================================
// Error checking

#ifdef _NAN_CHECK_
  for(int i=0;i<msg->size ;i++){
      CkAssert(isnan(msg->data[i].re)==0);
      CkAssert(isnan(msg->data[i].im)==0);
    }
#endif
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
    default: CkAbort("impossible iopt"); break;
  }//end switch

  if(rhoRsubplanes==1){
    for(int i=0,j=offset; i< size; i++,j+=sizeZ){chunk[j] = partlyIFFTd[i];}
  }else{
    for(int i=0; i< size; i++){
      chunk[(offset+index_pack[ix][isub][i])] = partlyIFFTd[i];
    }//endif
  }//endif

  delete msg;

//============================================================================
// If all chares for a gradient report, perform final FFT.

  countWhiteByrd[iopt]++;
  if(countWhiteByrd[iopt]==sizeZ*rhoRsubplanes){
    countWhiteByrd[iopt]=0;
    doneWhiteByrd++;

#ifndef CMK_OPTIMIZE
    double StartTime=CmiWallTimer();
#endif    
    FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
    fftcache->doRhoFFTRtoG_Gchare(chunk,chunk,rho_gs.numFull,rho_gs.numPoints,
                                  rho_gs.numLines,rho_gs.numRuns,rho_gs.runs,
                                  rho_gs.sizeZ,0,iplane_ind);
#ifndef CMK_OPTIMIZE
    traceUserBracketEvent(BwFFTRtoG_, StartTime, CmiWallTimer());    
#endif

  }//endif

//============================================================================
// When all 3 gradients are in g-space, then you will ready for the next step.

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
// Set the timer and some options

#ifndef CMK_OPTIMIZE
  double StartTime=CmiWallTimer();
#endif    

  int ioptFFTWhite  = 0; 
  int ioptSendWhite = 4; 

//============================================================================
// Compute my whiteByrd : store it in divrhox 

  double tpi,*hmati;
  CPXCFNCTS::CP_fetch_hmati(&hmati,&tpi);
  rho_gs.createWhiteByrd(hmati,tpi);

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();    
#endif

//============================================================================
// FFT my whitebyrd : which is stored in divRhox 

  FFTcache *fftcache = fftCacheProxy.ckLocalBranch();  
  complex *white = rho_gs.divRhoX;
  fftcache->doRhoFFTGtoR_Gchare(white,white,rho_gs.numFull,rho_gs.numPoints,
                                rho_gs.numLines,rho_gs.numRuns,rho_gs.runs,
                                rho_gs.sizeZ,0,iplane_ind);
#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();    
#endif

//============================================================================
// End tracing and send whitebyrd back to Rhoreal space

#ifndef CMK_OPTIMIZE
  traceUserBracketEvent(ByrdanddoFwFFTGtoR_, StartTime, CmiWallTimer());    
#endif

  RhoGSendRhoR(ioptSendWhite);

//---------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_Rho_GSpacePlane::ResumeFromSync(){

  if (config.useGIns0RhoRP) 
    ComlibResetProxy(&rhoRealProxy0_com);
  if (config.useGIns1RhoRP) 
    ComlibResetProxy(&rhoRealProxy1_com);
  if (config.useGIns2RhoRP)
    ComlibResetProxy(&rhoRealProxy2_com);
  if (config.useGIns3RhoRP) 
    ComlibResetProxy(&rhoRealProxy3_com);
  if (config.useGByrdInsRhoRBP) 
    ComlibResetProxy(&rhoRealProxyByrd_com);

}
//============================================================================



//============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
// Glenn's RhoReal exit 
//============================================================================
void CP_Rho_GSpacePlane::exitForDebugging(){
  countDebug++;  
  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;
  int nchareG          = sim->nchareRhoG;
  if(countDebug==nchareG){
    countDebug=0;
    CkPrintf("I am in the exitfordebuging rhoG puppy. Bye-bye\n");
    CkExit();
  }//endif
}
//============================================================================
