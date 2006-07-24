//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file fftCache.C
 * Add functions to allow application programmers to initialize these and
 * the corresponding functions in Charm++ to call these functions with
 * appropriate parameters 
 */
//==============================================================================

#include "util.h"
#include "../../include/debug_flags.h"
#include <math.h>
#include "cpaimd.h"
#include "groups.h"
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"

//==============================================================================

extern CProxy_FFTcache fftCacheProxy;
extern Config config;
extern int nstates;
extern int sizeX;
extern CProxy_AtomsGrp atomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;

//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
//FFTcache - sets up a fftcache on each processor
//==============================================================================
FFTcache::FFTcache(size2d planeSIZE, int _ngridaEext, int _ngridbEext, int _ngridcEext, 
                   int _ees_eext_on, int _ngridaNL, int _ngridbNL, int _ngridcNL, 
                   int _ees_NL_on){
//==============================================================================
// Local Variables

    int size[3];
    int sizeZ;

//==============================================================================
// Copy out

    planeSize   = planeSIZE;    //contains nfftz,nfftx for non-ees stuff

    ngridaEext  = _ngridaEext; 
    ngridbEext  = _ngridbEext;  
    ngridcEext  = _ngridcEext; 
    ees_eext_on = _ees_eext_on; 

    ngridaNL    = _ngridaNL; 
    ngridbNL    = _ngridbNL; 
    ngridcNL    = _ngridcNL; 
    ees_NL_on   = _ees_NL_on;

//==============================================================================
// Density and State Scratch

    int nlines_max = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_max;
    int psize      = nlines_max*planeSize[0];

//==============================================================================
// Density and State plans/scratch : funky names : non-cubic broken

   //------------------------------------------------
   // Double Pack Plans
    size[0] = planeSize[1]; size[1] = planeSize[0]; size[2] = 1;
    sizeZ   = planeSize[0];

    // really y plans
    fwdZ1DdpPlan = fftw_create_plan(planeSize[0], FFTW_FORWARD, 
                    FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdZ1DdpPlan = fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                    FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    // x is correct!
    fwdX1DdpPlan = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL, 
                    FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdX1DdpPlan = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                    FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    // really z plans
    fwdYPlan = fftw_create_plan(sizeZ,FFTW_FORWARD, 
                    FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdYPlan = fftw_create_plan(sizeZ,FFTW_BACKWARD, 
                    FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM); 

//==============================================================================
// Euler spline scratch and sizes

    int nlines_max_rho = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_max_rho;
    int psizeEext      = nlines_max_rho*ngridcEext;
    int psizeNL        = nlines_max*ngridcNL;

    int pmax = psize;
    if(ees_eext_on==1){pmax = (pmax>psizeEext ? pmax : psizeEext);}
    if(ees_NL_on==1)  {pmax = (pmax>psizeNL   ? pmax : psizeNL);}

    tmpData = (complex*) fftw_malloc(pmax*sizeof(complex)); 

//==============================================================================
// Euler spline plans : better names : non-cubic broken

   //-----------------------------------
   // Double Pack Plan NL : 
    size[0] = ngridaNL; size[1] = ngridbNL; size[2] = 1;

    fwdY1DdpPlanNL = fftw_create_plan(ngridbNL, FFTW_FORWARD, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdY1DdpPlanNL = fftw_create_plan(ngridbNL, FFTW_BACKWARD, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    fwdX1DdpPlanNL = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdX1DdpPlanNL = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    fwdZPlanNL = fftw_create_plan(ngridcNL,FFTW_FORWARD, 
                       FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdZPlanNL = fftw_create_plan(ngridcNL,FFTW_BACKWARD, 
                       FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);

   //-----------------------------------
   // Double Pack Plan Eext : 
    size[0] = ngridaEext; size[1] = ngridbEext; size[2] = 1;

    fwdY1DdpPlanEext = fftw_create_plan(ngridbEext, FFTW_FORWARD, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdY1DdpPlanEext = fftw_create_plan(ngridbEext, FFTW_BACKWARD, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    fwdX1DdpPlanEext = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdX1DdpPlanEext = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
    fwdZPlanEext = fftw_create_plan(ngridcEext,FFTW_FORWARD, 
                       FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
    bwdZPlanEext = fftw_create_plan(ngridcEext,FFTW_BACKWARD, 
                       FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);

//------------------------------------------------------------------------------
   }//end routine
//==============================================================================



//=============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void FFTcache::expandGSpace(complex* data, complex *packedData, 
                            RunDescriptor *runs, int numRuns, int numFull,
                            int numPoints, int nfftz)

//==============================================================================
  {//begin routine
//==============================================================================
// Expand out lines of z so that an FFT can be performed on them
//   The ``run'' stores each line in two parts
//     -3,-2,-1 have offset, joff = r*nffz + nfftz
//      0,1,2,3 have offset, joff = r*nfftz 
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0
//      The negative guys go 1st
//   Total size is nlines*nfftz where nlines=numRuns/2

  bzero(data,numFull*sizeof(complex));
  int nsub = (nfftz-runs[0].nz);
  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z + nsub; // k < 0
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      data[j] = packedData[k];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z; // k >= 0
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      data[j] = packedData[k];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void FFTcache::packGSpace(complex* data, complex *packedData, 
                            RunDescriptor *runs, int numRuns, int numFull,
                            int numPoints, int nfftz)

//==============================================================================
  {//begin routine
//==============================================================================
// Contract lines of z after FFT has been performed
//   The ``run'' stores each line in two parts
//     -3,-2,-1 have offset, joff = r*nffz + nfftz
//      0,1,2,3 have offset, joff = r*nfftz 
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0
//      The negative guys go 1st
//   Total size is nlines*nfftz where nlines=numRuns/2

  int koff = 0;
  int nsub = (nfftz-runs[0].nz);
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z + nsub;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      packedData[k] = data[j];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      packedData[k] = data[j];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// non-local : Gchare : data(gx,gy,z) -> data(gx,gy,gz) : backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doNlFftRtoG_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,
                                  int numLines, int numRuns, RunDescriptor *runs, 
                                  int nfftz){
//==============================================================================
// FFT in expanded form

  fft_split(
          bwdZPlanNL,             // Z-direction backward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data_in,   //input data
	  1,                      //stride
	  nfftz,                  //distance between z-data sets
	  NULL, 0, 0,             // input is ouput
          config.fftprogresssplit // split parameter
         ); 

//==============================================================================
// Pack for computing

  packGSpace(data_in,data_out,runs,numRuns,numFull,numPoints,nfftz); // data_out is upacked

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// non-local : Gchare : data(gx,gy,gz) -> data(gx,gy,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doNlFftGtoR_Gchare(complex *data_in,complex *data_out,
          int numFull, int numPoints,int numLines, int numRuns, 
          RunDescriptor *runs, int nfftz){
//==============================================================================
// Expand for ffting : Trickery

  expandGSpace(data_out,data_in,runs,numRuns,numFull,numPoints,nfftz);//data_out is expanded

//==============================================================================
// FFT in expanded form

  fft_split(
          fwdZPlanNL,               // Z-direction forward plan
          numLines,                 // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data_out, // data
	  1,                        //stride
	  nfftz,                    //distance between z-data sets
	  NULL, 0, 0,               // input is ouput
          config.fftprogresssplit   // split parameter
         ); 

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================
 

//=============================================================================
// non-local : Rchare : data(x,y,z) -> data(gx,gy,z) : backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doNlFFTRtoG_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                  int sizeX,int sizeY){
//==============================================================================
// FFT along x first

  rfftwnd_real_to_complex_split(
   	     bwdX1DdpPlanNL,             // backward plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}

//==============================================================================
// FFT along y

  fft_split(
            bwdY1DdpPlanNL,            // backward plan
	    nplane_x,                  // these many 1D ffts
	    (fftw_complex *)dataC,     // data set
            stride,                    // stride
            1,                         // spacing between data sets
            NULL,0,0,                  // input is output
            config.fftprogresssplitReal// progress splitting for BG/L
           );

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// non-local : Rchare : data(gx,gy,z) -> data(x,y,z) : Forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doNlFFTGtoR_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                  int sizeX,int sizeY){
//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory
//                       : nplane_x is spherical cutoff <= sizeX/2+1

  int stride = sizeX/2+1;
  fft_split(
        fwdY1DdpPlanNL,               // y-plan for NL Ees method
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	stride,                       // stride betwen elements (x is inner)
  	1,                            // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
  );

  // fftw only gives you one sign for real to complex : so do it yourself
  for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(
              fwdX1DdpPlanNL,               // x-plan for NL Ees method
	      sizeY,                        // how many
	      (fftw_complex *)dataC,        // input data
              1,                            // stride (x is inner)
              stride,                       // array separation
              NULL,0,0,                     // output = input = dataR real
  	      config.fftprogresssplitReal   // splitting parameter
   );

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================



//=============================================================================
// Eext : Gchare : data(gx,gy,z) -> data(gx,gy,gz) : backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTRtoG_Gchare(complex *data,int numFull, int numPoints,
                                    int numLines, int numRuns, RunDescriptor *runs, 
                                    int nfftz){
//==============================================================================
// FFT in expanded form

  fft_split(
          bwdZPlanEext,           // Z-direction backward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data,   //input data
	  1,                      //stride
	  nfftz,                  //distance between z-data sets
	  NULL, 0, 0,             // input is ouput
          config.fftprogresssplit // split parameter
         ); 

//==============================================================================
// pack for computing

  packGSpace(data,tmpData,runs,numRuns,numFull,numPoints,nfftz); // data is packed
  memcpy(data,tmpData,sizeof(complex)*numPoints); // data is expanded; cp to tmp; 

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// Eext : Gchare : data(gx,gy,gz) -> data(gx,gy,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTGtoR_Gchare(complex *data, int numFull, int numPoints,
                                    int numLines, int numRuns, RunDescriptor *runs,
                                    int nfftz){
//==============================================================================
// expand for ffting

  memcpy(tmpData,data,sizeof(complex)*numPoints); // data is packed; cp to tmp; 
  expandGSpace(data,tmpData,runs,numRuns,numFull,numPoints,nfftz);// data is expanded

//==============================================================================
// FFT

  fft_split(
          fwdZPlanEext,           // Z-direction forward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data,   //input data
	  1,                      //stride
	  nfftz,                  //distance between z-data sets
	  NULL, 0, 0,             // input is ouput
          config.fftprogresssplit // split parameter
         ); 

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// Eext : Rchare : data(x,y,z) -> data(gx,gy,z) : Backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTRtoG_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                    int sizeX,int sizeY){
//==============================================================================
// FFT along x first

  rfftwnd_real_to_complex_split(
   	     bwdX1DdpPlanEext,             // backward plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}

//==============================================================================
// FFT along y

  fft_split(
            bwdY1DdpPlanEext,          // backward plan
	    nplane_x,                  // these many 1D ffts
	    (fftw_complex *)dataC,     // data set
            stride,                    // stride
            1,                         // spacing between data sets
            NULL,0,0,                  // input is output
            config.fftprogresssplitReal// progress splitting for BG/L
           );

//------------------------------------------------------------------------------
 }//end routine
//==============================================================================


//=============================================================================
// Eext : Rchare : data(gx,gy,z) -> data(x,y,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTGtoR_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                    int sizeX,int sizeY){
//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory
//                       : nplane_x is spherical cutoff <= sizeX/2+1

  int stride = sizeX/2+1;
  fft_split(
        fwdY1DdpPlanEext,             // y-plan for NL Ees method
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	stride,                       // stride betwen elements (x is inner)
  	1,                            // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
       );

  // fftw only gives you one sign for real to complex : so do it yourself
  for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(
              fwdX1DdpPlanEext,             // x-plan for NL Ees method
	      sizeY,                        // how many
	      (fftw_complex *)dataC,        // input data
              1,                            // stride (x is inner)
              stride,                       // array separation
              NULL,0,0,                     // output = input = dataR real
  	      config.fftprogresssplitReal   // splitting parameter
             );

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// StatePlane : Gchare : data(gx,gy,z) -> data(gx,gy,gz) : backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doStpFftRtoG_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,
                                  int numLines, int numRuns, RunDescriptor *runs, 
                                  int nfftz){
//==============================================================================
// FFT in expanded form

  fft_split(
          bwdYPlan,               // Z-direction backward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data_in,//input data
	  1,                      //stride
	  nfftz,                  //distance between z-data sets
	  NULL, 0, 0,             // input is ouput
          config.fftprogresssplit // split parameter
         ); 

//==============================================================================
// Pack for computing

  packGSpace(data_in,data_out,runs,numRuns,numFull,numPoints,nfftz);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// StatePlane : Gchare : data(gx,gy,gz) -> data(gx,gy,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doStpFftGtoR_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,
                                  int numLines, int numRuns, RunDescriptor *runs, 
                                  int nfftz){
//==============================================================================
// Expand for ffting

  // data_out is expanded
  // data_in is contracted
  expandGSpace(data_out,data_in,runs,numRuns,numFull,numPoints,nfftz);

//==============================================================================
// FFT in expanded form

  fft_split(
          fwdYPlan,                // Z-direction forward plan
          numLines,                // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data_out,//input data
	  1,                       //stride
	  nfftz,                   //distance between z-data sets
	  NULL, 0, 0,              // input is ouput
          config.fftprogresssplit  // split parameter
         ); 

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================
 



//==============================================================================
// After the transpose : nfftz*nffty array nicely expanded out.
//                       Only nplane_x FFTs along Y need to be performed
//                       due to spherical truncation
//                       data is |psi|^2 : z is chare array index
//                                         x is inner index 
//                                         y is outer index
//                       chare[z].data[x+y*sizeX]
//                       chare[z].planeArr[y+x*sizeY] is input
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

double* FFTcache::doRealFwFFT(complex *planeArr)

//==============================================================================
   {//begin routine    
//==============================================================================

   double *data;
   int pSize    = planeSize[0] * sizeX;
   int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;

//==============================================================================
// Case doublePack && config.inPlaceFFT
// 

  if(config.doublePack) {
      data = (double *)fftw_malloc(pSize*sizeof(double));
      int stride = sizeX/2+1;
      // FFT along Y direction : Y moves with stride sizex/2+1 through memory
      fft_split(fwdZ1DdpPlan, // y-plan (label lies)
  	   nplane_x,          // how many < sizeX/2 + 1
	   (fftw_complex *)(planeArr),//input data
 	   stride,            // stride betwen elements (y is inner)
  	   1,                 // array separation (nffty elements)
  	   NULL,0,0,          // output data is input data 
           config.fftprogresssplitReal);
      // fftw only gives you one sign for real to complex : so do it yourself
      for(int i=0;i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}
      rfftwnd_complex_to_real_split(fwdX1DdpPlan,
			      planeSize[0],    // how many
			      (fftw_complex *)planeArr, 
                              1,               // stride (x is inner)
                              stride,          // array separation
			      NULL,0,0,        // output array is real 
                              config.fftprogresssplitReal);      
      // x is now the inner index as fftw has transposed for us
      double *realArr = reinterpret_cast<double*> (planeArr);
      for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
       for(int j=i*sizeX;j<(i+1)*sizeX;j++){
         data[j] = realArr[(j+i2)]*realArr[(j+i2)];
       }//endfor
      }//endfor

    return data;

  }//endif : 

//==============================================================================
// Case !doublePack && config.inPlaceFFT

  if(!config.doublePack){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Fix the non-double pack FFT\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

  return data;                

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoRealtoRhoG(double *realArr)
//==============================================================================
   {//begin routine 
//==============================================================================

   int nplane_rho_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_rho_x;

//==============================================================================
// Case : doublePack  and inplace no options thats how it is
   
   rfftwnd_real_to_complex_split(bwdX1DdpPlan,
			   planeSize[0],          // these many 1D ffts
			   realArr,1,(sizeX+2),   // x is inner here
			   NULL,0,0,              // output is input
                           config.fftprogresssplitReal);            
   complex *planeArr = reinterpret_cast<complex*> (realArr);
   int stride = sizeX/2+1;

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

   for (int i=0; i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}

   fft_split(bwdZ1DdpPlan,
	nplane_rho_x, // these many 1D ffts
	(fftw_complex *)planeArr, 
	 stride,1,
         NULL,0,0,
         config.fftprogresssplitReal);

//------------------------------------------------------------------------------
 }//end routine
//==============================================================================


//==============================================================================
//                             vks : z is chare array index
//                                   x is inner index 
//                                   y is outer index
//                                   chare[z].data[x+y*sizeX]
//
//          vks is from a message and should not be modified
//
//                           planeArr : chare[z].data[y+x*sizeY]
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRealBwFFT(const double *vks, complex *planeArr,
                           int ind_state, int ind_plane)
//==============================================================================
 {//begin routine 
//==============================================================================

   int pSize    = planeSize[0] * sizeX;
   int nplane_x = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;

//==============================================================================
// Output

#ifdef _CP_DEBUG_FFTR_VKSR_
   if(config.doublePack){
     if(ind_state==0 && ind_plane==0){
        double *realArr = reinterpret_cast<double*> (planeArr);
        FILE *fp = fopen("vks_state0_plane0_real.out","w");
          for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
            for(int j=i*sizeX;j<(i+1)*sizeX;j++){
              fprintf(fp,"%g %g\n",vks[j],realArr[(j+i2)]);
            }//endfor
          }//endfor
        fclose(fp);
     }//endif
   }//endif
#endif

//==============================================================================
// Case : doublePack  and inplace 

  if(config.doublePack){
     double *realArr = reinterpret_cast<double*> (planeArr);
     int stride = sizeX/2+1;
     for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
       for(int j=i*sizeX;j<(i+1)*sizeX;j++){
         realArr[(j+i2)] = realArr[(j+i2)]*vks[j];
       }//endfor
     }//endfor
     rfftwnd_real_to_complex_split(bwdX1DdpPlan,
			planeSize[0],          // these many 1D ffts
			realArr,
                        1,(sizeX+2),   // x is inner here
  		        NULL,0,0, 
                        config.fftprogresssplitReal);            
#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

     for (int i=0; i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}
     fft_split(bwdZ1DdpPlan,
	     nplane_x, // these many 1D ffts
	     (fftw_complex *)planeArr, 
	     stride,1,
             NULL,0,0,
             config.fftprogresssplitReal);
  }//endif : double pack and in-place

//==============================================================================

  if(!config.doublePack){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Fix the non-double pack FFT\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//------------------------------------------------------------------------------
   }//end routine
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void fft_split(fftw_plan plan, int howmany, fftw_complex *in, int istride,
	       int idist, fftw_complex *out, int ostride, int odist, int split)
//============================================================================
  {// begin routine 
//============================================================================

/*  int sizefft=howmany*idist*istride;
  CkPrintf("plan n %d, howmany %d istride %d idist %d ostride %d odist %d sizefft %d\n",
            plan->n,howmany, istride, idist, ostride, odist, sizefft);
  // for debugging 
  fftw_complex *scratch1=(fftw_complex *) fftw_malloc(sizefft*sizeof(complex));
  memcpy(scratch1,in,sizeof(complex)*sizefft);
  fftw(plan,         // da plan
       howmany,     // how many 
       scratch1,           //input data
       istride,      // stride betwen elements 
       idist,        // array separation
       out,          // output data (null if inplace)
       ostride,      // stride betwen elements  (0 inplace)
       odist);       //array separation (0 inplace)


  fftw_complex *scratch2=(fftw_complex *) fftw_malloc(sizefft*sizeof(complex));
  memcpy(scratch2,in,sizeof(complex)*sizefft);
*/

//============================================================================
//    fftw_complex *myin=scratch2;

  fftw_complex *myin=in;
  int thismany=split;
  int inleft=howmany;
  int numsplits=howmany/split;
  if(numsplits<=0)
    {numsplits=1;}
  if(howmany>split && howmany%split!=0)
    {numsplits++;}
  
  int inoff=0;
  int outoff=0;
  //  CkPrintf("numsplits %d\n",numsplits);
  for(int i=0;i<numsplits;i++){
      thismany=split;
      if(inleft<split)
	{thismany=inleft;}
      //      CkPrintf("split %d thismany %d inleft %d\n",i, thismany,inleft);
      fftw(plan,           // da plan
	   thismany,       // how many 
	   &(myin[inoff]), //input data
	   istride,        // stride betwen elements 
	   idist,          // array separation
	   (out==NULL) ? NULL : &(out[outoff]), // output data (null if inplace)
	   ostride,        // stride betwen elements  (0 inplace)
	   odist);         //array separation (0 inplace)

      inoff+=idist*(thismany);  
      outoff+=odist*(thismany);
      inleft-=thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
  }//endfor

//============================================================================

  /*  for(int j=0;j<sizefft;j++)
    {
      CkPrintf("%g %g\n",scratch1[j].re,scratch1[j].im);
    }
  CkPrintf("----\n");
  for(int j=0;j<sizefft;j++)
    {
      CkPrintf("%g %g\n",scratch2[j].re,scratch2[j].im);
    }
  CkAbort("fftd\n");
  */

//------------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void rfftwnd_complex_to_real_split(rfftwnd_plan plan, int howmany, fftw_complex *in, 
         int istride,int idist, fftw_real *out, int ostride, int odist, int split)
//============================================================================
   {//begin routine
//============================================================================

  int thismany=split;
  int inleft=howmany;
  int numsplits=howmany/split;
  if(numsplits<=0)
    {numsplits=1;}
  if(howmany>split && howmany%split!=0)
    {numsplits++;}
  
  int inoff=0;
  int outoff=0;
  //  CkPrintf("numsplits %d\n",numsplits);
  for(int i=0;i<numsplits;i++){
      thismany=split;
      if(inleft<split)
	{thismany=inleft;}
      rfftwnd_complex_to_real(plan,// da plan
	   thismany,               // how many 
	   &(in[inoff]),           //input data
	   istride,                // stride betwen elements 
	   idist,                  // array separation
	   (out==NULL) ? NULL : &(out[outoff]),          // output data (null if inplace)
	   ostride,                // stride betwen elements  (0 inplace)
	   odist);                 //array separation (0 inplace)

      inoff+=idist*(thismany);  
      outoff+=odist*(thismany);
      inleft-=thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
    }//endfor

//------------------------------------------------------------------------------
   }// end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void  rfftwnd_real_to_complex_split(rfftwnd_plan plan, int howmany, fftw_real *in, 
      int istride, int idist, fftw_complex *out, int ostride, int odist, int split)
//============================================================================
   {// begin routine
//============================================================================

  int thismany=split;
  int inleft=howmany;
  int numsplits=howmany/split;
  if(numsplits<=0)
    {numsplits=1;}
  if(howmany>split && howmany%split!=0)
    {numsplits++;}
  
  int inoff=0;
  int outoff=0;
  //  CkPrintf("numsplits %d\n",numsplits);
  for(int i=0;i<numsplits;i++){
      thismany=split;
      if(inleft<split)
	{thismany=inleft;}
      //      CkPrintf("split %d thismany %d inleft %d\n",i, thismany,inleft);
      rfftwnd_real_to_complex(plan,         // da plan
	   thismany,     // how many 
	   &(in[inoff]), //input data
	   istride,      // stride betwen elements 
	   idist,        // array separation
	   (out==NULL) ? NULL : &(out[outoff]), // output data (null if inplace)
	   ostride,      // stride betwen elements  (0 inplace)
	   odist);       //array separation (0 inplace)

      inoff+=idist*(thismany);  
      outoff+=odist*(thismany);
      inleft-=thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
    }// endfor

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


