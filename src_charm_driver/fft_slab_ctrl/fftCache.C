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
#include "../../src_mathlib/mathlib.h"

//#define TEST_ALIGN
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
                   int _ees_NL_on, int _nlines_max, int _nlines_max_rho){
//==============================================================================
// Local Variables

    int size[3];
    complex *cin; complex *cout; 
    double  *din; double *dout; 

//==============================================================================
// Copy out

    planeSize    = planeSIZE;    //contains nfftz,nfftx for non-ees stuff

    ngridaEext   = _ngridaEext; 
    ngridbEext   = _ngridbEext;  
    ngridcEext   = _ngridcEext; 
    ees_eext_on  = _ees_eext_on; 

    ngridaNL     = _ngridaNL; 
    ngridbNL     = _ngridbNL; 
    ngridcNL     = _ngridcNL; 
    ees_NL_on    = _ees_NL_on;
    int nlines_max     = _nlines_max;
    int nlines_max_rho = _nlines_max_rho;

    int sizeY    = planeSIZE[0];
    int sizeZ    = planeSIZE[1];

    int iopt     = config.fftopt;

    int skipR    =  1;
    int skipC    =  1;
    int nf_max   =  0;
    int nwork1   =  0; 
    int nwork2   =  0;
    int plus     =  1;
    int mnus     = -1;
    int unit     =  1;
    double scale =  1.0;
  

    if(iopt<0 || iopt>1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Bad fft option, dude! %d\n",iopt);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif

//==============================================================================
// Density, State and EES Scratch


    int pGsize         = nlines_max*planeSize[1];     // z-collections
    int pGsizeHart     = nlines_max_rho*planeSize[1]; // z-collections
    int pGsizeNL       = nlines_max*ngridcNL;         // z-collections
    int pGsizeEext     = nlines_max_rho*ngridcEext;   // z-collections
    int pRsize         = (sizeX/2+1)*planeSize[0];    // x-y plane
    int pRsizeNL       = (ngridaNL/2+1)*ngridbNL;     // x-y plane

    int pmax = pGsize; 
    pmax     = (pmax > pRsize   )  ? pmax : pRsize;
    pmax     = (pmax > pGsizeHart) ? pmax : pGsizeHart;

    if(ees_eext_on==1){pmax = (pmax>pGsizeEext) ? pmax : pGsizeEext;}
    if(ees_NL_on  ==1){pmax = (pmax>pGsizeNL)   ? pmax : pGsizeNL;}
    if(ees_NL_on  ==1){pmax = (pmax>pRsizeNL)   ? pmax : pRsizeNL;}

    tmpData  = (complex*) fftw_malloc(pmax*sizeof(complex)); 
    tmpDataR = reinterpret_cast<double*> (tmpData);

//==============================================================================
// Density and State plans : funky names : non-cubic broken

   //------------------------------------------------
   // Double Pack Plans

    nf_max = (sizeX > sizeY)  ? sizeX : sizeY; 
    nf_max = (sizeZ > nf_max) ? sizeZ : nf_max;
    nwork1 = 0;    nwork2 = 0;
    if(iopt==1){    
      int iadd = (int) (2.3*(double)nf_max);
      nwork1   = 20000+iadd;
      nwork2   = (2*nf_max+256)*64;
      nwork2   = (nwork1 > nwork2) ? nwork1 : nwork2;
    }//endif
    skipR = sizeX+2;
    skipC = sizeX/2+1;
    
    initFFTholder  (&fwdYPlanRho, &iopt,&nwork1,&nwork2,&scale,&plus,&sizeY,&skipC,&unit);
    initFFTholder  (&bwdYPlanRho, &iopt,&nwork1,&nwork2,&scale,&mnus,&sizeY,&skipC,&unit);
    initFFTholder  (&fwdYPlanRhoS,&iopt,&nwork1,&nwork2,&scale,&plus,&sizeY,&unit, &sizeY);
    initFFTholder  (&bwdYPlanRhoS,&iopt,&nwork1,&nwork2,&scale,&mnus,&sizeY,&unit, &sizeY);
    initCRFFTholder(&fwdXPlanRho, &iopt,&nwork1,&nwork2,&scale,&plus,&sizeX,&skipR,&skipC);
    initRCFFTholder(&bwdXPlanRho, &iopt,&nwork1,&nwork2,&scale,&mnus,&sizeX,&skipR,&skipC);
    initFFTholder  (&fwdZPlanRho, &iopt,&nwork1,&nwork2,&scale,&plus,&sizeZ,&unit, &sizeZ);
    initFFTholder  (&bwdZPlanRho, &iopt,&nwork1,&nwork2,&scale,&mnus,&sizeZ,&unit, &sizeZ);

    if(iopt==0){
        size[0] = planeSize[1]; size[1] = planeSize[0]; size[2] = 1; 
        fwdXPlanRho.rfftwPlan = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL, 
                                   FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
        bwdXPlanRho.rfftwPlan = rfftwnd_create_plan(1, (const int*)size, FFTW_REAL_TO_COMPLEX,
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        fwdYPlanRho.fftwPlan  =  fftw_create_plan(planeSize[0], FFTW_FORWARD, 
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdYPlanRho.fftwPlan =   fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        fwdYPlanRhoS.fftwPlan =  fftw_create_plan(planeSize[0], FFTW_FORWARD, 
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdYPlanRhoS.fftwPlan  =  fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        fwdZPlanRho.fftwPlan  =  fftw_create_plan(sizeZ,FFTW_FORWARD, 
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdZPlanRho.fftwPlan  =  fftw_create_plan(sizeZ,FFTW_BACKWARD, 
                                   FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM); 
    }//endif

//==============================================================================
// Euler spline plans : better names : non-cubic broken

   //-----------------------------------
   // Double Pack Plan NL : 

    nf_max = (ngridaNL> ngridbNL) ? ngridaNL : ngridbNL; 
    nf_max = (ngridcNL > nf_max ) ? ngridcNL : nf_max; 
    nwork1 = 0;   nwork2 = 0;
    if(iopt==1){    
      int iadd = (int) (2.3*(double)nf_max);
      nwork1   = 20000+iadd;
      nwork2   = (2*nf_max+256)*64;
      nwork2   = (nwork1 > nwork2) ? nwork1 : nwork2;
    }//endif
    skipR = ngridaNL+2;
    skipC = ngridaNL/2+1;
    if(ees_NL_on)
      {
	initFFTholder  (&fwdYPlanNL,&iopt,&nwork1,&nwork2,&scale,&plus,&ngridbNL,&skipC,&unit);
	initFFTholder  (&bwdYPlanNL,&iopt,&nwork1,&nwork2,&scale,&mnus,&ngridbNL,&skipC,&unit);
	initCRFFTholder(&fwdXPlanNL,&iopt,&nwork1,&nwork2,&scale,&plus,&ngridaNL,&skipR,&skipC);
	initRCFFTholder(&bwdXPlanNL,&iopt,&nwork1,&nwork2,&scale,&mnus,&ngridaNL,&skipR,&skipC);
	initFFTholder  (&fwdZPlanNL,&iopt,&nwork1,&nwork2,&scale,&plus,&ngridcNL,&unit, &ngridcNL);
	initFFTholder  (&bwdZPlanNL,&iopt,&nwork1,&nwork2,&scale,&mnus,&ngridcNL,&unit, &ngridcNL);
      }
    if(iopt==0){
       size[0] = ngridaNL; size[1] = ngridbNL; size[2] = 1;
       fwdXPlanNL.rfftwPlan = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL,
                                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
       bwdXPlanNL.rfftwPlan = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                                     FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
       fwdYPlanNL.fftwPlan  =  fftw_create_plan(ngridbNL, FFTW_FORWARD, 
                                     FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
       bwdYPlanNL.fftwPlan  =  fftw_create_plan(ngridbNL, FFTW_BACKWARD, 
                                     FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
       fwdZPlanNL.fftwPlan  =  fftw_create_plan(ngridcNL,FFTW_FORWARD, 
                                     FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
       bwdZPlanNL.fftwPlan  =  fftw_create_plan(ngridcNL,FFTW_BACKWARD, 
                                     FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
    }//end if

//==============================================================================
// Euler spline plans : better names : non-cubic broken

   //-----------------------------------
   // Double Pack Plan Eext : 

    nf_max = (ngridaEext> ngridbEext) ? ngridaEext : ngridbEext; 
    nf_max = (ngridcEext > nf_max )   ? ngridcEext : nf_max; 
    nwork1 = 0;   nwork2 = 0;
    if(iopt==1){    
      int iadd = (int) (2.3*(double)nf_max);
      nwork1   = 20000+iadd;
      nwork2   = (2*nf_max+256)*64;
      nwork2   = (nwork1 > nwork2) ? nwork1 : nwork2;
    }//endif
    skipR = ngridaEext+2;
    skipC = ngridaEext/2+1;
    if(ees_eext_on)
      {
	initFFTholder  (&fwdYPlanEext ,&iopt,&nwork1,&nwork2,&scale,&plus,&ngridbEext,
			&skipC,&unit);
	initFFTholder  (&bwdYPlanEext ,&iopt,&nwork1,&nwork2,&scale,&mnus,&ngridbEext,
			&skipC,&unit);
	initFFTholder  (&fwdYPlanEextS,&iopt,&nwork1,&nwork2,&scale,&plus,&ngridbEext,
			&unit, &ngridbEext);
	initFFTholder  (&bwdYPlanEextS,&iopt,&nwork1,&nwork2,&scale,&mnus,&ngridbEext,
			&unit, &ngridbEext);
	initCRFFTholder(&fwdXPlanEext,& iopt,&nwork1,&nwork2,&scale,&plus,&ngridaEext,
			&skipR,&skipC);
	initRCFFTholder(&bwdXPlanEext,& iopt,&nwork1,&nwork2,&scale,&mnus,&ngridaEext,
			&skipR,&skipC);
	initFFTholder  (&fwdZPlanEext,& iopt,&nwork1,&nwork2,&scale,&plus,&ngridcEext,
			&unit, &ngridcEext);
	initFFTholder  (&bwdZPlanEext,& iopt,&nwork1,&nwork2,&scale,&mnus,&ngridcEext,
			&unit, &ngridcEext);
      }
    if(iopt==0){
        size[0] = ngridaEext; size[1] = ngridbEext; size[2] = 1;
        fwdXPlanEext.rfftwPlan = rfftwnd_create_plan(1, (const int*)size,FFTW_COMPLEX_TO_REAL,
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdXPlanEext.rfftwPlan = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        fwdYPlanEext.fftwPlan  =  fftw_create_plan(ngridbEext, FFTW_FORWARD, 
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdYPlanEext.fftwPlan  =  fftw_create_plan(ngridbEext, FFTW_BACKWARD, 
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        fwdYPlanEextS.fftwPlan =  fftw_create_plan(ngridbEext, FFTW_FORWARD, 
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdYPlanEextS.fftwPlan =  fftw_create_plan(ngridbEext, FFTW_BACKWARD, 
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        fwdZPlanEext.fftwPlan  =  fftw_create_plan(ngridcEext,FFTW_FORWARD, 
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
        bwdZPlanEext.fftwPlan  =  fftw_create_plan(ngridcEext,FFTW_BACKWARD, 
                                       FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
    }//end switch

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

#define _BZERO_METH_
#ifdef _BZERO_METH_
    bzero(data,sizeof(complex)*numFull);
#endif

  int nsub = (nfftz-runs[0].nz);
  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff1 = l*nfftz + runs[r].z + nsub; // k < 0
    for (int i=0,j=joff1,k=koff; i<runs[r].length; i++,j++,k++) {
      data[j] = packedData[k];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    int joff2 = l*nfftz + runs[r1].z; // k >= 0
    for (int i=0,j=joff2,k=koff; i<runs[r1].length; i++,j++,k++) {
      data[j] = packedData[k];
    }//endfor
    koff += runs[r1].length;

#ifndef _BZERO_METH_
    int joff3 = joff2+runs[r1].length;
    for(int j=joff3;j<joff1;j++){data[j]=0.0;}
#endif

#ifdef CMK_VERSION_BLUEGENE
    if(r % 40==0){
      CmiNetworkProgress();
    }//endif
#endif

  }//endfor

  CkAssert(numPoints == koff);

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

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

#ifdef CMK_VERSION_BLUEGENE
    if(r % 40==0){CmiNetworkProgress();}
#endif

  }//endfor

  CkAssert(numPoints == koff);

#ifdef CMK_VERSION_BLUEGENE
  CmiNetworkProgress();
#endif

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// non-local : Gchare : data(gx,gy,z) -> data(gx,gy,gz) : backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doNlFFTRtoG_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,
                                  int numLines, int numRuns, RunDescriptor *runs, 
                                  int nfftz){
//==============================================================================
// FFT in expanded form

  fft_split(
          &bwdZPlanNL,             // Z-direction backward plan
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
void FFTcache::doNlFFTGtoR_Gchare(complex *data_in,complex *data_out,
          int numFull, int numPoints,int numLines, int numRuns, 
          RunDescriptor *runs, int nfftz){
//==============================================================================
// Expand for ffting : Trickery

  expandGSpace(data_out,data_in,runs,numRuns,numFull,numPoints,nfftz);//data_out is expanded

//==============================================================================
// FFT in expanded form

  fft_split(
          &fwdZPlanNL,               // Z-direction forward plan
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
// Hartree : Gchare : data(gx,gy,gz) -> data(gx,gy,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doHartFFTGtoR_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,int numLines, int numRuns, 
                                  RunDescriptor *runs, int nfftz){
//==============================================================================
// Expand for ffting : Trickery

  expandGSpace(data_out,data_in,runs,numRuns,numFull,numPoints,nfftz);//data_out is expanded

//==============================================================================
// FFT in expanded form

  fft_split(
          &fwdZPlanRho,             // Z-direction forward plan 
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
   	     &bwdXPlanNL,             // backward plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  if(bwdXPlanNL.option==0){
    for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along y

  fft_split(
            &bwdYPlanNL,            // backward plan
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
        &fwdYPlanNL,               // y-plan for NL Ees method
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	stride,                       // stride betwen elements (x is inner)
  	1,                            // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
  );

  // fftw only gives you one sign for real to complex : so do it yourself
  if(fwdXPlanNL.option==0){
    for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(
              &fwdXPlanNL,               // x-plan for NL Ees method
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
          &bwdZPlanEext,           // Z-direction backward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data,   //input data
	  1,                      //stride
	  nfftz,                  //distance between z-data sets
	  NULL, 0, 0,             // input is ouput
          config.fftprogresssplit // split parameter
         ); 

//==============================================================================
// pack for computing

  getCacheMem("doEextFFTRtoG_Gchare");
  packGSpace(data,tmpData,runs,numRuns,numFull,numPoints,nfftz); // tmpdata returns packed
  memcpy(data,tmpData,sizeof(complex)*numPoints); // data is expanded; cp in tmpdata; 
  freeCacheMem("doEextFFTRtoG_Gchare");

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//=============================================================================
// Eext : Gchare : data(gx,gy,gz) -> data(gx,gy,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTGtoR_Gchare(complex *data_in, complex *data_out, 
                                    int numFull, int numPoints,int numLines, 
                                    int numRuns, RunDescriptor *runs,
                                    int nfftz){
//==============================================================================
// expand for ffting

  expandGSpace(data_out,data_in,runs,numRuns,numFull,numPoints,nfftz);// data_out expanded

//==============================================================================
// FFT

  fft_split(
          &fwdZPlanEext,           // Z-direction forward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data_out,//output data
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
   	     &bwdXPlanEext,             // backward plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;

  if(bwdXPlanEext.option==0){
   for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along y

  fft_split(
            &bwdYPlanEext,          // backward plan
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
// Eext : Rchare : data(x,y,z) -> data(gx,gy,z) : Backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTRxToGx_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                    int sizeX,int sizeY){
//==============================================================================
// FFT along x first : stride 1

  rfftwnd_real_to_complex_split(
   	     &bwdXPlanEext,             // backward plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  if(bwdXPlanEext.option==0){
    for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//------------------------------------------------------------------------------
 }//end routine
//==============================================================================



//=============================================================================
// Eext : Rchare : data(x,y,z) -> data(gx,gy,z) : Backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTRyToGy_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                    int sizeX,int sizeY){
//==============================================================================
// FFT along y : stride 1

  fft_split(
            &bwdYPlanEextS,            // backward plan
	    nplane_x,                  // these many 1D ffts
	    (fftw_complex *)dataC,     // data set
            1,                         // stride
            sizeY,                     // spacing between data sets
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
        &fwdYPlanEext,             // y-plan for NL Ees method
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	stride,                       // stride betwen elements (x is inner)
  	1,                            // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
       );

  // fftw only gives you one sign for real to complex : so do it yourself
  if(fwdXPlanEext.option==0){
   for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(
              &fwdXPlanEext,             // x-plan for NL Ees method
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
// Eext : Rchare : data(gx,gy,z) -> data(x,y,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTGxToRx_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                      int sizeX,int sizeY){
//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  int stride = sizeX/2+1;
  if(fwdXPlanEext.option==0){
    for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

  rfftwnd_complex_to_real_split(
              &fwdXPlanEext,             // x-plan for NL Ees method
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
// Eext : Rchare : data(gx,gy,z) -> data(x,y,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doEextFFTGyToRy_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                    int sizeX,int sizeY){
//==============================================================================
// FFT along Y direction : Y moves with stride 1 through memory
//                       : nplane_x is spherical cutoff <= sizeX/2+1

  fft_split(
        &fwdYPlanEextS,               // y-plan for NL Ees method
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	1,                            // stride betwen elements (x is inner)
  	sizeY,                        // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
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
void FFTcache::doStpFFTRtoG_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,
                                  int numLines, int numRuns, RunDescriptor *runs, 
                                  int nfftz){
//==============================================================================
// FFT in expanded form

  fft_split(
          &bwdZPlanRho,           // Z-direction backward plan
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
void FFTcache::doStpFFTGtoR_Gchare(complex *data_in,complex *data_out,
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
          &fwdZPlanRho,            // Z-direction forward plan
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
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doStpFFTGtoR_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                   int sizeX,int sizeY){
//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory
//                       : nplane_x is spherical cutoff <= sizeX/2+1

  int stride = sizeX/2+1;
  fft_split(
        &fwdYPlanRho,                 // y-plan 
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	stride,                       // stride betwen elements (x is inner)
  	1,                            // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
       );

  // fftw only gives you one sign for real to complex : so do it yourself
  if(fwdXPlanRho.option==0){
    for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(
              &fwdXPlanRho,                 // x-plan 
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
// 
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doStpFFTRtoG_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                    int sizeX,int sizeY){
//==============================================================================
// FFT along x first

  rfftwnd_real_to_complex_split(
   	     &bwdXPlanRho,               // backward X plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  if(bwdXPlanRho.option==0){
    for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along y

  fft_split(
            &bwdYPlanRho,              // backward Y plan
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
// StatePlane : Gchare : data(gx,gy,gz) -> data(gx,gy,z) : forward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTGtoR_Gchare(complex *data_in,complex *data_out,
                                    int numFull, int numPoints,
                                    int numLines, int numRuns, RunDescriptor *runs, 
                                    int nfftz,int iexpand){
//==============================================================================
// Expand for ffting : data_out is expanded

  if(iexpand==1){
    expandGSpace(data_out,data_in,runs,numRuns,numFull,numPoints,nfftz);
  }//endif

//==============================================================================
// FFT in expanded form

  fft_split(
          &fwdZPlanRho,            // Z-direction forward plan
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
 


//=============================================================================
// StatePlane : Gchare : data(gx,gy,z) -> data(gx,gy,gz) : backward
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTRtoG_Gchare(complex *data_in,complex *data_out,
                                  int numFull, int numPoints,
                                  int numLines, int numRuns, RunDescriptor *runs, 
                                  int nfftz, int ipack){
//==============================================================================
// FFT in expanded form

  fft_split(
          &bwdZPlanRho,           // Z-direction backward plan
          numLines,               // # of ffts : one for every line of z in the chare
	  (fftw_complex *)data_in,//input data
	  1,                      //stride
	  nfftz,                  //distance between z-data sets
	  NULL, 0, 0,             // input is ouput
          config.fftprogresssplit // split parameter
         ); 

//==============================================================================
// Pack for computing if necessary : data_in is expanded : dataout is contracted

  if(ipack==1){
    packGSpace(data_in,data_out,runs,numRuns,numFull,numPoints,nfftz);
  }//endif

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================





//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTRtoG_Rchare(complex *dataC,double *dataR,int nplane_x, 
				   int sizeX,int sizeY)
//==============================================================================
   {//begin routine 
//==============================================================================
// FFT along x first

  rfftwnd_real_to_complex_split(
   	     &bwdXPlanRho,               // backward X plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  if(bwdXPlanRho.option==0){
    for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along y

  fft_split(
            &bwdYPlanRho,              // backward Y plan
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


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTRxToGx_Rchare(complex *dataC,double *dataR,int nplane_x, 
  				     int sizeX,int sizeY)
//==============================================================================
   {//begin routine 
//==============================================================================
// FFT along x first

  rfftwnd_real_to_complex_split(
   	     &bwdXPlanRho,               // backward X plan
	     sizeY,                      // these many 1D ffts
	     (fftw_real *)dataR,         // data set
             1,                          // stride
             (sizeX+2),                  // spacing between data sets in doubles
  	     NULL,0,0,                   // input array is output array
             config.fftprogresssplitReal // 
           );            

  int stride = sizeX/2+1;
  if(bwdXPlanRho.option==0){
    for (int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//------------------------------------------------------------------------------
   }//end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTRyToGy_Rchare(complex *dataC,double *dataR,int nplane_x, 
				   int sizeX,int sizeY)
//==============================================================================
   {//begin routine 
//==============================================================================
// FFT along y

  fft_split(
            &bwdYPlanRhoS,             // backward Y plan
	    nplane_x,                  // these many 1D ffts
	    (fftw_complex *)dataC,     // data set
            1,                         // stride
            sizeY,                     // spacing between data sets
            NULL,0,0,                  // input is output
            config.fftprogresssplitReal// progress splitting for BG/L
           );

//------------------------------------------------------------------------------
   }//end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTGtoR_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                   int sizeX,int sizeY){
//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory
//                       : nplane_x is spherical cutoff <= sizeX/2+1

  int stride = sizeX/2+1;
  fft_split(
        &fwdYPlanRho,                 // y-plan
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	stride,                       // stride betwen elements (x is inner)
  	1,                            // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
       );

  // fftw only gives you one sign for real to complex : so do it yourself
  if(fwdXPlanRho.option==0){
    for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(
              &fwdXPlanRho,                 // x-plan 
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


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTGyToRy_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                     int sizeX,int sizeY){
//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory
//                       : nplane_x is spherical cutoff <= sizeX/2+1

  fft_split(
        &fwdYPlanRhoS,                // y-plan 
	nplane_x,                     // how many < sizeX/2 + 1
        (fftw_complex *)(dataC),      //input data
 	1,                            // stride betwen elements (x is inner)
  	sizeY,                        // array separation (nffty elements)
        NULL,0,0,                     // output data is input data
        config.fftprogresssplitReal   // splitting parameter
       );

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void FFTcache::doRhoFFTGxToRx_Rchare(complex *dataC,double *dataR,int nplane_x, 
                                     int sizeX,int sizeY){
//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  int stride = sizeX/2+1;
  if(fwdXPlanRho.option==0){
    for(int i=0;i<stride*sizeY;i++){dataC[i].im = -dataC[i].im;}
  }//endif

  rfftwnd_complex_to_real_split(
              &fwdXPlanRho,                 // x-plan 
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


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/**
 * split up an fft call into multiple invocations with cmiprogress
 * calls between them.
 */
//============================================================================
void fft_split(FFTplanHolder *fftplanholder, int howmany, 
               fftw_complex *in,     int istride_in,    int idist_in, 
               fftw_complex *out_in, int ostride_in, int odist_in, int split)
//============================================================================
  {// begin routine 
//============================================================================

  fftw_plan plan    = fftplanholder->fftwPlan;
  int iopt          = fftplanholder->option;
  int isign         = fftplanholder->isign;
  int nfft          = fftplanholder->nfft;
  int nwork1        = fftplanholder->nwork1;
  int nwork2        = fftplanholder->nwork2;
  double scale      = fftplanholder->scale;
  double *work1     = fftplanholder->work1;
  double *work2     = fftplanholder->work2;

  int zero          = 0;
  int istride       = istride_in;
  int idist         = idist_in;
  int ostride       = ostride_in;
  int odist         = odist_in;
  fftw_complex *out = out_in;
  if(iopt==1){
    ostride = fftplanholder->ostride;
    odist   = fftplanholder->odist;
    out     = (out==NULL) ? in : out ;
  }//endif
  
  int thismany  = split;
  int inleft    = howmany;
  int numsplits = howmany/split;
  if(numsplits<=0){numsplits=1;}
  if(howmany>split && howmany%split!=0){numsplits++;}

//============================================================================

  int inoff=0;
  int outoff=0;

  for(int i=0;i<numsplits;i++){
      thismany = split;
      if(inleft<split){thismany=inleft;}
#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(in[inoff]) % 16 ==0);
      CkAssert((unsigned int) &(out[outoff]) % 16 ==0);
#endif
      fftw_complex *dummy = (out==NULL) ? NULL : &(out[outoff]) ;
      switch(iopt){
        case 0:
                fftw(
                    plan,          // da plan
                    thismany,      // how many 
                    &(in[inoff]),  // input data
                    istride,       // stride betwen elements 
                    idist,         // array separation
                    dummy,         // output data (null if inplace)
                    ostride,       // stride betwen elements  (0 inplace)
                    odist);       // array separation        (0 inplace)
   	        break;
        case 1: dcftWrap(&zero,(complex *)&(in[inoff]),&istride,&idist,
                               (complex *)dummy,&ostride,&odist,
                         &nfft,&thismany,&isign,&scale,work1,&nwork1,work2, &nwork2);
                break;
        default : CkAbort("impossible fft iopt");  break;
      }//end switch
      inoff  += idist*(thismany);  
      outoff += odist*(thismany);
      inleft -= thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
  }//endfor

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
void rfftwnd_complex_to_real_split(RFFTplanHolder *rfftplanholder, int howmany, 
         fftw_complex *in,    int istride_in,    int idist_in, 
         fftw_real   *out_in, int ostride_in, int odist_in, int split)
//============================================================================
   {//begin routine
//============================================================================

  rfftwnd_plan plan = rfftplanholder->rfftwPlan;
  int iopt          = rfftplanholder->option;
  int isign         = rfftplanholder->isign;
  int nfft          = rfftplanholder->nfft;
  int nwork1        = rfftplanholder->nwork1;
  int nwork2        = rfftplanholder->nwork2;
  double scale      = rfftplanholder->scale;
  double *work1     = rfftplanholder->work1;
  double *work2     = rfftplanholder->work2;

  int zero          = 0;
  int istride       = istride_in;
  int idist         = idist_in;
  int ostride       = ostride_in;
  int odist         = odist_in;
  fftw_real *out    = out_in;
  if(iopt==1){
    ostride = rfftplanholder->ostride;
    odist   = rfftplanholder->odist;
    if(out==NULL){
      out     = reinterpret_cast<double*> (in);
    }//endif
  }//endif

  int thismany  = split;
  int inleft    = howmany;
  int numsplits = howmany/split;
  if(numsplits<=0){numsplits=1;}
  if(howmany>split && howmany%split!=0){numsplits++;}

//============================================================================  

  int inoff=0;
  int outoff=0;

  for(int i=0;i<numsplits;i++){
      thismany=split;
      if(inleft<split){thismany=inleft;}
#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(in[inoff]) % 16 ==0);
      CkAssert((unsigned int) &(out[outoff]) % 16 ==0);
#endif
      fftw_real *dummy = (out==NULL) ? NULL : &(out[outoff]);
      switch(iopt){
        case 0:
          rfftwnd_complex_to_real(
              plan,          // da plan
              thismany,      // how many 
              &(in[inoff]),  // input data
              istride,       // stride betwen elements 
              idist,         // array separation
              dummy,         // output data (null if inplace)
              ostride,       // stride betwen elements  (0 inplace)
              odist);       // array separation        (0 inplace)
	  break;
        case 1: dcrftWrap(&zero,(complex *)&(in[inoff]),&idist,(double *)dummy,&odist,  
                          &nfft,&thismany,&isign,&scale,work1,&nwork1,work2,&nwork2); break;
        default : CkAbort("impossible fft iopt"); break;
      }//end switch
      inoff  += idist*(thismany);  
      outoff += odist*(thismany);
      inleft -= thismany;      
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
void rfftwnd_real_to_complex_split(RFFTplanHolder *rfftplanholder, int howmany, 
        fftw_real *in,        int istride_in,    int idist_in, 
        fftw_complex *out_in, int ostride_in, int odist_in, int split)
//============================================================================
   {// begin routine
//============================================================================

  rfftwnd_plan plan = rfftplanholder->rfftwPlan;
  int iopt          = rfftplanholder->option;
  int isign         = rfftplanholder->isign;
  int nfft          = rfftplanholder->nfft;
  int nwork1        = rfftplanholder->nwork1;
  int nwork2        = rfftplanholder->nwork2;
  double scale      = rfftplanholder->scale;
  double *work1     = rfftplanholder->work1;
  double *work2     = rfftplanholder->work2;

  int zero          = 0;
  int istride       = istride_in;
  int idist         = idist_in;
  int ostride       = ostride_in;
  int odist         = odist_in;
  fftw_complex *out = out_in;
  if(iopt==1){
    ostride = rfftplanholder->ostride;
    odist   = rfftplanholder->odist;
    if(out==NULL){
      out     = reinterpret_cast<fftw_complex*> (in);
    }//endif
  }//endif

  int thismany  = split;
  int inleft    = howmany;
  int numsplits = howmany/split;
  if(numsplits<=0){numsplits=1;}
  if(howmany>split && howmany%split!=0){numsplits++;}

//============================================================================    

  int inoff  = 0;
  int outoff = 0;
  //  CkPrintf("inside rfft plan %p howmany %d in %p istride %d idist %d out %p ostride %d odist %d split %d dummy %p\n",plan,  howmany, in, istride_in, idist_in, out_in,  ostride_in, odist_in, split,(out==NULL) ? NULL : &(out[outoff]));
  for(int i=0;i<numsplits;i++){
      thismany=split;
      if(inleft<split){thismany=inleft;}
#ifdef TEST_ALIGN
      CkAssert((unsigned int) &(in[inoff]) % 16 ==0);
      CkAssert((unsigned int) &(out[outoff]) % 16 ==0);
#endif
      fftw_complex *dummy = (out==NULL) ? NULL : &(out[outoff]);
      switch(iopt){
        case 0:
          rfftwnd_real_to_complex(
              plan,         // da plan
              thismany,     // how many 
              &(in[inoff]), // input data
              istride,      // stride betwen elements 
              idist,        // array separation
              dummy,        // output data (null if inplace)
              ostride,      // stride betwen elements  (0 inplace)
              odist);       // array separation (0 inplace)
	  break;
        case 1: drcftWrap(&zero,(double *)&(in[inoff]),&idist,(complex *)dummy,&odist,
                          &nfft,&thismany,&isign,&scale,work1,&nwork1,work2,&nwork2); break;
        default : CkAbort("impossible fft iopt");  break;
      }//end switch
      inoff  += idist*(thismany);  
      outoff += odist*(thismany);
      inleft -= thismany;      
#ifdef CMK_VERSION_BLUEGENE
      CmiNetworkProgress();
#endif
    }// endfor

//----------------------------------------------------------------------------
  }//end routine
//============================================================================



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void initFFTholder(FFTplanHolder *plan, int *iopt,int *nwork1,int *nwork2, double *scale,
                   int *isign, int *nfft, int *stride, int *skip){
//============================================================================

  plan->option  = iopt[0];
  plan->nwork1  = nwork1[0];
  plan->nwork2  = nwork2[0];
  plan->scale   = scale[0];
  plan->isign   = isign[0];
  plan->nfft    = nfft[0];
  plan->ostride = stride[0];
  plan->odist   = skip[0];
  plan->work1   = NULL;
  plan->work2   = NULL;

  if(iopt[0] == 1){
    int unit = 1;
    complex x[2];
    double *work1 = (double*) fftw_malloc(nwork1[0]*sizeof(double)); 
    double *work2 = (double*) fftw_malloc(nwork2[0]*sizeof(double)); 
    dcftWrap(&unit,x,stride,skip,x,stride,skip,nfft,&unit,isign,scale,
             work1,nwork1,work2,nwork2);
    plan->work1   = work1;
    plan->work2   = work2;
  }//endif

//----------------------------------------------------------------------------
  }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void initRCFFTholder(RFFTplanHolder *plan, int *iopt,int *nwork1,int *nwork2, double *scale,
                    int *isign, int *nfft, int *skipR, int *skipC){
//============================================================================

  plan->option  = iopt[0];
  plan->nwork1  = nwork1[0];
  plan->nwork2  = nwork2[0];
  plan->scale   = scale[0];
  plan->isign   = isign[0];
  plan->nfft    = nfft[0];
  plan->ostride = 1;
  plan->odist   = skipC[0];
  plan->work1   = NULL;
  plan->work2   = NULL;

  if(iopt[0] == 1){
    int unit = 1;
    complex xc[2]; double xr[2];
    double *work1 = (double*) fftw_malloc(nwork1[0]*sizeof(double));
    double *work2 = (double*) fftw_malloc(nwork2[0]*sizeof(double));
    drcftWrap(&unit,xr,skipR,xc,skipC,nfft,&unit,isign,scale,work1,nwork1,work2,nwork2);
    plan->work1   = work1;
    plan->work2   = work2;
  }//endif

//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void initCRFFTholder(RFFTplanHolder *plan, int *iopt,int *nwork1,int *nwork2, double *scale,
                    int *isign, int *nfft, int *skipR, int *skipC){
//============================================================================

  plan->option  = iopt[0];
  plan->nwork1  = nwork1[0];
  plan->nwork2  = nwork2[0];
  plan->scale   = scale[0];
  plan->isign   = isign[0];
  plan->nfft    = nfft[0];
  plan->ostride = 1;
  plan->odist   = skipR[0];
  plan->work1   = NULL;
  plan->work2   = NULL;

  if(iopt[0] == 1){
    int unit = 1;
    complex xc[2]; double xr[2];
    double *work1 = (double*) fftw_malloc(nwork1[0]*sizeof(double));
    double *work2 = (double*) fftw_malloc(nwork2[0]*sizeof(double));
    dcrftWrap(&unit,xc,skipC,xr,skipR,nfft,&unit,isign,scale,work1,nwork1,work2,nwork2);
    plan->work1   = work1;
    plan->work2   = work2;
  }//endif

//----------------------------------------------------------------------------
  }//end routine
//============================================================================
