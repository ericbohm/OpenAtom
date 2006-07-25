//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file rhoslab.C
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
RhoRealSlab::~RhoRealSlab(){

    fftw_free(Vks); //
    fftw_free(density);
    complex *dummy;
	dummy  = reinterpret_cast<complex*> (doFFTonThis);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (rhoIRX);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (rhoIRY);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (rhoIRZ);
	fftw_free((fftw_complex *)dummy);
	dummy  = reinterpret_cast<complex*> (gradientCorrection);
	fftw_free((fftw_complex *)dummy);
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the RealSpaceDensity Constructor */
//==============================================================================
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim,
                     int numRealSpace, int numRhoG, int myIndexX, int myIndexY)
//==============================================================================
   {//begin routine
//==============================================================================

	rho_rs->sizeX = xdim;
	rho_rs->sizeY = ydim;
	rho_rs->sizeZ = zdim;

	rho_rs->xdim = rho_rs->sizeX;
	rho_rs->ydim = 1;
	rho_rs->zdim = rho_rs->sizeZ;

	rho_rs->startx = 0;
	rho_rs->starty = 0;
	rho_rs->startz = 0; 

	rho_rs->size     = (rho_rs->sizeX+2) * (rho_rs->sizeZ);
	rho_rs->trueSize = (rho_rs->sizeX) * (rho_rs->sizeZ);
        int sizenow      = rho_rs->size;

	rho_rs->Vks     =  (double *)fftw_malloc(sizenow*sizeof(double));
	bzero(rho_rs->Vks, sizenow*sizeof(double));
	rho_rs->density =  (double *)fftw_malloc(sizenow*sizeof(double));
	bzero(rho_rs->density, sizenow*sizeof(double));
        int csizenow    = sizenow/2;

        complex *dummy;
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->doFFTonThis = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->rhoIRX      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->rhoIRY      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->rhoIRZ      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
	rho_rs->gradientCorrection = reinterpret_cast<double*> (dummy);

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoGSlab::~RhoGSlab()
{
  if(Rho!=NULL)
    fftw_free(Rho);
  if(divRhoX!=NULL)
    fftw_free( divRhoX);
  if(divRhoY!=NULL)
    fftw_free( divRhoY);
  if(divRhoZ!=NULL)
    fftw_free( divRhoZ);
  if(packedRho!=NULL)
    fftw_free( packedRho);
  if(packedVks!=NULL)
    fftw_free( packedVks);
  if(runs!=NULL)
      delete [] runs;
  if(k_x!=NULL)
      fftw_free( k_x);
  if(k_y!=NULL)
      fftw_free( k_y);
  if(k_z!=NULL)
      fftw_free( k_z);
}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoGSlab::pup(PUP::er &p) {
  // local bools so we only bother with non null objects
  bool RhoMake=false;
  bool divRhoXMake=false;
  bool divRhoYMake=false;
  bool divRhoZMake=false;
  bool packedRhoMake=false;
  bool packedVksMake=false;
  bool VksMake=false;
  p|sizeX;
  p|sizeY;
  p|sizeZ;
  p|runsToBeSent;
  p|numRuns;
  p|numLines;
  p|numFull;
  p|numPoints;
  p|ehart_ret;
  p|eext_ret;
  p|ewd_ret;
  p|size;
  p|nPacked;
  if(!p.isUnpacking())
    {// create flags for each array in pup
      RhoMake= (Rho!=NULL) ? true :false;
      divRhoXMake= (divRhoX!=NULL) ? true :false;
      divRhoYMake= (divRhoY!=NULL) ? true :false;
      divRhoZMake= (divRhoZ!=NULL) ? true :false;
      packedRhoMake= (packedRho!=NULL) ? true :false;
      packedVksMake= (packedVks!=NULL) ? true :false;
      VksMake= (Vks!=NULL) ? true :false;
    }
  p|RhoMake;
  p|divRhoXMake;
  p|divRhoYMake;
  p|divRhoZMake;
  p|packedRhoMake;
  p|packedVksMake;
  p|VksMake;
  if(p.isUnpacking())
    {
      if(RhoMake)
	Rho       = (complex *)fftw_malloc(numFull*sizeof(complex));
      else
	Rho = NULL;
      if(divRhoXMake)
	divRhoX   = (complex *)fftw_malloc(numFull*sizeof(complex));
      else
	divRhoX   = NULL;
      if(divRhoYMake)
	divRhoY   = (complex *)fftw_malloc(numFull*sizeof(complex));
      else 
	divRhoY   = NULL;
      if(divRhoZMake)
	divRhoZ   = (complex *)fftw_malloc(numFull*sizeof(complex));
      else
	divRhoZ   = NULL;
      if(packedRhoMake)
	packedRho = (complex *)fftw_malloc(nPacked*sizeof(complex));
      else
	packedRho  = NULL;
      if(packedVksMake)
	packedVks = (complex *)fftw_malloc(nPacked*sizeof(complex));
      else
	packedVks = NULL;

      if(VksMake)
	Vks=(complex *)fftw_malloc(numFull*sizeof(complex));
      else
	Vks=NULL;
      runs     = new RunDescriptor[numRuns];
      k_x = (int *)fftw_malloc(numPoints*sizeof(int));
      k_y = (int *)fftw_malloc(numPoints*sizeof(int));
      k_z = (int *)fftw_malloc(numPoints*sizeof(int));

    }
  if(RhoMake)
    PUParray(p,Rho,numFull);
  if(divRhoXMake)
    PUParray(p,divRhoX,numFull);
  if(divRhoYMake)
    PUParray(p,divRhoY,numFull);
  if(divRhoZMake)
    PUParray(p, divRhoZ,numFull);
  if(packedRhoMake)
    PUParray(p, packedRho,nPacked);
  if(packedVksMake)
    PUParray(p, packedVks,nPacked);
  if(VksMake)
    PUParray(p,Vks,numFull);
  PUParray(p,k_x,numPoints);
  PUParray(p,k_y,numPoints);
  PUParray(p,k_z,numPoints);
  PUParray(p,runs,numRuns);

}
//==============================================================================


//==============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void RhoGSlab::divRhoGdot(double *hmati, double tpi,complex *tmpRho){

//==============================================================================

  int nfftz = sizeZ;

//==============================================================================
// 

  bzero(divRhoX,sizeof(complex)*numFull);  
  bzero(divRhoY,sizeof(complex)*numFull);  
  bzero(divRhoZ,sizeof(complex)*numFull);  
  double gx,gy,gz;

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      gx = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      complex tmp = (tmpRho[k].multiplyByi())*(-1.0);
      divRhoX[j] = tmp*gx;
      divRhoY[j] = tmp*gy;
      divRhoZ[j] = tmp*gz;
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      gx = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      complex tmp = (tmpRho[k].multiplyByi())*(-1.0);
      divRhoX[j] = tmp*gx;
      divRhoY[j] = tmp*gy;
      divRhoZ[j] = tmp*gz;
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


//==============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void RhoGSlab::createWhiteByrd(double *hmati, double tpi){

//==============================================================================

  int nfftz = sizeZ;
  double gx,gy,gz;
  complex tmp;

//==============================================================================

  complex *whitebyrd = divRhoX; // zeroing done carefully inside loop
                                // so that we can save memory by reusing divRhoX
  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff1 = l*nfftz + runs[r].z;
    for (int i=0,j=joff1,k=koff; i<runs[r].length; i++,j++,k++) {
      gx  = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy  = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz  = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      tmp = divRhoX[j]*gx + divRhoY[j]*gy + divRhoZ[j]*gz;
      whitebyrd[j] = tmp.multiplyByi()*(-1.0); 
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    int joff2 = l*nfftz + runs[r1].z;
    for (int i=0,j=joff2,k=koff; i<runs[r1].length; i++,j++,k++) {
      gx  = tpi*(k_x[k]*hmati[1] + k_y[k]*hmati[2] + k_z[k]*hmati[3]);
      gy  = tpi*(k_x[k]*hmati[4] + k_y[k]*hmati[5] + k_z[k]*hmati[6]);
      gz  = tpi*(k_x[k]*hmati[7] + k_y[k]*hmati[8] + k_z[k]*hmati[9]);
      tmp = divRhoX[j]*gx + divRhoY[j]*gy + divRhoZ[j]*gz;
      whitebyrd[j] = tmp.multiplyByi()*(-1.0); 
    }//endfor
    koff += runs[r1].length;

    int joff3 = joff2+runs[r1].length;
    for(int j=joff3;j<joff1;j++){whitebyrd[j]=0.0;}

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


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/*
 * In G-space our representation uses wrapping of coordinates.
 * When reading from the files, if x < 0, x += sizeX;
 */
//==============================================================================

void RhoGSlab::setKVectors(int *n){

//======================================================================
// Construct the k-vectors

  k_x = (int *)fftw_malloc(numPoints*sizeof(int));
  k_y = (int *)fftw_malloc(numPoints*sizeof(int));
  k_z = (int *)fftw_malloc(numPoints*sizeof(int));
  
  int r, i, dataCovered = 0;
  int x, y, z;
  for (r = 0; r < numRuns; r++) { // 2*number of lines z
    x = runs[r].x;
    if (x > sizeX/2) x -= sizeX;
    y = runs[r].y;
    if (y > sizeY/2) y -= sizeY;
    z = runs[r].z;
    if (z > sizeZ/2) z -= sizeZ;

    for (i = 0; i < runs[r].length; i++) { //pts in lines of z
      k_x[dataCovered] = x;
      k_y[dataCovered] = y;
      k_z[dataCovered] = (z+i);
      dataCovered++;
    }//endfor
  }//endfor

  CkAssert(dataCovered == numPoints);

//==============================================================================
// Set the return values

  *n    = numPoints;
  
//==============================================================================
  }//end routine
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoRealSlab::pup(PUP::er &p) {
  p|sizeZ;
  p|sizeY;
  p|sizeZ;
  p|exc_ret;
  p|muxc_ret;
  p|exc_gga_ret;
  p|size;
  p|trueSize;
  p|xdim;
  p|ydim;
  p|zdim;
  p|startx;
  p|starty;
  p|startz; 
  int csize    = size/2;
  if(p.isUnpacking())
    {
	Vks     =  (double *)fftw_malloc(size*sizeof(double));
	bzero(Vks, size*sizeof(double));
	density =  (double *)fftw_malloc(size*sizeof(double));
	bzero(density, size*sizeof(double));
        complex *dummy;
	dummy               = (complex*) fftw_malloc(csize*sizeof(complex));
	doFFTonThis = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csize*sizeof(complex));
	rhoIRX      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csize*sizeof(complex));
	rhoIRY      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csize*sizeof(complex));
	rhoIRZ      = reinterpret_cast<double*> (dummy);
	dummy               = (complex*) fftw_malloc(csize*sizeof(complex));
	gradientCorrection = reinterpret_cast<double*> (dummy);
    }
  PUParray(p,Vks, size);
  PUParray(p,density,size);
  PUParray(p, doFFTonThis,size);
  PUParray(p, rhoIRX,size);
  PUParray(p, rhoIRY,size);
  PUParray(p, rhoIRZ,size);
  PUParray(p, gradientCorrection,size);
}
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoRealSlab::doFwFFTGtoR(int iopt,double probScale){
//==============================================================================

   FFTcache *fsc             = fftCacheProxy.ckLocalBranch();
   fftw_plan    fwdZ1DdpPlan = fsc->fwdZ1DdpPlan; 
   rfftwnd_plan fwdX1DdpPlan = fsc->fwdX1DdpPlan; 

   CPcharmParaInfo *sim   = (scProxy.ckLocalBranch ())->cpcharmParaInfo;      
   int nplane_x           = sim->nplane_rho_x;
   int pSize              = (sizeX+2) * sizeY;
   int stride             = sizeX/2+1;

   double *data;
   switch(iopt){
     case 0: data = Vks; break;
     case 1: data = rhoIRX; break;
     case 2: data = rhoIRY; break;
     case 3: data = rhoIRZ; break;
     case 4: data = doFFTonThis; break;
   }//endif
  complex *planeArr = reinterpret_cast<complex*> (data);

//==============================================================================
// FFT along Y direction : Y moves with stride sizex/2+1 through memory

   fft_split(fwdZ1DdpPlan,    // y-plan (label lies)
	nplane_x,        // how many < sizeX/2 + 1
        (fftw_complex *)(planeArr),//input data
 	stride,          // stride betwen elements (x is inner)
  	1,               // array separation (nffty elements)
        NULL,0,0,        // output data is input data
        config.fftprogresssplitReal);

#ifdef CMK_VERSION_BLUEGENE
   CmiNetworkProgress();
#endif

  // fftw only gives you one sign for real to complex : so do it yourself
  for(int i=0;i<stride*sizeY;i++){planeArr[i].im = -planeArr[i].im;}

//==============================================================================
// FFT along X direction : X moves with stride 1 through memory

  rfftwnd_complex_to_real_split(fwdX1DdpPlan,
		          sizeY,    // how many
			  (fftw_complex *)planeArr, 
                          1,               // stride (x is inner)
                          stride,          // array separation
  			  NULL,0,0,       // output array is real 
                          config.fftprogresssplitReal); 

//==============================================================================
// copy out to a nice stride in memory

  double *temp = gradientCorrection;
  memcpy(temp,data,pSize*sizeof(double));

  for(int i=0,i2=0;i<sizeY;i++,i2+=2){
    for(int j=i*sizeX;j<(i+1)*sizeX;j++){
      data[j] = temp[(j+i2)]*probScale;
    }//endfor
  }//endfor

//==============================================================================
  }//end if
//==============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::uPackAndScale(double *uPackData, double *PackData, 
                                double scale)
//============================================================================
   {//begin routine
//============================================================================

   for(int i=0,i2=0;i<sizeZ;i++,i2+=2){
     for(int j=i*sizeX;j<(i+1)*sizeX;j++){
       uPackData[(j+i2)] = PackData[j]*scale;
     }//endfor
     uPackData[sizeX]   = 0.0;
     uPackData[sizeX+1] = 0.0;
   }//endfor

//============================================================================
   }
//============================================================================

