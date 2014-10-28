/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file rhoSlab.C
 * Add functions to allow application programmers to initialize these and
 * the corresponding functions in Charm++ to call these functions with
 * appropriate parameters 
 */
//==============================================================================

#include "utility/util.h"
#include "debug_flags.h"
#include <cmath>
#include "main/cpaimd.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_ctrl/CP_State_Plane.h"

//==============================================================================

extern CkVec <CProxy_FFTcache> UfftCacheProxy;
extern Config config;
extern int nstates;
extern int sizeX;

//==============================================================================
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoRealSlab::~RhoRealSlab(){

  fftw_free(densityC);
  fftw_free(VksC);
  fftw_free(rhoIRXC);
  fftw_free(rhoIRYC);
  fftw_free(rhoIRZC);
  fftw_free(VksHartC);
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the RealSpaceDensity Constructor */
//==============================================================================
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim, 
    int xdimA, int ydimA,  int myIndexX, int myIndexY,
    int rhoRsubplanes)
  //==============================================================================
{//begin routine
  //==============================================================================

  rho_rs->rhoRsubplanes = rhoRsubplanes;
  rho_rs->sizeX = xdim;
  rho_rs->sizeY = ydim;
  rho_rs->sizeZ = zdim;

  rho_rs->size     = (rho_rs->sizeX+2) * (rho_rs->sizeY);
  rho_rs->trueSize = (rho_rs->sizeX) * (rho_rs->sizeY);

  int sizenow      = rho_rs->size;
  int csizenow     = sizenow/2;

  complex *dummy;
  dummy            = (complex*) fftw_malloc(csizenow*sizeof(complex));
  rho_rs->VksC     = dummy;
  rho_rs->Vks      = reinterpret_cast<double*> (dummy);

  dummy            = (complex*) fftw_malloc(csizenow*sizeof(complex));
  rho_rs->densityC = dummy;
  rho_rs->density  = reinterpret_cast<double*> (dummy); 

  dummy            = (complex*) fftw_malloc(csizenow*sizeof(complex));
  rho_rs->rhoIRXC  = dummy;
  rho_rs->rhoIRX   = reinterpret_cast<double*> (dummy);

  dummy            = (complex*) fftw_malloc(csizenow*sizeof(complex));
  rho_rs->rhoIRYC  = dummy;
  rho_rs->rhoIRY   = reinterpret_cast<double*> (dummy);

  dummy            = (complex*) fftw_malloc(csizenow*sizeof(complex));
  rho_rs->rhoIRZC  = dummy;
  rho_rs->rhoIRZ   = reinterpret_cast<double*> (dummy);

  dummy            = (complex*) fftw_malloc(csizenow*sizeof(complex));
  rho_rs->VksHartC = dummy;
  rho_rs->VksHart  = reinterpret_cast<double*> (dummy);

  // if you have an extra transpose, you need a little more memory
  // to receive messages asynchronously from other elements

  rho_rs->csizeInt = 0;
  rho_rs->rsizeInt = 0;
  if(rhoRsubplanes>1){

    csizenow = xdimA*ydimA;
    rho_rs->csizeInt = csizenow;
    rho_rs->rsizeInt = 2*csizenow;

    dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
    rho_rs->rhoIRXCint  = dummy;
    rho_rs->rhoIRXint   = reinterpret_cast<double*> (dummy);

    dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
    rho_rs->rhoIRYCint  = dummy;
    rho_rs->rhoIRYint   = reinterpret_cast<double*> (dummy);

    dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
    rho_rs->rhoIRZCint  = dummy;
    rho_rs->rhoIRZint   = reinterpret_cast<double*> (dummy);

    dummy               = (complex*) fftw_malloc(csizenow*sizeof(complex));
    rho_rs->VksHartCint = dummy;
    rho_rs->VksHartint  = reinterpret_cast<double*> (dummy);

  }//endif

  //==============================================================================
}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoGSlab::~RhoGSlab(){

  if(Rho      !=NULL)fftw_free(Rho);
  if(divRhoX  !=NULL)fftw_free(divRhoX);
  if(divRhoY  !=NULL)fftw_free(divRhoY);
  if(divRhoZ  !=NULL)fftw_free(divRhoZ);
  if(packedRho!=NULL)fftw_free(packedRho);
  if(packedVks!=NULL)fftw_free(packedVks);

  if(runs!=NULL) delete [] runs;
  if(k_x !=NULL) fftw_free(k_x);
  if(k_y !=NULL) fftw_free(k_y);
  if(k_z !=NULL) fftw_free(k_z);
  if(perdCorr!=NULL)fftw_free(perdCorr);

}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoGSlab::pup(PUP::er &p) {
  // local bools so we only bother with non null objects
  p|iperd;
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

  bool RhoMake      =false;
  bool divRhoXMake  =false;
  bool divRhoYMake  =false;
  bool divRhoZMake  =false;
  bool packedRhoMake=false;
  bool packedVksMake=false;
  bool VksMake      =false;
  if(!p.isUnpacking()){// create flags for each array in pup
    RhoMake       = (Rho      !=NULL) ? true :false;
    divRhoXMake   = (divRhoX  !=NULL) ? true :false;
    divRhoYMake   = (divRhoY  !=NULL) ? true :false;
    divRhoZMake   = (divRhoZ  !=NULL) ? true :false;
    packedRhoMake = (packedRho!=NULL) ? true :false;
    packedVksMake = (packedVks!=NULL) ? true :false;
    VksMake       = (Vks      !=NULL) ? true :false;
  }//endif

  p|RhoMake;
  p|divRhoXMake;
  p|divRhoYMake;
  p|divRhoZMake;
  p|packedRhoMake;
  p|packedVksMake;
  p|VksMake;

  if(p.isUnpacking()){
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
      Vks       = (complex *)fftw_malloc(numFull*sizeof(complex));
    else
      Vks       = NULL;
    runs = new RunDescriptor[numRuns];
    k_x  = (int *)fftw_malloc(numPoints*sizeof(int));
    k_y  = (int *)fftw_malloc(numPoints*sizeof(int));
    k_z  = (int *)fftw_malloc(numPoints*sizeof(int));
    perdCorr = NULL;
    if(iperd!=3){perdCorr = (double *)fftw_malloc(numPoints*sizeof(double));}
  }//endif : unpacking malloc

  if(RhoMake)       PUParray(p,Rho,numFull);
  if(divRhoXMake)   PUParray(p,divRhoX,numFull);
  if(divRhoYMake)   PUParray(p,divRhoY,numFull);
  if(divRhoZMake)   PUParray(p,divRhoZ,numFull);
  if(packedRhoMake) PUParray(p,packedRho,nPacked);
  if(packedVksMake) PUParray(p,packedVks,nPacked);
  if(VksMake)       PUParray(p,Vks,numFull);

  PUParray(p,runs,numRuns);
  PUParray(p,k_x,numPoints);
  PUParray(p,k_y,numPoints);
  PUParray(p,k_z,numPoints);
  if(iperd!=3){PUParray(p,perdCorr,numPoints);}

  //------------------------------------------------------------------------------
}//end intense pupping experience
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

  }//endfor

  CkAssert(numPoints == koff);

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

  }//endfor

  CkAssert(numPoints == koff);
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

  p|sizeZ;   // fftX size
  p|sizeY;   // fftYZ size
  p|sizeZ;   // fftZ size
  p|size;    // plane size for fftw : bigger
  p|trueSize;// the real plane size
  p|exc_ret; // energies
  p|muxc_ret;
  p|exc_gga_ret;
  p|rhoRsubplanes;
  p|csizeInt;
  p|rsizeInt;

  if(p.isUnpacking()){
    int csize  = size/2;
    VksC       = (complex*) fftw_malloc(csize*sizeof(complex));
    Vks        =  reinterpret_cast<double*> (VksC);
    densityC   = (complex*) fftw_malloc(csize*sizeof(complex));
    density    =  reinterpret_cast<double*> (densityC);
    rhoIRXC    = (complex*) fftw_malloc(csize*sizeof(complex));
    rhoIRX     = reinterpret_cast<double*> (rhoIRXC);
    rhoIRYC    = (complex*) fftw_malloc(csize*sizeof(complex));
    rhoIRY     = reinterpret_cast<double*> (rhoIRYC);
    rhoIRZC    = (complex*) fftw_malloc(csize*sizeof(complex));
    rhoIRZ     = reinterpret_cast<double*> (rhoIRZC);
    VksHartC   = (complex*) fftw_malloc(csize*sizeof(complex));
    VksHart    = reinterpret_cast<double*> (VksHartC);
    if(rhoRsubplanes>1){
      rhoIRXCint    = (complex*) fftw_malloc(csizeInt*sizeof(complex));
      rhoIRXint     = reinterpret_cast<double*> (rhoIRXCint);
      rhoIRYCint    = (complex*) fftw_malloc(csizeInt*sizeof(complex));
      rhoIRYint     = reinterpret_cast<double*> (rhoIRYCint);
      rhoIRZCint    = (complex*) fftw_malloc(csizeInt*sizeof(complex));
      rhoIRZint     = reinterpret_cast<double*> (rhoIRZCint);
      VksHartCint   = (complex*) fftw_malloc(csizeInt*sizeof(complex));
      VksHartint    = reinterpret_cast<double*> (VksHartCint);
    }
  }
  PUParray(p,Vks,    size);
  PUParray(p,density,size);
  PUParray(p,rhoIRX, size);
  PUParray(p,rhoIRY, size);
  PUParray(p,rhoIRZ, size);
  PUParray(p,VksHart,size);

  if(rhoRsubplanes==1){
    PUParray(p,rhoIRXint, rsizeInt);
    PUParray(p,rhoIRYint, rsizeInt);
    PUParray(p,rhoIRZint, rsizeInt);
    PUParray(p,VksHartint,rsizeInt);
  }//endif

}
//==============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::uPackScale(double *uPackData, double *PackData,double scale){

  for(int i=0;i<size;i++){uPackData[i] = PackData[i]*scale;}

}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::scale(double *uPackData,double scale){

  for(int i=0;i<size;i++){uPackData[i] *= scale;}

}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::uPackScaleShrink(double *uPackData,double *PackData,double scale){

  for(int i=0,i2=0;i<sizeY;i++,i2+=2){
    for(int j=i*sizeX;j<(i+1)*sizeX;j++){
      uPackData[j] =PackData[(j+i2)]*scale;
    }//endfor
  }//endfor

}//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::uPackShrink(double *uPackData,double *PackData){

  for(int i=0,i2=0;i<sizeY;i++,i2+=2){
    for(int j=i*sizeX;j<(i+1)*sizeX;j++){
      uPackData[j] =PackData[(j+i2)];
    }//endfor
  }//endfor

}//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void RhoRealSlab::uPackScaleGrow(double *uPackData,double *PackData,double scale){
  for(int i=0,i2=0;i<sizeY;i++,i2+=2){
    for(int j=i*sizeX;j<(i+1)*sizeX;j++){
      uPackData[(j+i2)] =PackData[j]*scale;
    }//endfor
    int k = (i+1)*sizeX+i2;
    uPackData[k]    =0.0;
    uPackData[(k+1)]=0.0;
  }//endfor

}//end routine
//============================================================================
