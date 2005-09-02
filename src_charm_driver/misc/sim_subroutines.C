//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* Add functions to allow application programmers to initialize these and */
/* the corresponding functions in Charm++ to call these functions with */
/* appropriate parameters */
//==============================================================================

#include "util.h"
#include "../../include/debug_flags.h"
#include <math.h>
#include "groups.h"
#include "cpaimd.h"
#include "sim_subroutines.h"
#include "CP_State_Plane.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"

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
GStateSlab::~GStateSlab() 
{

    if(packedPlaneData    !=NULL) delete [] packedPlaneData;
    if(packedPlaneDataTemp!=NULL) delete [] packedPlaneDataTemp;
    if(packedForceData    !=NULL) delete [] packedForceData;
    if(packedPlaneDataCG  !=NULL) delete [] packedPlaneDataCG;
    if(runs               !=NULL) delete [] runs;

    packedPlaneData     = NULL;
    packedPlaneDataTemp = NULL;
    packedForceData     = NULL;
    packedPlaneDataCG   = NULL;
    runs                = NULL;

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::pup(PUP::er &p) {
//==============================================================================
// Dont have to pup fftw plans - they reside in CPcharmParaInfoGrp.
// no they live in the fft cache group

        p|numNonZeroPlanes;
	p|numRuns;
	p|numLines;
	p|numPoints;
	p|numFull;
	p|mysizeX; 
	p|planeSize;
        p|fftReqd;
	p|S_grainSize;
	p|xdim;
	p|ydim;
	p|zdim;
	p|iplane_ind;
        p|istate_ind;
        p|eke_ret;
        p|fovlap_loc;
        p|ihave_kx0;
        p|kx0_strt;
        p|kx0_end;

	if (p.isUnpacking()) {runs = new RunDescriptor[numRuns];}
	for (int i = 0; i < numRuns; i++){runs[i].pup(p);}

	if (p.isUnpacking()) {
           packedPlaneData     = new complex[numPoints];
	   packedPlaneDataCG   = new complex[numPoints];
	   packedPlaneDataTemp = new complex[numPoints];
	   packedForceData     = new complex[numPoints];
	}//endif
	p((void *) packedPlaneData, numPoints*sizeof(complex));
	p((void *) packedPlaneDataCG, numPoints*sizeof(complex));
	p((void *) packedPlaneDataTemp, numPoints*sizeof(complex));
	p((void *) packedForceData, numPoints*sizeof(complex));

//==============================================================================
  }//end routine
//==============================================================================



//=============================================================================
// packed g-space of size numPoints is expanded to numFull =numRuns/2*nfftz
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================

void FFTcache::expandGSpace(complex* data, complex *packedPlaneData, 
                            RunDescriptor *runs, int numRuns, int numFull,
                            int numPoints, int nfftz)

//==============================================================================
  {//begin routine
//==============================================================================
// Expand out lines of z so that an FFT can be performed on them
//   The ``run'' stores each line in two parts
//      0,1,2,3 have offset, joff = r*nfftz 
//     -3,-2,-1 have offset, joff = r*nffz + nfftz-3
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0
//   Total size is nlines*nfftz where nlines=numRuns/2

  complex *origdata = packedPlaneData;
  memset(data, 0, sizeof(complex)*numFull);

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      data[j] = origdata[k];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      data[j] = origdata[k];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//------------------------------------------------------------------------------
  }//end routine
//==============================================================================


//==============================================================================
// Expanded gspace of size numFull =numRuns/2*nfftz to packed size numPoints
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::compressGSpace(const complex *expnddata, int type) {
//==============================================================================
// Contract lines of z after the backFFT has been performed.
//   The ``run'' stores each line in two parts
//      0,1,2,3 have offset, joff = r*nfftz 
//     -3,-2,-1 have offset, joff = r*nffz + nfftz-3
//      runs[r].z stores kz if kz>0 and nfftz-kz if kz<0

  complex *origdata;
  if(type == 1){origdata = packedPlaneData;}
  if(type == 2){origdata = packedForceData;}
  if(type!=2){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Unpacking invalid data, aborting %d\n",type);
    CkPrintf("Only the forces should be contracted\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  int nfftz = planeSize[1];

  int koff = 0;
  for (int r = 0,l=0; r < numRuns; r+=2,l++) {

    int joff = l*nfftz + runs[r].z;
    for (int i=0,j=joff,k=koff; i<runs[r].length; i++,j++,k++) {
      origdata[k]= expnddata[j];
    }//endfor
    koff += runs[r].length;

    int r1=r+1;
    joff = l*nfftz + runs[r1].z;
    for (int i=0,j=joff,k=koff; i<runs[r1].length; i++,j++,k++) {
      origdata[k] = expnddata[j];
    }//endfor
    koff += runs[r1].length;

  }//endfor

  CkAssert(numPoints == koff);

//==============================================================================
  }//end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::addForces(const complex *points, const int *k_x){
//==============================================================================
    int i;
    double wght;

    if(!config.doublePack){

      for(i = 0; i < numPoints; i++){
	 packedForceData[i] += points[i];
      }//endfor

    }else{

      for(i = 0; i < numPoints; i++){
        wght = (k_x[i] == 0 ? 1.0 : 2.0);
#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
	packedForceData[i] = 0.0; 
#endif
#ifdef _CP_DEBUG_VKS_ONLY_
	packedForceData[i] = 0.0; 
#endif
	packedForceData[i] *= wght; 
	packedForceData[i] += points[i];
      }//endfor

    }//endif : doublePack

//==============================================================================
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

void GStateSlab::setKVectors(int *n, int **kk_x, int **kk_y, int **kk_z){

//======================================================================
// Construct the k-vectors

  int *k_x = new int[numPoints];
  int *k_y = new int[numPoints];
  int *k_z = new int[numPoints];
  
  int r, i, dataCovered = 0;
  int x, y, z;
  for (r = 0; r < numRuns; r++) { // 2*number of lines z
    x = runs[r].x;
    if (x > sizeX/2) x -= sizeX;
    y = runs[r].y;
    if (y > planeSize[0]/2) y -= planeSize[0];
    z = runs[r].z;
    if (z > planeSize[1]/2) z -= planeSize[1];

    for (i = 0; i < runs[r].length; i++) { //pts in lines of z
      k_x[dataCovered] = x;
      k_y[dataCovered] = y;
      k_z[dataCovered] = (z+i);
      dataCovered++;
    }//endfor
  }//endfor

  CkAssert(dataCovered == numPoints);

//======================================================================
// Find pts with k_x==0 then check the layout

  ihave_kx0 = 0;
  kx0_strt = 0;
  int nkx0  = 0;
  for(i=0;i<numPoints;i++){
    if(k_x[i]==0){
      if(ihave_kx0==0){kx0_strt=i;}
      ihave_kx0=1;
      nkx0++;
    }//endif
  }//endif
  kx0_end = kx0_strt + nkx0;

  for(i=1;i<numPoints;i++){  
    if(k_x[i]<k_x[(i-1)]){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("kx should be stored in increasing order\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

//==============================================================================
// Set the return values

  *n    = numPoints;
  *kk_x = k_x;
  *kk_y = k_y;
  *kk_z = k_z;
  
//==============================================================================
  }//end routine
//==============================================================================



//==============================================================================
// Expand the data set to line form : numLines*nfftz
// Perform the fft along z direction : numLines ffts 
// to get started towards realspace
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
complex* FFTcache::doGSRealFwFFT(complex *packedPlaneData, RunDescriptor *runs, 
                          int numRuns, int numLines,int numFull, int numPoints,
                          int nfftz, bool fftReqd)
//==============================================================================
   {//begin routine
//==============================================================================

    FFTcache *sc   = fftCacheProxy.ckLocalBranch();
    fftw_plan plan = sc->fwdYPlan; // for historic reasons its the `y' plan.
                                   // but sorry, its the z plan.

    expandGSpace(fftData,packedPlaneData,runs,numRuns,numFull,numPoints,nfftz);

    if(fftReqd){
        fftw(plan,        // direction Z now
	     numLines,    // these many ffts : one for every line of z in the chare
	     (fftw_complex *)fftData, //input data
	     1,           //stride
	     nfftz,       //distance between z-data sets
	     NULL, 0, 0); // junk because input array stores the output (in-place)
    }//endif
    return fftData;

//==============================================================================
   }//end routine
//==============================================================================



//==============================================================================
// Perform the back fft along z direction : numLines ffts
// then compress the data to non-zero points as you are back in g-space
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::doBwFFT(complex *fftData) {

   FFTcache *sc   = fftCacheProxy.ckLocalBranch();
   fftw_plan plan = sc->bwdYPlan;  // the Y is a misnomer : its Z
   int expandtype = 2;
   int nfftz      = planeSize[1];

   if(fftReqd){
        fftw(plan,        // direction Z
	     numLines,    // these many ffts : lines of z
	     (fftw_complex *)fftData, //input data
	     1,           //stride
	     nfftz,       //distance between arrays
	     NULL, 0, 0); // junk because input array stores output (in-place)
   }//endfor
   compressGSpace(fftData, expandtype);

   delete [] fftData;
   fftData = NULL;

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the GStatePlane constructor */
//==============================================================================
void initGStateSlab(GStateSlab *gs, int sizeX, size2d size, int gSpaceUnits, 
          int realSpaceUnits, int s_grain, int iplane_ind,int istate_ind) 
//==============================================================================
   {//begin routine
//==============================================================================
// Explanation of organization of data: A point is psi(kx,ky,kz)  
//
//   Each gStateSlab is a collection of pts in lines of constant kx,ky.
//   The maximum number of pts in any line is size[1] = nfftz
//   Some lines are longer/shorter to spherical trunction (|k|<k_cut).
//   Also, kx>=0 when at the Gamma point or doublePack=1, the only
//   thing that has been recently tested. Collections of lines are parallelized. 
//   Each collection does not correspond to a unique plane of kx.
//   The parameter gSpaceUnits is obsolete and must be unity.
//   In order to create the state in real space, psi(x,y,z), an fft 
//   of all kz lines is performed followed by a transpose. The result of the transpose
//   is parallelized by planes of z. This is different than it used to be.
//   The number of points AFTER the fft is nlines*nfftz will be different
//   for each chare array

   if(gSpaceUnits!=1 || realSpaceUnits!=1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("gspacePPC==1 obsolete. We no longer parallelize g-space by plane\n");
     CkPrintf("realSpacePPC==1 although real space is parallelized by plane\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif
   gs->numNonZeroPlanes=1;
   gs->mysizeX   = sizeX;
   gs->planeSize = size;  // fftsizes (sizeY,sizeZ)  sizeX is a global
   gs->numPoints = 0;     // number of packed points (data size before fft)
   gs->numRuns   = 0;     // 2*(number of lines) in the collection
   gs->numLines  = 0;     // (number of lines) in the collection
   gs->fftReqd   = false; // false if this chare has no packed pts
   gs->numFull   = 0;     // number of pts : numLines*nfftz

   gs->S_grainSize = s_grain;   // PC grainsize
   gs->xdim = 1;                // may need some love
   gs->ydim = gs->planeSize[0];
   gs->zdim = gs->planeSize[1];

   gs->iplane_ind   = iplane_ind;
   gs->istate_ind   = istate_ind; 

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RealStateSlab::~RealStateSlab() {
//==============================================================================

  if(planeArr != NULL) {delete [] planeArr;planeArr = NULL; }

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
//FFTcache - sets up a fftcache on each processor
//==============================================================================
FFTcache::FFTcache(size2d planeSIZE, int ArraySize){
//==============================================================================
	
        planeSize = planeSIZE;  //contains nfftz,nfftx

        int nlines_max = scProxy.ckLocalBranch()->cpcharmParaInfo->nlines_max;
        int psize      = nlines_max*planeSize[0];

	fftData = new complex[psize]; 

	//Non Double Pack plans : 2 Y planes and 2 X plans.
	fwdZ1DPlan = fftw_create_plan(planeSize[0], FFTW_FORWARD,  
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	bwdZ1DPlan = fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	fwdX1DPlan = fftw_create_plan(sizeX, FFTW_FORWARD, 
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	bwdX1DPlan = fftw_create_plan(sizeX, FFTW_BACKWARD, 
                     FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	int size[3];
	size[0] = planeSize[1]; size[1] = planeSize[0]; size[2] = 1;
	//Double Pack Plans
        fwdZ1DdpPlan = fftw_create_plan(planeSize[0], FFTW_FORWARD, 
                          FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	
	bwdZ1DdpPlan = fftw_create_plan(planeSize[0], FFTW_BACKWARD, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
	fwdX1DdpPlan = rfftwnd_create_plan(1, (const int*)size, FFTW_COMPLEX_TO_REAL, 
                       FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
        bwdX1DdpPlan = rfftwnd_create_plan(1, (const int*)size,FFTW_REAL_TO_COMPLEX,
                        FFTW_MEASURE | FFTW_IN_PLACE|FFTW_USE_WISDOM);
        int sizeZ = planeSize[0];
	fwdYPlan = fftw_create_plan(sizeZ,FFTW_FORWARD, 
                                 FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);
	bwdYPlan = fftw_create_plan(sizeZ,FFTW_BACKWARD, 
                                 FFTW_MEASURE|FFTW_IN_PLACE|FFTW_USE_WISDOM);

//==============================================================================
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

  if(config.doublePack && config.inPlaceFFT) {
      data = new double[pSize];
      int stride = sizeX/2+1;
      // FFT along Y direction : Y moves with stride sizex/2+1 through memory
      fftw(fwdZ1DdpPlan,    // y-plan (label lies)
  	   nplane_x,        // how many < sizeX/2 + 1
	   (fftw_complex *)(planeArr),//input data
 	   stride,          // stride betwen elements (y is inner)
  	   1,               // array separation (nffty elements)
	   NULL,0,0);       // output data is input data
      // fftw only gives you one sign for real to complex : so do it yourself
      for(int i=0;i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}
      rfftwnd_complex_to_real(fwdX1DdpPlan,
			      planeSize[0],    // how many
			      (fftw_complex *)planeArr, 
                              1,               // stride (x is inner)
                              stride,          // array separation
    		              NULL,0,0);       // output array is real 
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

 if(!config.doublePack && config.inPlaceFFT){
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
void RealStateSlab::zeroOutPlanes() {    
  memset(planeArr,0,size*sizeof(complex));
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::allocate() {
  if(planeArr == NULL) {planeArr = new complex [size];}
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::destroy() {
  if(planeArr != NULL) {delete [] planeArr; planeArr=NULL;}
}
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

#ifdef _CP_DEBUG_VKS_RSPACE_
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

  if(config.doublePack && config.inPlaceFFT){
     double *realArr = reinterpret_cast<double*> (planeArr);
     int stride = sizeX/2+1;
     for(int i=0,i2=0;i<planeSize[0];i++,i2+=2){
       for(int j=i*sizeX;j<(i+1)*sizeX;j++){
         realArr[(j+i2)] = realArr[(j+i2)]*vks[j];
       }//endfor
     }//endfor
     rfftwnd_real_to_complex(bwdX1DdpPlan,
			planeSize[0],          // these many 1D ffts
			realArr,1,(sizeX+2),   // x is inner here
                        NULL,0,0);            
     for (int i=0; i<stride*planeSize[0];i++){planeArr[i].im = -planeArr[i].im;}
     fftw(bwdZ1DdpPlan,
	     nplane_x, // these many 1D ffts
	     (fftw_complex *)planeArr, 
  	     stride,1,NULL,0,0);
  }//endif : double pack and in-place

//==============================================================================

  if(!config.doublePack && config.inPlaceFFT){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Fix the non-double pack FFT\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the CP_State_RealSpacePlane constructor */
//==============================================================================
void initRealStateSlab(RealStateSlab *rs, size2d planeSize, int gSpaceUnits, 
                       int realSpaceUnits, int stateIndex, int planeIndex)
//==============================================================================
   {//begin routine
//==============================================================================
// Explanation of the organization:
//   Each CP_State_RealSpacePlane instance is actually a planes of psi_I(x,y,z).
//   We index into planeArr using the z-coordinate
//   and each chare is an x-y plane.
//   Data arrives from the gspaceplane of total size : nlines_tot*nfftz

   if(gSpaceUnits!=1 || realSpaceUnits!=1){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("gspacePPC==1 obsolete. We no longer parallelize g-space by plane\n");
     CkPrintf("realSpacePPC==1 although real space is parallelized by plane\n");
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
   }//endif

   rs->planeSize  = planeSize;         // FFTSize Y and Z, sizeX is global
   rs->thisState  = stateIndex;        // my state (I)
   rs->thisPlane  = planeIndex;        // my plane (z)
   rs->nsize      = planeSize[0]*sizeX; // when fft is completed
   rs->numPlanesToExpect = scProxy.ckLocalBranch()->cpcharmParaInfo->nplane_x;
   int rsize_a    = planeSize[0]*(planeSize[1]/2+1);
   int rsize_b    = planeSize[1]*(planeSize[0]/2+1);
   rs->rsize      = (rsize_a > rsize_b ? rsize_a : rsize_b);

   if(config.doublePack)  {
      rs->size = rs->rsize;
   }else{
      rs->size = rs->nsize;
   }//endif

   rs->planeArr = new complex[(rs->size)];

//==============================================================================    
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoRealSlab::~RhoRealSlab()
{
	fftwnd_destroy_plan(fft2dFwPlan);
	fftwnd_destroy_plan(fft2dBwPlan);
	//delete [] Vks;
	delete density;
	delete [] doFFTonThis;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
complex *RhoRealSlab::doBwFFT()
{

        memcpy(density, doFFTonThis,  sizeX * sizeY * sizeof(complex));
     
	int planeSize       = sizeX * sizeZ;
        complex *fftResults = new complex[planeSize];

        // do x-y plane fft in-place
        fftwnd_one(fft2dBwPlan, (fftw_complex *)fftResults, NULL);
        memcpy(fftResults,doFFTonThis, planeSize * sizeof(complex));

return doFFTonThis;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoRealSlab::doFwFFT()
{
	int j, planeSize = sizeX * sizeZ;
#ifndef CMK_OPTIMIZE
	double StartTime=CmiWallTimer();
#endif

        fftwnd_one(fft2dFwPlan,(fftw_complex *)doFFTonThis, NULL);
#ifndef CMK_OPTIMIZE
	traceUserBracketEvent(VksofRFFT_, StartTime, CmiWallTimer());    
#endif

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the RealSpaceDensity Constructor */
//==============================================================================
void initRhoRealSlab(RhoRealSlab *rho_rs, int xdim, int ydim, int zdim,
                     int numRealSpace, int numRhoG, int myIndex)
//==============================================================================
   {//begin routine
//==============================================================================

	rho_rs->sizeX = xdim;
	rho_rs->sizeY = ydim;
	rho_rs->sizeZ = zdim;

	rho_rs->fft2dFwPlan = fftw2d_create_plan(rho_rs->sizeZ, rho_rs->sizeX, 
                FFTW_FORWARD, FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
	rho_rs->fft2dBwPlan = fftw2d_create_plan(rho_rs->sizeZ, rho_rs->sizeX, 
                FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);

	rho_rs->size = rho_rs->sizeX * rho_rs->sizeZ;
	rho_rs->xdim = rho_rs->sizeX;
	rho_rs->Vks     = new complex[rho_rs->size]; 
	rho_rs->density = new complex[rho_rs->size];

	//we copy the real densities into the real parts of these and do our stuff.
	rho_rs->doFFTonThis = new complex[rho_rs->size];
	rho_rs->ydim = 1;
	rho_rs->zdim = rho_rs->sizeZ;
	rho_rs->startx = 0;
	rho_rs->starty = 0;
	rho_rs->startz = 0; 
	// rho_rs->type = XY_PLANE; 

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RhoGSlab::~RhoGSlab()
{
	fftw_destroy_plan(fft1dFwPlan);
	fftw_destroy_plan(fft1dBwPlan);
	delete [] chunk;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
// complete the transform of rho(r) to rho(g) : last FFT after tranpose
// 
//==============================================================================
void RhoGSlab::doBwFFT(int index){

    // do the line inv-fft
    complex *temp1 = new complex[sizeY];
    complex *temp2 = new complex[sizeY];

    int x, y;
    for (x = 0; x < sizeX; x++) {
      for (y = 0; y < sizeY; y++){temp1[y] = chunk[y * sizeX + x];}
      fftw_complex *temp_in  = reinterpret_cast<fftw_complex *>(temp1);
      fftw_complex *temp_out = reinterpret_cast<fftw_complex *>(temp2);
      fftw_one(fft1dBwPlan,temp_in,temp_out);
      for (y = 0; y < sizeY; y++){chunk[y * sizeX + x] = temp1[y];}
    }//endfor

    delete [] temp1;
    delete [] temp2;
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RhoGSlab::doFwFFT() {

    complex *temp1 = new complex[sizeY];
    complex *temp2 = new complex[sizeY];

    // do the fwd-fft
    int y, x;
    for (x = 0; x < sizeX; x++) {
      for (y = 0; y < sizeY; y++){temp1[y] = chunk[y * sizeX + x];}
      fftw_complex *temp_in  = reinterpret_cast<fftw_complex *>(temp1);
      fftw_complex *temp_out = reinterpret_cast<fftw_complex *>(temp2);
      fftw_one(fft1dFwPlan,temp_in,temp_out);
      for (y = 0; y < sizeY; y++){chunk[y * sizeX + x] = temp1[y];}
    }//endfor x

    delete [] temp1;
    delete [] temp2;

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the CP_Rho_GSpacePlane constructor  */
//==============================================================================
void initRhoGSlab(RhoGSlab *rho_gs, int xdim, int ydim, int zdim, int numRealSpace, 
                            int numRealSpaceDensity, int myIndex)
//==============================================================================
   {//begin routine
//==============================================================================

	rho_gs->sizeX = xdim;
	rho_gs->sizeY = ydim;
	rho_gs->sizeZ = zdim;
	rho_gs->fft1dFwPlan = fftw_create_plan(
          rho_gs->sizeY,FFTW_FORWARD,FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);
	rho_gs->fft1dBwPlan = fftw_create_plan(
          rho_gs->sizeY,FFTW_BACKWARD,FFTW_MEASURE | FFTW_IN_PLACE | FFTW_USE_WISDOM);

	// NOTE on in-memory  handling of the chunk:
	// The chunk is handled as a series of x-z planes, with pencils in x direction
	rho_gs->chunk = new complex[rho_gs->sizeY * rho_gs->sizeX];

	rho_gs-> size = rho_gs->sizeX * rho_gs->sizeY;
	/* complex *rho_g = rho_gs.chunk;*/
	rho_gs->xdim = rho_gs->sizeX;
	rho_gs->ydim = rho_gs->sizeY;
	rho_gs->zdim = 1;
	rho_gs->startx = 0;
	rho_gs->starty = 0;
	rho_gs->startz = myIndex;
	//   rho_gs->type = ZX_PLANE;

	return;

//==============================================================================
   }//end routine
//==============================================================================
