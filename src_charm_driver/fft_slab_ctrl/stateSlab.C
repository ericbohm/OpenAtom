//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file stateSlab.C
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
   gs->xNHC=0.0;   //to satisfy valgrind
   gs->iplane_ind   = iplane_ind;
   gs->istate_ind   = istate_ind; 
   gs->initNHC();

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
GStateSlab::~GStateSlab() {

    if(packedPlaneData    !=NULL) delete [] packedPlaneData;
    if(packedPlaneDataTemp!=NULL) delete [] packedPlaneDataTemp;
    if(packedForceData    !=NULL) delete [] packedForceData;
    if(packedPlaneDataScr !=NULL) delete [] packedPlaneDataScr;
    if(packedVelData      !=NULL) delete [] packedVelData;
    if(runs               !=NULL) delete [] runs;
    destroyNHC();

    packedPlaneData     = NULL;
    packedPlaneDataTemp = NULL;
    packedForceData     = NULL;
    packedPlaneDataScr  = NULL;
    packedVelData       = NULL;
    runs                = NULL;

}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::pup(PUP::er &p) {
//==============================================================================
// Dont have to pup fftw plans - they live in the fft cache group

//  CkPrintf("gs pup\n:");
        p|ees_nonlocal;
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
        p|ngridaNL;
        p|ngridbNL;
        p|ngridcNL;
	p|iplane_ind;
        p|istate_ind;
        p|ihave_kx0;
        p|ihave_g000;
        p|ind_g000;
        p|kx0_strt;
        p|kx0_end;
        p|nkx0; p|nkx0_uni; p|nkx0_red; p|nkx0_zero;
        p|eke_ret;
        p|fictEke_ret;
        p|ekeNhc_ret;
        p|degfree;
        p|degfreeNHC;
        p|gammaNHC;

	if (p.isUnpacking()) {runs = new RunDescriptor[numRuns];}
	for (int i = 0; i < numRuns; i++){runs[i].pup(p);}

	if (p.isUnpacking()) {
           packedPlaneData     = new complex[numPoints];
	   packedPlaneDataScr  = new complex[numPoints];
	   packedPlaneDataTemp = new complex[numPoints];
	   packedForceData     = new complex[numPoints];
	   packedVelData       = new complex[numPoints];
           packedRedPsi        = new complex[nkx0];
	}//endif
	p((char *) packedPlaneData, numPoints*sizeof(complex));
	p((char *) packedPlaneDataScr, numPoints*sizeof(complex));
	p((char *) packedPlaneDataTemp, numPoints*sizeof(complex));
	p((char *) packedForceData, numPoints*sizeof(complex));
	p((char *) packedVelData, numPoints*sizeof(complex));
	p((char *) packedRedPsi, nkx0*sizeof(complex));

        p|len_nhc_cp;
        p|num_nhc_cp;
        p|kTCP;
        p|tauNHCCP;
        p|xNHC;
        p|mNHC;
        int nsize  = num_nhc_cp*len_nhc_cp;
        double *vt = new double[nsize];
        double *ft = new double[nsize];
	if (!p.isUnpacking()) {
          int iii=0;
          for(int i =0;i<num_nhc_cp;i++){
          for(int j =0;j<len_nhc_cp;j++){
            vt[iii] = vNHC[i][j];
            ft[iii] = fNHC[i][j];
            iii++;
          }}
	}//endif sending
        p(vt,nsize);
        p(ft,nsize);
	if (p.isUnpacking()) {
          initNHC();
          int iii=0;
          for(int i =0;i<num_nhc_cp;i++){
          for(int j =0;j<len_nhc_cp;j++){
            vNHC[i][j] = vt[iii];
            fNHC[i][j] = ft[iii];
            iii++;
          }}
	}//endif receiving 
        delete []vt;
        delete []ft;
	//  CkPrintf("end gs pup\n:");

//==============================================================================
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
// Find pts with k_x==0 then check the layout : kx=0 first

  ihave_g000 = 0;
  ind_g000   = -1;
  ihave_kx0 = 0;
  nkx0  = 0;
  nkx0_uni=0;
  nkx0_red=0;
  nkx0_zero=0;
  kx0_strt =0;
  for(i=0;i<numPoints;i++){
    if(k_x[i]==0 && k_y[i]>0){nkx0_uni++;}
    if(k_x[i]==0 && k_y[i]<0){nkx0_red++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]>=0){nkx0_uni++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]<0){nkx0_red++;}
    if(k_x[i]==0 && k_y[i]==0 && k_z[i]==0){nkx0_zero++;ihave_g000=1;ind_g000=i;}
    if(k_x[i]==0){
      if(ihave_kx0==0){kx0_strt=i;}
      ihave_kx0=1;
      nkx0++;
    }//endif
  }//endif
  kx0_end = kx0_strt + nkx0;

  if(kx0_strt!=0){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("kx=0 should be stored first | kx_srt !=0\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif
  
  if(nkx0!=nkx0_uni+nkx0_red){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Incorrect count of redundant guys\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkExit();
  }//endif

  for(i=0;i<nkx0;i++){  
    if(k_x[i]!=0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("kx should be stored consecutively and first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endif

  for(i=0;i<nkx0_red;i++){  
    if(k_y[i]>0 || (k_y[i]==0 && k_z[i]>=0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endfor

  for(i=nkx0_red;i<nkx0_uni;i++){  
    if(k_y[i]<0 || (k_y[i]==0 && k_z[i]<0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
  }//endfor

  packedRedPsi  = new complex[nkx0];
  memset(packedRedPsi, 0, sizeof(complex)*nkx0);

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
        fft_split(plan,     // direction Z
	     numLines,      // these many ffts : lines of z
	     (fftw_complex *)fftData, //input data
	     1,             //stride
	     nfftz,         //distance between arrays
	     NULL, 0, 0,    // junk because input array stores output (in-place)
             config.fftprogresssplit);
   }//endfor
   compressGSpace(fftData, expandtype);

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
   rs->numPlanesToExpect = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
   int rsize_a    = planeSize[0]*(planeSize[1]/2+1);
   int rsize_b    = planeSize[1]*(planeSize[0]/2+1);
   rs->rsize      = (rsize_a > rsize_b ? rsize_a : rsize_b);

   if(config.doublePack)  {
      rs->size = rs->rsize;
   }else{
      rs->size = rs->nsize;
   }//endif
   rs->planeArr = (complex *) fftw_malloc(rs->size * sizeof(complex));

//==============================================================================    
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RealStateSlab::~RealStateSlab() {
//==============================================================================

  if(planeArr != NULL) {fftw_free(planeArr);planeArr = NULL; }

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::pup(PUP::er &p) {
//==============================================================================
  p|planeSize;
  p|size;
  if (p.isUnpacking()) 
    {
      planeArr = (complex *) fftw_malloc(size*sizeof(complex));
    }
  PUParray(p, planeArr, size);
  p|thisState;    
  p|thisPlane;      
  p|numPlanesToExpect;
  p|nsize;
  p|rsize;
  p|size;
  p|e_gga; 
}
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
  if(planeArr == NULL) {planeArr = (complex *) fftw_malloc(size*sizeof(complex));}
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::destroy() {
  if(planeArr != NULL) {fftw_free(planeArr); planeArr=NULL;}
}
//==============================================================================


