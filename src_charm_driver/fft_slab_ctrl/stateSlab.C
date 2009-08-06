//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file stateSlab.C
 * Add functions to allow application programmers to initialize these and
 * the corresponding functions in Charm++ to call these functions with
 * appropriate parameters 
 */
//==============================================================================

#include "utility/util.h"
#include "debug_flags.h"
#include <cmath>
#include "main/cpaimd.h"
#include "main/groups.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_ctrl/CP_State_Plane.h"
#include "main/eesCache.h"

//==============================================================================

extern CkVec <CProxy_FFTcache>           UfftCacheProxy;
extern Config                    config;
extern CkVec <CProxy_AtomsGrp>           UatomsGrpProxy;
extern CProxy_CPcharmParaInfoGrp scProxy;
extern CkVec <CProxy_eesCache>           UeesCacheProxy;
extern int nstates;
extern int sizeX;

//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the GStatePlane constructor */
//==============================================================================
void initGStateSlab(GStateSlab *gs, int sizeX, int sizeY, int sizeZ, int gSpaceUnits, 
          int realSpaceUnits, int s_grain, int iplane_ind,int istate_ind,
          int len_nhc_cp, int num_nhc_cp, int nck_nhc_cp) 
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
   gs->planeSize = size2d(sizeY,sizeZ);  // fftsizes (sizeY,sizeZ)  sizeX is a global
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
   gs->initNHC(len_nhc_cp,num_nhc_cp,nck_nhc_cp);

//==============================================================================
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
GStateSlab::~GStateSlab() {

    if(packedPlaneData    !=NULL) fftw_free( packedPlaneData);
    if(packedForceData    !=NULL) fftw_free( packedForceData);
    if(packedVelData      !=NULL) fftw_free( packedVelData);
    if(cp_min_opt==0){
      if(packedPlaneDataScr !=NULL) fftw_free( packedPlaneDataScr);
    }//endif
#ifdef  _CP_DEBUG_UPDATE_OFF_
    if(cp_min_opt==1){
      if(packedPlaneDataTemp!=NULL) fftw_free( packedPlaneDataTemp);
    }//endif
#endif
    destroyNHC();

    packedPlaneData     = NULL;
    packedForceData     = NULL;
    packedVelData       = NULL;
    packedPlaneDataScr  = NULL;
    packedPlaneDataTemp = NULL;

}
//==============================================================================

//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::pup(PUP::er &p) {
//==============================================================================
// Dont have to pup fftw plans - they live in the fft cache group

//  CkPrintf("gs pup\n:");
        p|cp_min_opt;
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
        p|potNHC_ret;
        p|degfree;
        p|degfreeNHC;

	if (p.isUnpacking()) {
           packedPlaneData     = (complex *)fftw_malloc(numPoints*sizeof(complex));
	   packedForceData     = (complex *)fftw_malloc(numFull*sizeof(complex));
	   packedVelData       = (complex *)fftw_malloc(numPoints*sizeof(complex));
           packedRedPsi        = (complex *)fftw_malloc(nkx0*sizeof(complex));
           if(cp_min_opt==0){
  	     packedPlaneDataScr  = (complex *)fftw_malloc(numPoints*sizeof(complex));
	   }//endif
#ifdef  _CP_DEBUG_UPDATE_OFF_
           if(cp_min_opt==1){
  	     packedPlaneDataTemp = (complex *)fftw_malloc(numPoints*sizeof(complex));
	   }//endif
#endif
	}//endif
        p((char *) packedPlaneData, numPoints*sizeof(complex));
	p((char *) packedForceData, numFull*sizeof(complex));
	p((char *) packedVelData, numPoints*sizeof(complex));   //cg under min
	p((char *) packedRedPsi, nkx0*sizeof(complex));
        if(cp_min_opt==0){
  	  p((char *) packedPlaneDataScr, numPoints*sizeof(complex));
	}//endif
#ifdef  _CP_DEBUG_UPDATE_OFF_
        if(cp_min_opt==1){
  	  p((char *) packedPlaneDataTemp, numPoints*sizeof(complex));
	}//endif
#endif

        p|nck_nhc_cp;
        p|len_nhc_cp;
        p|num_nhc_cp;
        p|kTCP;
        p|tauNHCCP;
	if (p.isUnpacking()) {
          initNHC(len_nhc_cp,num_nhc_cp,nck_nhc_cp);
	}//endif recving
        int nsize   = num_nhc_cp*len_nhc_cp*nck_nhc_cp;
        double *xt  = new double[nsize];
        double *xtp = new double[nsize];
        double *vt  = new double[nsize];
        double *ft  = new double[nsize];
	if(p.isPacking()){
          int iii=0;
          for(int k =0;k<num_nhc_cp;k++){
          for(int i =0;i<num_nhc_cp;i++){
          for(int j =0;j<len_nhc_cp;j++){
            xt[iii]  = xNHC[k][i][j];
            xtp[iii] = xNHCP[k][i][j];
            vt[iii]  = vNHC[k][i][j];
            ft[iii]  = fNHC[k][i][j];
            iii++;
          }}}
	}//endif sending
        p(xt,nsize);
        p(xtp,nsize);
        p(vt,nsize);
        p(ft,nsize);
        p(degFreeSplt,nck_nhc_cp);
        p(istrNHC,nck_nhc_cp);
        p(iendNHC,nck_nhc_cp);
        p(mNHC,len_nhc_cp);
        p(v0NHC,num_nhc_cp);
        p(a2NHC,num_nhc_cp);
        p(a4NHC,num_nhc_cp);
	if (p.isUnpacking()) {
          int iii=0;
          for(int k =0;k<nck_nhc_cp;k++){
          for(int i =0;i<num_nhc_cp;i++){
          for(int j =0;j<len_nhc_cp;j++){
            xNHC[k][i][j]  = xt[iii];
            xNHCP[k][i][j] = xtp[iii];
            vNHC[k][i][j]  = vt[iii];
            fNHC[k][i][j]  = ft[iii];
            iii++;
          }}}
	}//endif receiving 
        delete []xt;
        delete []xtp;
        delete []vt;
        delete []ft;

//==============================================================================
  }//end routine
//==============================================================================




//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GStateSlab::addForces(complex *points, const int *k_x){
//==============================================================================
    int i;
    double wght;

    if(!config.doublePack){

      for(i = 0; i < numPoints; i++){
	 packedForceData[i] += points[i];
      }//endfor

    }else{

      int nfreq=1000;
      for(i = 0; i < nkx0; i++){
#ifdef _CP_DEBUG_VKS_OFF_  // only non-local, no vks forces
	packedForceData[i].re = 0.0; packedForceData[i].im = 0.0; 
#endif
	packedForceData[i] += points[i];
#ifdef CMK_BLUEGENEL
        if(i%nfreq==0){CmiNetworkProgress();}
#endif
      }//endfor

      for(i = nkx0; i < numPoints; i++){
#ifdef _CP_DEBUG_VKS_OFF_       // only non-local, no vks forces
	packedForceData[i].re = 0.0; packedForceData[i].im = 0.0; 
#endif
	packedForceData[i] *= 2.0; 
	packedForceData[i] += points[i];
#ifdef CMK_BLUEGENEL
        if(i%nfreq==0){CmiNetworkProgress();}
#endif
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

void GStateSlab::setKRange(int n, int *k_x, int *k_y, int *k_z){

//======================================================================
// Construct the k-vectors

  CkAssert(n == numPoints);

//======================================================================
// Find pts with k_x==0 then check the layout : kx=0 first

  int i;
  ihave_g000 = 0;
  ind_g000   = -1;
  ihave_kx0  = 0;
  nkx0       = 0;
  nkx0_uni   = 0;
  nkx0_red   = 0;
  nkx0_zero  = 0;
  kx0_strt   = 0;
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
    CkAbort("kx=0 should be stored first | kx_srt !=0\n");
  }//endif
  
  if(nkx0!=nkx0_uni+nkx0_red){
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkPrintf("Incorrect count of redundant guys\n");
    CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
    CkAbort("Incorrect count of redundant guys\n");
  }//endif

  for(i=0;i<nkx0;i++){  
    if(k_x[i]!=0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("kx should be stored consecutively and first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkAbort("kx should be stored consecutively and first\n");
    }//endif
  }//endif

  for(i=0;i<nkx0_red;i++){  
    if(k_y[i]>0 || (k_y[i]==0 && k_z[i]>=0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkAbort("ky <0 should be stored first\n");
    }//endif
  }//endfor

  for(i=nkx0_red;i<nkx0_uni;i++){  
    if(k_y[i]<0 || (k_y[i]==0 && k_z[i]<0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkAbort("ky <0 should be stored first\n");
    }//endif
  }//endfor

  packedRedPsi   = (complex *)fftw_malloc(nkx0*sizeof(complex));
  packedRedPsiV  = (complex *)fftw_malloc(nkx0*sizeof(complex));
  memset(packedRedPsi,  0, sizeof(complex)*nkx0);
  memset(packedRedPsiV, 0, sizeof(complex)*nkx0);

//==============================================================================
  }//end routine
//==============================================================================



//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/* This gets called at the end of the CP_State_RealSpacePlane constructor */
//==============================================================================
void initRealStateSlab(RealStateSlab *rs, int ngrid_a,int ngrid_b, int ngrid_c,
                       int gSpaceUnits, int realSpaceUnits, int stateIndex, int planeIndex)
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

   rs->thisState  = stateIndex;                // my state (I)
   rs->thisPlane  = planeIndex;                // my plane (z)
   rs->nsize      = ngrid_a*ngrid_b;             // when fft is completed
   rs->numPlanesToExpect = scProxy.ckLocalBranch()->cpcharmParaInfo->nchareG;
   int rsize_a    = ngrid_a*(ngrid_b/2+1);
   int rsize_b    = ngrid_b*(ngrid_a/2+1);
   rs->rsize      = (rsize_a > rsize_b ? rsize_a : rsize_b);

   if(config.doublePack)  {
      rs->size = rs->rsize;
   }else{
      rs->size = rs->nsize;
   }//endif

   rs->planeArr  = (complex *) fftw_malloc(rs->size * sizeof(complex));
   rs->planeArrR = reinterpret_cast<double*> (rs->planeArr);

   rs->ngrid_a = ngrid_a;
   rs->ngrid_b = ngrid_b;

//==============================================================================    
   }//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
RealStateSlab::~RealStateSlab() {
//==============================================================================

  if(planeArr != NULL) {fftw_free(planeArr);planeArr = NULL; planeArrR=NULL;}

}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::pup(PUP::er &p) {
//==============================================================================
  p|ngrid_a;
  p|ngrid_b;
  p|size;
  if (p.isUnpacking()) {
      planeArr = (complex *) fftw_malloc(size*sizeof(complex));
      planeArrR = reinterpret_cast<double*> (planeArr);
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
void RealStateSlab::allocate() {
  if(planeArr == NULL) {
     planeArr = (complex *) fftw_malloc(size*sizeof(complex));
     planeArrR = reinterpret_cast<double*> (planeArr);
  }//endif
}
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RealStateSlab::destroy() {
  if(planeArr != NULL) {fftw_free(planeArr); planeArr=NULL; planeArrR = NULL;}
}
//==============================================================================


