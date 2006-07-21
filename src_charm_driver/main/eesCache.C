//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file eesCache.C
 * 
 *
 */
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================


//==============================================================================
#include "charm++.h"
#include "util.h"
#include "cpaimd.h"
#include "groups.h"
#include "eesCache.h"
#include <math.h>
#include "fftCacheSlab.h"
#include "CP_State_Plane.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "../../src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"


//----------------------------------------------------------------------------
extern CProxy_CP_State_ParticlePlane     particlePlaneProxy;
extern CProxy_CP_State_RealParticlePlane realParticlePlaneProxy;
extern CProxy_CP_Rho_RHartExt            rhoRHartExtProxy;
extern CProxy_CP_Rho_GHartExt            rhoGHartExtProxy;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CProxy_AtomsGrp                   atomsGrpProxy;
extern CProxy_eesCache                   eesCacheProxy;

#define _EESCACHE_VERBOSE_OFF_

//==============================================================================



//==============================================================================
// Group constructor
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
eesCache::eesCache(int _nchareRPP, int _nchareGPP, int _nchareRHart,
                   int _nchareGHart, int _nstates, int _nchareRhoG)
//==============================================================================
   {//begin rotuine
//==============================================================================

   itimeRPP        = 0;
   itimeRHart      = 0;
   nMallSize       = 100;

   nchareRPP       = _nchareRPP;
   nchareGPP       = _nchareGPP;
   nchareRHart     = _nchareRHart;
   nchareGHart     = _nchareGHart;
   nchareGSP       = _nchareGPP;
   nchareRhoG      = _nchareRhoG; 
   nstates         = _nstates;

   nchareGSPProcT  = 0;
   nchareGSPProc   = 0;
   nchareRPPProc   = 0;
   nchareGPPProc   = 0;
   nchareRHartProc = 0;
   nchareGHartProc = 0;

   gspStateInd     = new int [nMallSize];
   gspPlaneInd     = new int [nMallSize];

   allowedRppChares      = new int[nchareRPP];
   allowedGppChares      = new int[nchareGPP];
   allowedRhoRHartChares = new int[nchareRHart];
   allowedRhoGHartChares = new int[nchareGHart];
   allowedGspChares      = new int[nchareGSP];
   allowedRhoGChares     = new int[nchareRhoG];

   for(int i=0;i<nchareRPP  ;i++){allowedRppChares[i]      = 0;}
   for(int i=0;i<nchareGPP  ;i++){allowedGppChares[i]      = 0;}
   for(int i=0;i<nchareRHart;i++){allowedRhoRHartChares[i] = 0;}
   for(int i=0;i<nchareGHart;i++){allowedRhoGHartChares[i] = 0;}
   for(int i=0;i<nchareGSP;  i++){allowedGspChares[i]      = 0;}
   for(int i=0;i<nchareRhoG ;i++){allowedRhoGChares[i]     = 0;}

   indGppChares      = new int[nchareGPP];  // over dimensioned
   indRppChares      = new int[nchareRPP];
   indRhoRHartChares = new int[nchareRHart];
   indRhoGHartChares = new int[nchareGHart];
   indGspChares      = new int[nchareGSP];
   indRhoGChares     = new int[nchareRhoG];

   GppData           = new GPPDATA     [nchareGPP]; // over dimensioned
   RppData           = new RPPDATA     [nchareRPP];
   RhoGHartData      = new RHOGHARTDATA[nchareGHart];
   RhoRHartData      = new RHORHARTDATA[nchareRHart];
   GspData           = new GSPDATA[nchareGSP];
   RhoGData          = new RHOGDATA[nchareRhoG];

}// end constructor
//==============================================================================


//==============================================================================
// realParticlePlane Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::registerCacheRPP  (int index){
//==============================================================================
// Cache warm up

  if(allowedRppChares[index]==0){
    nchareRPPProc          += 1;
    allowedRppChares[index] = 1;
    RppData[index].init(index);
  }//endif

}//end routine
//==============================================================================


//==============================================================================
// GParticlePlane Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::registerCacheGPP  (int index, int ncoef, int *ka, int *kb, int *kc){
//==============================================================================
// Cache warm up

   if(allowedGppChares[index]==0){
     nchareGPPProc          += 1;
     allowedGppChares[index] = 1;
     GppData[index].init(index,ncoef,ka,kb,kc);
   }//endif

}//end routine
//==============================================================================


//==============================================================================
// RhoRhart Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::registerCacheRHart(int index){
//==============================================================================
// Cache warm up

  if(allowedRhoRHartChares[index]==0){
#ifdef  _EESCACHE_VERBOSE_
    CkPrintf("Registering Rhart %d\n",index);
#endif
    nchareRHartProc             += 1;
    allowedRhoRHartChares[index] = 1;
    RhoRHartData[index].init(index);
  }//endif

}//end routine
//==============================================================================


//==============================================================================
// RhoGhart Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::registerCacheGHart(int index, int ncoef, int *ka, int *kb, int *kc){
//==============================================================================
// Cache warm up

  if(allowedRhoGHartChares[index]==0){
    nchareGHartProc             += 1;
    allowedRhoGHartChares[index] = 1;
    RhoGHartData[index].init(index,ncoef,ka,kb,kc);
  }//endif

}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::registerCacheGSP(int is ,int ip){

  if(allowedGspChares[ip]==0){
    nchareGSPProc       += 1;
    allowedGspChares[ip] = 1;
    GspData[ip].init(ip);
  }//endif

  int i = nchareGSPProcT;
  gspStateInd[i]  = is;
  gspPlaneInd[i]  = ip;
  nchareGSPProcT += 1;

  if( (nchareGSPProcT % nMallSize)==0){
    int *tempS = new int [(nchareGSPProcT+nMallSize)];
    int *tempP = new int [(nchareGSPProcT+nMallSize)];
    for(int j=0;j<nchareGSPProcT;j++){
      tempS[j] = gspStateInd[j];
      tempP[j] = gspPlaneInd[j];
    }//endfor
    delete [] gspStateInd;
    delete [] gspPlaneInd;
    gspStateInd = tempS;
    gspPlaneInd = tempP;
  }//endif
}
//==============================================================================



//==============================================================================
// Initialize the RealParticlePlane Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RPPDATA::init(int index_in){

   CPNONLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);

   index          = index_in;
   int n_interp21 = n_interp*n_interp+1; // extra space for piny

   plane_index    = (int *)fftw_malloc(natm*sizeof(int));

   igrid = (int **)fftw_malloc(natm*sizeof(int*));
   mn    = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_x = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_y = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_z = (double **)fftw_malloc(natm*sizeof(double*));
   for(int i=0;i<natm;i++){
     igrid[i] = (int *)fftw_malloc(n_interp21*sizeof(int));
     mn[i]    = (double *)fftw_malloc(n_interp21*sizeof(double));
     dmn_x[i] = (double *)fftw_malloc(n_interp21*sizeof(double));
     dmn_y[i] = (double *)fftw_malloc(n_interp21*sizeof(double));
     dmn_z[i] = (double *)fftw_malloc(n_interp21*sizeof(double));
   }//endfor

}//end routine
//==============================================================================


//==============================================================================
// Initialize the GParticlePlane Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GPPDATA::init(int index_in,int ncoef_in, int *ka, int *kb, int *kc){

  index    = index_in;
  ncoef    = ncoef_in;
  b_re     = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  b_im     = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  h_gspl   = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  ind_gspl = (int *)fftw_malloc((ncoef+1)*sizeof(int));

  CPNONLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);
  CPNONLOCAL::eesSetEesWghtGgrp(ncoef,ka,kb,kc,b_re,b_im,ngrid_a,ngrid_b,ngrid_c,
                                n_interp);
  CPNONLOCAL::eesSplProjectorGgrp(ncoef,ka,kb,kc,h_gspl,ind_gspl);

}
//==============================================================================


//==============================================================================
// Initialize the RhoRhart Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RHORHARTDATA::init(int index_in){

   CPLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);

   index          = index_in;
   int n_interp21 = n_interp*n_interp+1;

   plane_index    = (int *)fftw_malloc(natm*sizeof(int));

   igrid = (int **)fftw_malloc(natm*sizeof(int*));
   mn    = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_x = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_y = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_z = (double **)fftw_malloc(natm*sizeof(double*));
   for(int i=0;i<natm;i++){
     igrid[i] = (int *)fftw_malloc(n_interp21*sizeof(int));
     mn[i]    = (double *)fftw_malloc(n_interp21*sizeof(double));
     dmn_x[i] = (double *)fftw_malloc(n_interp21*sizeof(double));
     dmn_y[i] = (double *)fftw_malloc(n_interp21*sizeof(double));
     dmn_z[i] = (double *)fftw_malloc(n_interp21*sizeof(double));
   }//endfor

}//end routine
//==============================================================================


//==============================================================================
// Initialize the RhoGhart Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RHOGHARTDATA::init(int index_in,int ncoef_in, int *ka, int *kb, int *kc){

  index    = index_in;
  ncoef    = ncoef_in;
  b_re     = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  b_im     = (double *)fftw_malloc((ncoef+1)*sizeof(double));

  CPLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);
  CPLOCAL::eesSetEesWghtGgrp(ncoef,ka,kb,kc,b_re,b_im,ngrid_a,ngrid_b,ngrid_c,
                             n_interp);
}
//==============================================================================


//==============================================================================
// realParticlePlane Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::queryCacheRPP  (int index,int itime,int iter){
//==============================================================================
// Cache compute

  if(itime != itimeRPP && iter==1){

    if(itime!=itimeRPP+1 || iter != 1 || allowedRppChares[index]!=1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Broken RPP cache query by %d at %d %d %d\n",index,itime,itimeRPP,iter);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    itimeRPP= itime;

#ifdef  _EESCACHE_VERBOSE_
    CkPrintf("HI, I am rPP %d in query : %d\n",index,iter);
#endif

    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
    Atom *atoms  = ag->atoms;
    CPNONLOCAL::eesAtmBsplineRgrp(atoms,allowedRppChares,RppData);

  }//endif : time to update the B-splines

}//end routine
//==============================================================================


//==============================================================================
// RhoGhart Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::queryCacheRHart(int index,int itime,int iter){
//==============================================================================
// Cache compute : 1st guy in does the job

#ifdef  _EESCACHE_VERBOSE_
  CkPrintf("Querying Rhart by %d at t= %d %d\n",index,itime,itimeRHart);
#endif

  if(itime != itimeRHart){

    if(itime!=itimeRHart+1 || iter != 1 || allowedRhoRHartChares[index] != 1){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Broken HartR cache query by %d at %d %d %d\n",index,itime,itimeRHart,iter);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    itimeRHart= itime;

#ifdef  _EESCACHE_VERBOSE_
    CkPrintf("Computing eesAtmBspline\n");
#endif

    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
    Atom *atoms  = ag->atoms;
    CPLOCAL::eesAtmBsplineRgrp(atoms,allowedRhoRHartChares,RhoRHartData);

  }//endif : time to update the B-splines

}//end routine
//==============================================================================


//==============================================================================
// Initialize the GParticlePlane Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void GSPDATA::init(int index_in){
//==============================================================================

  CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo; 

 //------------------------------------------------------------
 // Set the variables from the generic parainfo group

  index    = index_in;
  ngrid_a  = sim->sizeX;
  ngrid_b  = sim->sizeY;
  ngrid_c  = sim->sizeZ;
  ncoef    = sim->npts_per_chareG[index];  
  numLines = sim->nlines_per_chareG[index];
  numRuns  = 2*numLines;

 //------------------------------------------------------------
 // Set the runs from the generic group : 

  CkVec <RunDescriptor> *sortedRunDescriptors = sim->sortedRunDescriptors;
  runs  = new RunDescriptor[numRuns];
  int x = index;
  if(sortedRunDescriptors[x].size()!=numRuns){
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkPrintf("Broken run descriptor for plane %d : %d %d\n",
                x,numRuns,sortedRunDescriptors[x].size());
     CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
     CkExit();
  }//endif
  for(int j=0;j<numRuns;j++){runs[j] = sortedRunDescriptors[x][j];}

 //------------------------------------------------------------
 // set the k-vectors

  ka = (int *)fftw_malloc(ncoef*sizeof(int));
  kb = (int *)fftw_malloc(ncoef*sizeof(int));
  kc = (int *)fftw_malloc(ncoef*sizeof(int));
  genericSetKvector(ncoef,ka,kb,kc,numRuns,runs,&gCharePkg,1,ngrid_a,ngrid_b,ngrid_c);

//==============================================================================
  }//end routine
//==============================================================================


//==============================================================================
//  A generic routine to set kvectors and indices from runs
//==============================================================================
void genericSetKvector(int numPoints, int *k_x, int *k_y, int *k_z,
                       int numRuns, RunDescriptor *runs, GCHAREPKG *gCharePkg,
                       int checkFill, int ngrid_a, int ngrid_b, int ngrid_c){
//======================================================================
// Construct the k-vectors
  
  int dataCovered = 0;
  for (int r = 0; r < numRuns; r++) { // 2*number of lines z
    int x, y, z;
    x = runs[r].x;
    if (x > ngrid_a/2){x -= ngrid_a;}
    y = runs[r].y;
    if (y > ngrid_b/2){y -= ngrid_b;}
    z = runs[r].z;
    if (z > ngrid_c/2){z -= ngrid_c;}
    for(int i = 0; i < runs[r].length; i++) { //pts in lines of z
      k_x[dataCovered] = x;
      k_y[dataCovered] = y;
      k_z[dataCovered] = (z+i);
      dataCovered++;
    }//endfor
  }//endfor

  CkAssert(dataCovered == numPoints);

//======================================================================
// Find pts with k_x==0 then check the layout : kx=0 first

  int ihave_g000 =  0;
  int ind_g000   = -1;
  int ihave_kx0  = 0;
  int nkx0       = 0;
  int nkx0_uni   = 0;
  int nkx0_red   = 0;
  int nkx0_zero  = 0;
  int kx0_strt   = 0;
  for(int i=0;i<numPoints;i++){
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
  int kx0_end = kx0_strt + nkx0;

  if(checkFill==1){
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

   for(int i=0;i<nkx0;i++){  
    if(k_x[i]!=0){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("kx should be stored consecutively and first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
   }//endfor

   for(int i=0;i<nkx0_red;i++){  
    if(k_y[i]>0 || (k_y[i]==0 && k_z[i]>=0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
     }//endif
   }//endfor

   for(int i=nkx0_red;i<nkx0_uni;i++){  
    if(k_y[i]<0 || (k_y[i]==0 && k_z[i]<0)){
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("ky <0 should be stored first\n");
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
   }//endfor

  }//endif : check the ordering : states only

//==============================================================================
// Set the return values

  gCharePkg->ihave_g000 = ihave_g000;
  gCharePkg->ind_g000   = ind_g000;
  gCharePkg->ihave_kx0  = ihave_kx0;
  gCharePkg->nkx0       = nkx0;
  gCharePkg->nkx0_uni   = nkx0_uni;
  gCharePkg->nkx0_red   = nkx0_red;
  gCharePkg->nkx0_zero  = nkx0_zero;
  gCharePkg->kx0_strt   = kx0_strt;
  gCharePkg->kx0_end    = kx0_end;

//==============================================================================
  }//end routine
//==============================================================================
