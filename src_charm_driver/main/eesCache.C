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
#include "utility/util.h"
#include "cpaimd.h"
#include "eesCache.h"
#include "AtomsCache.h"
#include "src_piny_physics_v1.0/include/class_defs/Interface_ctrl.h"
#include "PhysScratchCache.h"
#include <cmath>
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_ctrl/CP_State_Plane.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpintegrate.h"


//----------------------------------------------------------------------------
extern Config config;
extern CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>            UrhoGHartExtProxy;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CProxy_PhysScratchCache         pScratchProxy;
extern CkVec <CProxy_AtomsCache>                   UatomsCacheProxy;
extern CkVec <CProxy_eesCache>                   UeesCacheProxy;

#define _EESCACHE_VERBOSE_OFF_

//==============================================================================



//==============================================================================
// Group constructor
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
eesCache::eesCache(int _nchareRPP, int _nchareGPP, int _nchareRHart,
                   int _nchareGHart, int _nstates, int _nchareRhoG, UberCollection _thisInstance): thisInstance(_thisInstance)
//==============================================================================
   {//begin rotuine
//==============================================================================

   itimeRPP        = 0;
   itimeRHart      = 0;
   rpp_on          = 0;
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

   GppData           = new GPPDATA *   [nchareGPP]; // over dimensioned
   RppData           = new RPPDATA *   [nchareRPP];
   RhoGHartData      = new RHOGHARTDATA*[nchareGHart];
   RhoRHartData      = new RHORHARTDATA*[nchareRHart];
   GspData           = new GSPDATA*[nchareGSP];
   //RhoGData          = new RHOGDATA[nchareRhoG];

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
    RppData[index] = new RPPDATA();
    RppData[index]->init(index);
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
     GppData[index] = new GPPDATA();
     GppData[index]->init(index,ncoef,ka,kb,kc);
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
    RhoRHartData[index] = new RHORHARTDATA();
    RhoRHartData[index]->init(index);
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
    RhoGHartData[index] = new RHOGHARTDATA();
    RhoGHartData[index]->init(index,ncoef,ka,kb,kc);
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
    GspData[ip] = new GSPDATA();
    GspData[ip]->init(ip);
  }//endif

  int i = nchareGSPProcT;
  if(nMallSize<=i ){
    CkPrintf("Bad Mall size in registerCacheGSP %d %d\n",nMallSize,i);
    CkExit();
  }//endif

  gspStateInd[i]  = is;
  gspPlaneInd[i]  = ip;
  nchareGSPProcT += 1;

  if(2*nchareGSPProcT >= nMallSize){
    nMallSize *= 2;
    int *tempS = new int [nMallSize];
    int *tempP = new int [nMallSize];
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
#define _CP_BRK_BETTER_OFF  // must flip cp_ees_nonlocal.h too

   CPNONLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);

   index         = index_in;
   int n_interp2 = n_interp*n_interp;

   plane_index    = (int *)fftw_malloc(natm*sizeof(int));
   nBreakJ        = (int *)fftw_malloc(natm*sizeof(int));

   igrid   = (int **)fftw_malloc(natm*sizeof(int*));
#ifdef _CP_BRK_BETTER_
   sBreakJ = (int **)fftw_malloc(natm*sizeof(int*));
#endif
   mn    = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_x = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_y = (double **)fftw_malloc(natm*sizeof(double*));
   dmn_z = (double **)fftw_malloc(natm*sizeof(double*));
   for(int i=0;i<natm;i++){
     igrid[i]    = (int *)fftw_malloc(n_interp2*sizeof(int))-1;    
#ifdef _CP_BRK_BETTER_
     sBreakJ[i]  = (int *)fftw_malloc((n_interp2+1)*sizeof(int))-1;    
#endif
     double *tmp = (double *)fftw_malloc(4*n_interp2*sizeof(double));
     int ioff    = 0;
     mn[i]       = &tmp[ioff]-1;  ioff+=n_interp2;
     dmn_x[i]    = &tmp[ioff]-1;  ioff+=n_interp2;
     dmn_y[i]    = &tmp[ioff]-1;  ioff+=n_interp2;
     dmn_z[i]    = &tmp[ioff]-1;  
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
//==============================================================================
// General intialization : Set parameters

   index          = index_in;
   rhoRsubplanes  = config.rhoRsubplanes;

   CPLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);

//====================================================
// The subplane decomposition vectors

   subStr = (int *)fftw_malloc(rhoRsubplanes*sizeof(int));
   subEnd = (int *)fftw_malloc(rhoRsubplanes*sizeof(int));
   subSiz = (int *)fftw_malloc(rhoRsubplanes*sizeof(int));
   ntemp  = (int *)fftw_malloc(rhoRsubplanes*sizeof(int));
   itemp  = (int **)fftw_malloc(rhoRsubplanes*sizeof(int*));
   for(int s=0;s<rhoRsubplanes;s++){itemp[s] = (int *)fftw_malloc(n_interp*sizeof(int))-1;}

   for(int s=0; s<rhoRsubplanes; s++){
     int div   = (ngrid_b / rhoRsubplanes); // parallelize y
     int rem   = (ngrid_b % rhoRsubplanes);
     int add   = (s < rem ? 1 : 0);
     int max   = (s < rem ? s : rem);
     subStr[s] = div*s + max;              // start of y desired by chare s
     subSiz[s] = div + add;                // total of y desired by chare s
     subEnd[s] = subStr[s] + subSiz[s];    // end   of y desired by chare s
   }//endfor

//====================================================
// Mallocs 

   plane_index = (int *)fftw_malloc(natm*sizeof(int));

   nSub  = (int     **)fftw_malloc(rhoRsubplanes*sizeof(int    *));
   igrid = (int    ***)fftw_malloc(rhoRsubplanes*sizeof(int   **));
   mn    = (double ***)fftw_malloc(rhoRsubplanes*sizeof(double**));
   dmn_x = (double ***)fftw_malloc(rhoRsubplanes*sizeof(double**));
   dmn_y = (double ***)fftw_malloc(rhoRsubplanes*sizeof(double**));
   dmn_z = (double ***)fftw_malloc(rhoRsubplanes*sizeof(double**));

   int n_interp21 = n_interp*n_interp+1;
   for(int s=0;s<rhoRsubplanes;s++){
     nSub[s]  = (int    *)fftw_malloc (natm*sizeof(int   *));
     igrid[s] = (int    **)fftw_malloc(natm*sizeof(int   *));
     mn[s]    = (double **)fftw_malloc(natm*sizeof(double*));
     dmn_x[s] = (double **)fftw_malloc(natm*sizeof(double*));
     dmn_y[s] = (double **)fftw_malloc(natm*sizeof(double*));
     dmn_z[s] = (double **)fftw_malloc(natm*sizeof(double*));
     for(int i=0;i<natm;i++){
       igrid[s][i] = (int    *)fftw_malloc(n_interp21*sizeof(int   ))-1;
       mn[s][i]    = (double *)fftw_malloc(n_interp21*sizeof(double))-1;
       dmn_x[s][i] = (double *)fftw_malloc(n_interp21*sizeof(double))-1;
       dmn_y[s][i] = (double *)fftw_malloc(n_interp21*sizeof(double))-1;
       dmn_z[s][i] = (double *)fftw_malloc(n_interp21*sizeof(double))-1;
     }//endfor : atoms
   }//endfor : subplanes

//-----------------------------------------------------------------------------
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

    AtomsCache *ag = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    CkAssert(ag!=NULL);
    FastAtoms *fastAtoms = &(ag->fastAtoms);
    CkAssert(fastAtoms!=NULL);
    CkAssert(fastAtoms->x!=NULL);
#if CMK_TRACE_ENABLED
   double  StartTime=CmiWallTimer();
#endif    

   CPNONLOCAL::eesAtmBsplineRgrp(fastAtoms,allowedRppChares,RppData, pScratchProxy.ckLocalBranch()->psscratch);

#if CMK_TRACE_ENABLED
  traceUserBracketEvent(eesAtmBspline_, StartTime, CmiWallTimer());    
#endif

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

    AtomsCache *ag = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    FastAtoms *fastAtoms = &(ag->fastAtoms);
#if CMK_TRACE_ENABLED
    double  StartTime=CmiWallTimer();
#endif    

    CPLOCAL::eesAtmBsplineRgrp(fastAtoms,allowedRhoRHartChares,RhoRHartData, pScratchProxy.ckLocalBranch()->psscratch);

#if CMK_TRACE_ENABLED
    traceUserBracketEvent(eesAtmBspline_, StartTime, CmiWallTimer());    
#endif

  }//endif : time to update the B-splines

}//end routine
//==============================================================================


//==============================================================================
// Initialize the GstatePlane Cache Data class
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

  g  = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  g2 = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  ka = (int *)fftw_malloc(ncoef*sizeof(int));
  kb = (int *)fftw_malloc(ncoef*sizeof(int));
  kc = (int *)fftw_malloc(ncoef*sizeof(int));

  CPNONLOCAL::genericSetKvector(ncoef,ka,kb,kc,g,g2,numRuns,runs,&gCharePkg,1,
                                ngrid_a,ngrid_b,ngrid_c);

 //------------------------------------------------------------
 // set the coef masses

  int mydoublePack = config.doublePack;
  coef_mass        = (double *)fftw_malloc(ncoef*sizeof(double));
  CPINTEGRATE::CP_create_mass(ncoef,ka,kb,kc,coef_mass,mydoublePack);

//==============================================================================
  }//end routine
//==============================================================================
