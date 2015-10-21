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
extern CPcharmParaInfo                      simReadOnly;
extern CkVec <CProxy_CP_State_ParticlePlane>     UparticlePlaneProxy;
extern CkVec <CProxy_CP_State_RealParticlePlane> UrealParticlePlaneProxy;
extern CkVec <CProxy_CP_Rho_RHartExt>            UrhoRHartExtProxy;
extern CkVec <CProxy_CP_Rho_GHartExt>            UrhoGHartExtProxy;
extern CkVec <CProxy_PhysScratchCache>           UpScratchProxy;
extern CkVec <CProxy_AtomsCache>                 UatomsCacheProxy;
extern CkVec <CProxy_eesCache>                   UeesCacheProxy;

#define _EESCACHE_VERBOSE_OFF_

//==============================================================================



//==============================================================================
// Group constructor
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
eesCache::eesCache(int _nchareRPP, int _nchareGPP, int _nstates, int _nkpoint,
    UberCollection _thisInstance): thisInstance(_thisInstance)
                                   //==============================================================================
{//begin rotuine
  //==============================================================================

  itimeRPP        = 0;
  itimeRHart      = 0;
  rpp_on          = 0;
  nMallSize       = 100;

  nkpoint         = _nkpoint;
  nchareRPP       = _nchareRPP;
  nchareGPP       = _nchareGPP;
  nchareGSP       = _nchareGPP;
  nstates         = simReadOnly.nstates;

  nchareGSPProcT  = 0;
  nchareGSPProc   = 0;
  nchareRPPProc   = 0;
  nchareGPPProc   = 0;

  allowedRppChares      = new int[nchareRPP];
  allowedGppChares      = new int[nchareGPP];
  allowedGspChares      = new int[nchareGSP];

  for(int i=0;i<nchareRPP  ;i++){allowedRppChares[i]      = 0;}
  for(int i=0;i<nchareGPP  ;i++){allowedGppChares[i]      = 0;}
  for(int i=0;i<nchareGSP;  i++){allowedGspChares[i]      = 0;}

  indGppChares      = new int[nchareGPP];  // over dimensioned
  indRppChares      = new int[nchareRPP];
  indGspChares      = new int[nchareGSP];

  GppData           = new GPPDATA *   [nchareGPP]; // over dimensioned
  RppData           = new RPPDATA *   [nchareRPP];
  GspData           = new GSPDATA*[nchareGSP];

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
    GppData[index]->init(nkpoint,index,ncoef,ka,kb,kc);
  }//endif

}//end routine
//==============================================================================


//==============================================================================
// RhoRhart Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int eesCache::registerCacheRHart(int xindex, int yindex, int gridc_start,
    int gridc_end, int gridb_start, int gridb_end) {
  //==============================================================================

  for(int i = 0; i < RhoRHartData.size(); i++) {
    if(RhoRHartData[i].xindex == xindex && RhoRHartData[i].yindex == yindex) {
      return i;
    }
  }

  int returnIndex = RhoRHartData.size();
  RHORHARTDATA tempData;
  tempData.xindex = xindex;
  tempData.yindex = yindex;
  tempData.gridc_start = gridc_start;
  tempData.gridb_start = gridb_start;
  tempData.gridc_end = gridc_end;
  tempData.gridb_end = gridb_end;
  tempData.mygridc = gridc_end - gridc_start;
  tempData.mygridb = gridb_end - gridb_start;
  RhoRHartData.push_back(tempData);
  RhoRHartData[returnIndex].init();
  return returnIndex;

}//end routine
//==============================================================================


//==============================================================================
// RhoGhart Cache Management tool
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
int eesCache::registerCacheGHart(int index, std::vector< gridPoint > * myPoints)
{
  //==============================================================================
  for(int i = 0; i < RhoGHartData.size(); i++) {
    if(RhoGHartData[i].index == index) {
      return i;
    }
  }
  int returnIndex = RhoGHartData.size();
  RHOGHARTDATA tempData;
  tempData.index = index;
  tempData.points = myPoints;
  RhoGHartData.push_back(tempData);
  RhoGHartData[returnIndex].init();
  return returnIndex;
}//end routine
//==============================================================================


//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void eesCache::registerCacheGSP(int kpt, int spin, int is ,int ip){

  if(allowedGspChares[ip]==0){
    nchareGSPProc       += 1;
    allowedGspChares[ip] = 1;
    GspData[ip] = new GSPDATA();
    GspData[ip]->init(ip,nkpoint);
  }//endif

  int i = nchareGSPProcT;

  gspKptSpinStatePlaneVec.push_back(KthSpinStatePlaneTuple(kpt, spin, is, ip));
  //  CkPrintf("[%d] register kpt %d spin %d state %d plane %d\n", CkMyPe(), kpt, spin, is, ip);
  nchareGSPProcT += 1;

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
void GPPDATA::init(int nkpoint_in,int index_in,int ncoef_in, int *ka, int *kb, int *kc){

  nkpoint  = nkpoint_in;
  index    = index_in;
  ncoef    = ncoef_in;
  b_re     = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  b_im     = (double *)fftw_malloc((ncoef+1)*sizeof(double));

  h_gspl    = cmall_mat(0,nkpoint,0,ncoef,"CPPDATA:init");
  ind_gspl  = cmall_int_mat(0,nkpoint,0,ncoef,"CPPDATA:init");

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
void RHORHARTDATA::init() {
  //==============================================================================
  // General intialization : Set parameters

  CPLOCAL::getEesPrms(&ngrid_a, &ngrid_b, &ngrid_c, &n_interp, &natm);

  //====================================================
  // Mallocs
  itemp  = (int     **)fftw_malloc(mygridc * sizeof(int    *));
  nSub   = (int     **)fftw_malloc(mygridc * sizeof(int    *));
  igrid  = (int    ***)fftw_malloc(mygridc * sizeof(int   **));
  mn     = (double ***)fftw_malloc(mygridc * sizeof(double**));
  dmn_x  = (double ***)fftw_malloc(mygridc * sizeof(double**));
  dmn_y  = (double ***)fftw_malloc(mygridc * sizeof(double**));
  dmn_z  = (double ***)fftw_malloc(mygridc * sizeof(double**));
  plane_index = (int **)fftw_malloc(mygridc * sizeof(int*));

  int n_interp21 = n_interp * n_interp + 1;
  for(int c = 0; c < mygridc; c++) {
    plane_index[c] = (int *)fftw_malloc(natm * sizeof(int));
    nSub[c]  = (int     *)fftw_malloc(natm * sizeof(int));
    itemp[c] = (int     *)fftw_malloc(n_interp * sizeof(int)) - 1;
    igrid[c] = (int    **)fftw_malloc(natm * sizeof(int *));
    mn[c]    = (double **)fftw_malloc(natm * sizeof(double*));
    dmn_x[c] = (double **)fftw_malloc(natm * sizeof(double*));
    dmn_y[c] = (double **)fftw_malloc(natm * sizeof(double*));
    dmn_z[c] = (double **)fftw_malloc(natm * sizeof(double*));
    for(int i = 0; i < natm; i++){
      igrid[c][i] = (int    *)fftw_malloc(n_interp21 * sizeof(int   )) - 1;
      mn[c][i]    = (double *)fftw_malloc(n_interp21 * sizeof(double)) - 1;
      dmn_x[c][i] = (double *)fftw_malloc(n_interp21 * sizeof(double)) - 1;
      dmn_y[c][i] = (double *)fftw_malloc(n_interp21 * sizeof(double)) - 1;
      dmn_z[c][i] = (double *)fftw_malloc(n_interp21 * sizeof(double)) - 1;
    }//endfor : atoms
  }
  //-----------------------------------------------------------------------------
}//end routine
//==============================================================================


//==============================================================================
// Initialize the RhoGhart Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RHOGHARTDATA::init() {

  ncoef    = (*points).size();
  b_re     = (double *)fftw_malloc((ncoef+1)*sizeof(double));
  b_im     = (double *)fftw_malloc((ncoef+1)*sizeof(double));

  CPLOCAL::getEesPrms(&ngrid_a, &ngrid_b, &ngrid_c, &n_interp, &natm);
  CPLOCAL::eesSetEesWghtGgrp(ncoef, (*points), b_re, b_im, ngrid_a, ngrid_b,
      ngrid_c, n_interp);
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

    CPNONLOCAL::eesAtmBsplineRgrp(fastAtoms,allowedRppChares,RppData, UpScratchProxy[thisInstance.proxyOffset].ckLocalBranch()->psscratch);

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
void eesCache::queryCacheRHart(int index, int itime, int iter){
  //==============================================================================
  // Cache compute : 1st guy in does the job

#ifdef  _EESCACHE_VERBOSE_
  CkPrintf("Querying Rhart by %d at t= %d %d\n", index, itime, itimeRHart);
#endif

  if(itime != itimeRHart) {

    if(itime != itimeRHart+1 || iter != 1) {
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkPrintf("Broken RHart cache query by %d at %d %d %d\n",index, itime,
          itimeRHart, iter);
      CkPrintf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@\n");
      CkExit();
    }//endif
    itimeRHart = itime;

#ifdef  _EESCACHE_VERBOSE_
    CkPrintf("Computing eesAtmBspline\n");
#endif

    AtomsCache *ag = UatomsCacheProxy[thisInstance.proxyOffset].ckLocalBranch();
    FastAtoms *fastAtoms = &(ag->fastAtoms);

#if CMK_TRACE_ENABLED
    double  StartTime=CmiWallTimer();
#endif

    CPLOCAL::eesAtmBsplineRgrp(fastAtoms, RhoRHartData,
        UpScratchProxy[thisInstance.proxyOffset].ckLocalBranch()->psscratch);

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
void GSPDATA::init(int index_in,int nkpoint_in){
  //==============================================================================

  CPcharmParaInfo *sim = CPcharmParaInfo::get();

  //------------------------------------------------------------
  // Set the variables from the generic parainfo group

  index    = index_in;
  ngrid_a  = sim->sizeX;
  ngrid_b  = sim->sizeY;
  ngrid_c  = sim->sizeZ;
  ncoef    = sim->npts_per_chareG[index];
  numLines = sim->nlines_per_chareG[index];
  numRuns  = 2*numLines;
  nkpoint  = nkpoint_in;

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
  g2 = cmall_mat(0,nkpoint,0,ncoef,"eesCache.C");
  g  = cmall_mat(0,nkpoint,0,ncoef,"eesCache.C");

  CPNONLOCAL::genericSetKvector(ncoef,ka,kb,kc,numRuns,runs,&gCharePkg,1,
      ngrid_a,ngrid_b,ngrid_c,g2,g);

  //------------------------------------------------------------
  // set the coef masses

  int mydoublePack = config.doublePack;
  coef_mass        = (double *)fftw_malloc(ncoef*sizeof(double));
  CPINTEGRATE::CP_create_mass(ncoef,ka,kb,kc,coef_mass,mydoublePack);

  //==============================================================================
}//end routine
//==============================================================================
