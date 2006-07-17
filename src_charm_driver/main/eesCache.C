//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
/** \file groups.C
 * 
 *           Processor group class Functions : Atoms and parainfo
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
//==============================================================================

#define _EESCACHE_VERBOSE_OFF_


//==============================================================================
// Group constructor
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
eesCache::eesCache(int _nchareRPP, int _nchareGPP, int _nchareRHart,
                   int _nchareGHart)
//==============================================================================
   {//begin rotuine
//==============================================================================

   itimeRPP        = 0;
   itimeRHart      = 0;

   nchareRPP       = _nchareRPP;
   nchareGPP       = _nchareGPP;
   nchareRHart     = _nchareRHart;
   nchareGHart     = _nchareGHart;
   nchareRPPProc   = 0;
   nchareGPPProc   = 0;
   nchareRHartProc = 0;
   nchareGHartProc = 0;

   allowedRppChares      = new int[nchareRPP];
   allowedGppChares      = new int[nchareGPP];
   allowedRhoRHartChares = new int[nchareRHart];
   allowedRhoGHartChares = new int[nchareGHart];
   for(int i=0;i<nchareRPP;i++)  {allowedRppChares[i]      = 0;}
   for(int i=0;i<nchareGPP;i++)  {allowedGppChares[i]      = 0;}
   for(int i=0;i<nchareRHart;i++){allowedRhoRHartChares[i] = 0;}
   for(int i=0;i<nchareGHart;i++){allowedRhoGHartChares[i] = 0;}

   indGppChares      = new int[nchareGPP];  // over dimensioned
   indRppChares      = new int[nchareRPP];
   indRhoRHartChares = new int[nchareRHart];
   indRhoGHartChares = new int[nchareGHart];

   GppData           = new GPPDATA     [nchareGPP]; // over dimensioned
   RppData           = new RPPDATA     [nchareRPP];
   RhoGHartData      = new RHOGHARTDATA[nchareGHart];
   RhoRHartData      = new RHORHARTDATA[nchareRHart];

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
    nchareRHartProc             += 1;
    allowedRhoRHartChares[index] = 1;
#ifdef  _EESCACHE_VERBOSE_
    CkPrintf("Registering Rhart %d\n",index);
#endif
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
// Initialize the RealParticlePlane Cache Data class
//==============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//==============================================================================
void RPPDATA::init(int index_in){

   CPNONLOCAL::getEesPrms(&ngrid_a,&ngrid_b,&ngrid_c,&n_interp,&natm);

   index          = index_in;
   int n_interp21 = n_interp*n_interp+1; // extra space for piny

   plane_index    = new int[natm];

   igrid = new int    *[natm];
   mn    = new double *[natm];
   dmn_x = new double *[natm];
   dmn_y = new double *[natm];
   dmn_z = new double *[natm];
   for(int i=0;i<natm;i++){
     igrid[i] = new int    [n_interp21];
     mn[i]    = new double [n_interp21];
     dmn_x[i] = new double [n_interp21];
     dmn_y[i] = new double [n_interp21];
     dmn_z[i] = new double [n_interp21];
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
  b_re     = new double[ncoef+1];
  b_im     = new double[ncoef+1];
  h_gspl   = new double[ncoef+1];
  ind_gspl = new int[ncoef+1];

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

   plane_index    = new int[natm];

   igrid = new int    *[natm];
   mn    = new double *[natm];
   dmn_x = new double *[natm];
   dmn_y = new double *[natm];
   dmn_z = new double *[natm];
   for(int i=0;i<natm;i++){
     igrid[i] = new int    [n_interp21];
     mn[i]    = new double [n_interp21];
     dmn_x[i] = new double [n_interp21];
     dmn_y[i] = new double [n_interp21];
     dmn_z[i] = new double [n_interp21];
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
  b_re     = new double[ncoef+1];
  b_im     = new double[ncoef+1];

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

    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
    Atom *atoms  = ag->atoms;
    CkPrintf("HI, I am rPP %d in query : %d\n",index,iter);
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

    AtomsGrp *ag = atomsGrpProxy.ckLocalBranch();
    Atom *atoms  = ag->atoms;
#ifdef  _EESCACHE_VERBOSE_
    CkPrintf("Computing eesAtmBspline\n");
#endif
    CPLOCAL::eesAtmBsplineRgrp(atoms,allowedRhoRHartChares,RhoRHartData);

  }//endif : time to update the B-splines

}//end routine
//==============================================================================
