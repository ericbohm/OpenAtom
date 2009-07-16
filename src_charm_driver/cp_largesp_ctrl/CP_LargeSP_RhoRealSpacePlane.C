/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_LargeSP_RhoRealSpacePlane.C
 * This is the description of the "life" of a CP_LargeSP_RhoRealSpacePlane object.
 *
 * At the start of the program, the constructor CP_LargeSP_RhoRealSpacePlane() is called.
 * 
 * The CP_Rho_RealSpacePlanes interpolate onto the larger grid by
 * sending their charge to CP_LargeSP_RhoRealSpacePlane using the
 * acceptRhoR() method. 
 *
 * CP_LargeSP_RhoRealSpacePlane fft transposes this data to CP_LargeSP_RhoGSpacePlane.
 * When CP_LargeSP_RhoGSpacePlane is done interacting with NAMD it fft
 * transposes the updated grid back to CP_LargeSP_RhoRealSpacePlane
 * via acceptLSPRhoG.  CP_LargeSP_RhoRealSpacePlane then sends the
 * updated charge grad back to CP_Rho_RealSpacePlane for VKS et al.
 */
//============================================================================

#include "charm++.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "main/groups.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_ctrl/CP_State_Plane.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================
extern CProxy_TimeKeeper                 TimeKeeperProxy;
extern CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
extern CkVec <CProxy_CP_LargeSP_RhoGSpacePlane>      UlsRhoGProxy;
extern CProxy_CPcharmParaInfoGrp         scProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;


extern CkGroupID            mCastGrpId;

extern Config    config;
extern int       nstates;

bool is_pow2(int );



//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the Rho_Reals.
// Interpolates onto the larger grid, FFT transpose to LargeSP_RhoG,
// receive the ifft from LargeSP_RhoG and ungrid back to Rho_Reals.
//
//============================================================================
CP_LargeSP_RhoRealSpacePlane::CP_LargeSP_RhoRealSpacePlane(
					     UberCollection _instance) :
  thisInstance(_instance)
//============================================================================
   {//begin routine
//============================================================================

#ifdef _CP_DEBUG_LARGESP_RHOREAL_VERBOSE_
    CkPrintf("[%d %d] LargeSP_RhoReal constructs \n",thisIndex.x, thisIndex.y);
#endif

//============================================================================
// Get parameters from the globals/groups

    CPcharmParaInfo *sim = (scProxy.ckLocalBranch ())->cpcharmParaInfo;



//============================================================================
// handle class member local initialization



//============================================================================
// Migration

    usesAtSync = CmiTrue;
    setMigratable(false);


//---------------------------------------------------------------------------
   }//end routine
//============================================================================

/** post constructor initialization */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_LargeSP_RhoRealSpacePlane::init(){

// place holder

//---------------------------------------------------------------------------
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_LargeSP_RhoRealSpacePlane::~CP_LargeSP_RhoRealSpacePlane(){

//---------------------------------------------------------------------------
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//  Pup my variables for migration
//============================================================================
void CP_LargeSP_RhoRealSpacePlane::pup(PUP::er &p){
  ArrayElement2D::pup(p);
//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Data comes from Rho_Real once an algorithm step.
/** 
 * Here the density is interpolated onto the large grid.
 */
//============================================================================
void CP_LargeSP_RhoRealSpacePlane::acceptRhoR() {
//============================================================================

#ifdef _CP_DEBUG_LARGE_RHOREAL_VERBOSE_
    CkPrintf("[%d,%d] LSP_RhoReal accepting Density %d \n",
              thisIndex.x,thisIndex.y,CkMyPe());
#endif



//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Data comes from LSP_RhoG once an algorithm step.
/** 
 * Here the density is interpolated onto the large grid.
 */
//============================================================================
void CP_LargeSP_RhoRealSpacePlane::acceptLSPRhoG() {
//============================================================================

#ifdef _CP_DEBUG_LARGE_RHOREAL_VERBOSE_
    CkPrintf("[%d, %d] LSP_RhoReal accepting LSPRhoG %d\n",
              thisIndex.x,thisIndex.y,CkMyPe());
#endif



//----------------------------------------------------------------------------
  }//end routine
//============================================================================


