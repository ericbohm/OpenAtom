/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_LargeSP_RhoGSpacePlane.C
 * @defgroup LargeSparse LargeSparse
 * This is the description of the "life" of a CP_LargeSP_RhoGSpacePlane object.
 *
 * At the start of the program, the constructor CP_LargeSP_RhoGSpacePlane() is called.
 * 
 * CP_LargeSP_RhoGSpacePlanes receive the FFT transpose from
 * CP_LargeSP_RhoRealSpacePlane.  They then receive S(g) from the
 * large box MD code (NAMD) and combine that with the charge grid to
 * compute the QM force contribution.  The QM force contribution is
 * then returned to the MD code (using a ckcallback?). Finally, the
 * grid is FFT transposed back to CP_LargeSP_RhoRealSpacePlane.
 */
//============================================================================

#include "charm++.h"
#include <iostream>
#include <fstream>
#include <cmath>

#include "debug_flags.h"
#include "utility/util.h"
#include "main/cpaimd.h"
#include "fft_slab_ctrl/fftCacheSlab.h"
#include "cp_state_ctrl/CP_State_Plane.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================
extern CProxy_TimeKeeper                 TimeKeeperProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;


extern CkGroupID            mCastGrpId;

extern Config    config;
extern int       nstates;

bool is_pow2(int );


/** @addtogroup LargeSparse
    @{
*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the Rho_Gs.
// Interpolates onto the larger grid, FFT transpose to LargeSP_RhoG,
// receive the ifft from LargeSP_RhoG and ungrid back to Rho_Gs.
//
//============================================================================
CP_LargeSP_RhoGSpacePlane::CP_LargeSP_RhoGSpacePlane(
					     UberCollection _instance) :
  thisInstance(_instance)
//============================================================================
   {//begin routine
//============================================================================

#ifdef _CP_DEBUG_LARGESP_RHOG_VERBOSE_
    CkPrintf("[%d] LargeSP_RhoG constructs \n",thisIndex.x);
#endif

//============================================================================
// Get parameters from the globals/groups

    CPcharmParaInfo *sim = CPcharmParaInfo::get();



//============================================================================
// handle class member local initialization



//============================================================================
// Migration

    usesAtSync = true;
    setMigratable(false);


//---------------------------------------------------------------------------
   }//end routine
//============================================================================

/** post constructor initialization */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_LargeSP_RhoGSpacePlane::init(){

// place holder

//---------------------------------------------------------------------------
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_LargeSP_RhoGSpacePlane::~CP_LargeSP_RhoGSpacePlane(){

//---------------------------------------------------------------------------
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//  Pup my variables for migration
//============================================================================
void CP_LargeSP_RhoGSpacePlane::pup(PUP::er &p){
  ArrayElement1D::pup(p);
//---------------------------------------------------------------------------
   }//end routine
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Data comes from LargeSP_RhoReal once an algorithm step.
/** 
 * Here the density is interpolated onto the large grid.
 */
//============================================================================
void CP_LargeSP_RhoGSpacePlane::acceptLSPRhoR() {
//============================================================================

#ifdef _CP_DEBUG_LARGE_RHOREAL_VERBOSE_
    CkPrintf("LSP_RhoG accepting LSPRhoR %d %d\n",
              thisIndex.x,CkMyPe());
#endif



//----------------------------------------------------------------------------
  }//end routine
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
//  Data comes from MD once an algorithm step.
/** 
 * Here the density is interpolated onto the large grid.
 */
//============================================================================
void CP_LargeSP_RhoGSpacePlane::acceptMDSg() {
//============================================================================

#ifdef _CP_DEBUG_LARGE_MDSG_VERBOSE_
    CkPrintf("LSP_RhoG accepting MDSg %d %d \n",
              thisIndex.x,CkMyPe());
#endif



//----------------------------------------------------------------------------
  }//end routine
//============================================================================

/*@}*/
