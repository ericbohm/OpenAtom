/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_VanderWaalsR.C
 * This is the description of the "life" of a CP_VanderWaals object.
 *
 * --Marcello please insert a description--
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
extern CkVec <CProxy_CP_Rho_RealSpacePlane>      UrhoRealProxy;
extern CkVec <CProxy_CP_VanderWaalsR>      UVdWRealProxy;
extern CkVec <CProxy_CP_VanderWaalsG>      UVdWGProxy;
extern CkVec <CProxy_FFTcache>                   UfftCacheProxy;
extern CkVec <CProxy_GSpaceDriver>               UgSpaceDriverProxy;

extern CkGroupID            mCastGrpId;

extern Config    config;
extern int       nstates;


//#define _CP_DEBUG_VDWR_VERBOSE_

/** @addtogroup VanderWaals
    @{
*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//
// This class (array) accepts the real space densities from all the states
// Performs lots of operations to compute the VanderWaals forces
//
//============================================================================
CP_VanderWaalsR::CP_VanderWaalsR(UberCollection _instance) :
  thisInstance(_instance)
//============================================================================
   {//begin routine
//============================================================================

#ifdef _CP_DEBUG_VDWR_VERBOSE_
    CkPrintf("[%d %d %d] VdWR constructs \n",thisIndex.x, thisIndex.y, thisIndex.z);
#endif

//============================================================================
   }//end routine
//============================================================================

/** post constructor initialization */
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_VanderWaalsR::init(){

// make sections or other fancy post launch stuff
}

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_VanderWaalsR::~CP_VanderWaalsR(){
}
//============================================================================


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
//  Pup my variables for migration
//============================================================================
void CP_VanderWaalsR::pup(PUP::er &p){
  ArrayElement3D::pup(p);
  //  p|thisInstance; # how do you pup a const?
//---------------------------------------------------------------------------
   }//end routine
//============================================================================
/*@}*/
