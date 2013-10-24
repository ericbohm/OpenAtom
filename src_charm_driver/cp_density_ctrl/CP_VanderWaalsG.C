//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
/** \file CP_VanderWaalsG.C
 * @defgroup VanderWaals VanderWaals
 *  This is a description of the "life" of a CP_VanderWaalsG  object
 *
 * --Marcello please insert a description--
 *
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

#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cplocal.h"
#include "src_piny_physics_v1.0/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h"

//============================================================================

extern Config                               config;
extern CkVec <CProxy_CP_Rho_RealSpacePlane> UrhoRealProxy;
extern CkVec <CProxy_CP_Rho_GSpacePlane> UrhoGProxy;
extern CkVec <CProxy_CP_VanderWaalsR> UVdWRealProxy;
extern CkVec <CProxy_GSpaceDriver>          UgSpaceDriverProxy;
extern CkVec <CProxy_FFTcache>           UfftCacheProxy;

//#define _CP_DEBUG_VDWG_VERBOSE_

//============================================================================

/** @addtogroup VanderWaals
    @{
*/

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_VanderWaalsG::CP_VanderWaalsG(UberCollection _instance) :
  thisInstance(_instance)
//============================================================================
    {//begin routine
//============================================================================
#ifdef _CP_DEBUG_VDWG_VERBOSE_
    CkPrintf("[%d %d] VdW GS constructor\n",thisIndex.x,thisIndex.y);
#endif
//============================================================================
// Set counters local variables

    CPcharmParaInfo *sim = CPcharmParaInfo::get();
//---------------------------------------------------------------------------
   }//end routine
//============================================================================

//post constructor initialization
//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_VanderWaalsG::init()
{
    CPcharmParaInfo *sim = CPcharmParaInfo::get();

}


//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
CP_VanderWaalsG::~CP_VanderWaalsG(){

}
//============================================================================

//============================================================================
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//============================================================================
void CP_VanderWaalsG::pup(PUP::er &p){
  //  ArrayElement2D::pup(p);
  //  p|thisInstance;
//--------------------------------------------------------------------------
   }//end routine
//============================================================================

/*@}*/
