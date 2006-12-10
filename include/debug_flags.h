//=============================================================================
//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
//=============================================================================
//                    Useful debugging flags
//=============================================================================

//=============================================================================
// Major League Debug Controllers

//-----------------------------------------------------------------------
//#define  _CP_DEBUG_SCALC_ONLY_  // REALLY scalc only NOTHING ELSE
#ifdef _CP_DEBUG_SCALC_ONLY_
#define _CP_DEBUG_SFNL_OFF_
#define _CP_DEBUG_VKS_OFF_
#define _CP_DEBUG_RHO_OFF_  // not really necesary as no path to rho
#endif

//-----------------------------------------------------------------------
//#define  _CP_DEBUG_VKS_ONLY_  //what we called scalc only : states make rho but don't send
#ifdef _CP_DEBUG_VKS_ONLY_
#define _CP_DEBUG_SFNL_OFF_
#define _CP_DEBUG_RHO_OFF_     // blocks send to rho from state chares
#endif

//-----------------------------------------------------------------------
//#define  _CP_DEBUG_NON_LOCAL_ONLY_  // SFNL + scalc, no density creation 
#ifdef _CP_DEBUG_NON_LOCAL_ONLY_
#define _CP_DEBUG_VKS_OFF_
#define _CP_DEBUG_RHO_OFF_  // not really necessary as there is no path to rho
#endif

//-----------------------------------------------------------------------
//#define _CP_DEBUG_HARTEEXT_OFF_   // this leaves everything on but eext/hartree

//=============================================================================
// src_charm_driver/paircalc/ckPairCalculator.C and a dozen other places
//#define _NAN_CHECK_

//=============================================================================
// src_charm_driver/main/groups.C
//#define _CP_DEBUG_PSI_OFF_

//=============================================================================
// src_charm_driver/fft_slab_ctrl
//#define _CP_DEBUG_NON_LOCAL_ONLY_
//#define _CP_DEBUG_FFTR_VKSR_

//=============================================================================
// src_charm_driver/utilities/util.C
//#define _CP_DEBUG_UTIL_VERBOSE_
//#define _CP_DEBUG_LINE_

//=============================================================================
// include/CPcharmParaInfo.h
//#define _CP_DEBUG_PARAINFO_VERBOSE_

//=============================================================================
// src_charm_driver/cp_density_ctrl/CP_Rho_RealSpacePlane.C : 
//#define _CP_DEBUG_RHOR_VERBOSE_
//#define _CP_DEBUG_RHOR_RHO_
//#define _CP_DEBUG_RHOR_VKSA_
//#define _CP_DEBUG_RHOR_VKSB_
//#define _CP_DEBUG_RHOR_VKSC_
//#define _CP_DEBUG_RHOR_VKSD_
//#define _CP_DEBUG_RHOR_VKSE_
//#define _CP_DEBUG_HARTEEXT_OFF_

//=============================================================================
// src_charm_driver/cp_density_ctrl/CP_Rho_GSpacePlane.C : 
//#define _CP_DEBUG_VKS_GSPACE_
//#define _CP_DEBUG_VKS_FORCES_
//#define _CP_DEBUG_RHOG_VERBOSE_
//#define _CP_DEBUG_RHOG_RHOG_
//#define _CP_DEBUG_RHOG_VKSA_
//#define _CP_DEBUG_HARTEEXT_OFF_

//=============================================================================
// src_charm_driver/cp_state_ctrl/CP_State_GSpacePlane.C    
//#define _CP_DEBUG_NON_LOCAL_ONLY_
//#define _CP_DEBUG_OLDFORCE_
//#define _CP_DEBUG_LMAT_
//#define _CP_DEBUG_LMAT_
//#define _CP_DEBUG_STATEG_VERBOSE_
//#define _CP_DEBUG_DYNAMICS_
//#define _CP_DEBUG_ORTHO_OFF_
//#define _CP_DEBUG_SFNL_OFF_
//#define _CP_DEBUG_PSI_OFF_
#define _CP_DEBUG_COEF_SCREEN_
#define _CP_DEBUG_UPDATE_OFF_
//#define _CP_DEBUG_NONLOC_BARRIER_

//=============================================================================
// src_charm_driver/cp_state_ctrl/CP_State_RealSpacePlane.C : 
//#define _CP_DEBUG_STATER_VKS_
//#define _CP_DEBUG_STATER_VERBOSE_
//#define _CP_DEBUG_RHO_OFF_

//=============================================================================
// src_charm_driver/structure_factor/StructureFactor.C : 
//#define _CP_DEBUG_SF_CALC_

//=============================================================================
// src_charm_driver/structure_factor/StructFactorCache.C : 
//#define _CP_DEBUG_SF_CACHE_

//=============================================================================
// src_charm_driver/orthog_ctrl/ortho.C
//#define _CP_DEBUG_SMAT_
//#define _CP_DEBUG_TMAT_
//#define _CP_ORTHO_DUMP_SMAT_
//#define _CP_ORTHO_DUMP_TMAT_
//#define _CP_ORTHO_DEBUG_COMPARE_SMAT_
//#define _CP_ORTHO_DEBUG_COMPARE_TMAT_ 

//=============================================================================
// src_charm_driver/load_balance/MapTable.C
//#define _MAP_DEBUG_
//#define _MAP_VERBOSE_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_orthog/cp_lowdin.C
//#define _CP_DEBUG_NODIAG_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_integrate/cp_min_std.C
//#define _CP_DEBUG_NEWFORCE_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_ions/cp_rspace_ion.C
// src_piny_physics_v1.0/classical_physics/integrate/control_integrate.C
//#define _CP_DEBUG_ATM_FORC_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_local/cp_hart_ext.C
//#define _CP_DEBUG_VKS_HART_EEXT_

//=============================================================================
// src_piny_physics_v1.0/abinitio_physics/cp_nonlocal/cp_eke.C
//#define _CP_DEBUG_OLDFORCE_

//=============================================================================
// src_piny_physics_v1.0/piny_to_driver/Parainfoinit.C
//#define _CP_DEBUG_LINES_PER_PLANE_
