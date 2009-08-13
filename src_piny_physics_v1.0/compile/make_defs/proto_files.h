#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
#
#  Abreviations for include files  PI_MD.
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================

#-------------------------------------------------------
#              STANDARD INCLUDE FILE

STANDARD      = standard_include.h
DEFINES       = $(physics)/include/class_defs/piny_constants.h 

#-------------------------------------------------------
#              TYPEDEFS 

CLASS_GEN     = $(physics)/include/class_defs/allclass_gen.h
CLASS_MDINT   = $(physics)/include/class_defs/allclass_mdintegrate.h
CLASS_MDATM   = $(physics)/include/class_defs/allclass_mdatoms.h
CLASS_MDINTER = $(physics)/include/class_defs/allclass_mdinter.h
CLASS_MDINTRA = $(physics)/include/class_defs/allclass_mdintra.h
CLASS_CP      = $(physics)/include/class_defs/allclass_cp.h
TYP_STAT      = $(physics)/include/class_defs/typedefs_stat.h
TYP_PAR       = $(physics)/include/class_defs/typedefs_par.h
DBG_FLAGS     = $(physics)/../include/debug_flags.h

CLASS_PhysicsAtomPosinit = $(physics)/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h
CLASS_PARAINIT  = $(physics)/include/class_defs/PINY_INIT/PhysicsParamTrans.h
CLASS_CP_NONLOC = $(physics)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h
CLASS_CP_LOC    = $(physics)/include/class_defs/CP_OPERATIONS/class_cplocal.h
CLASS_CP_XCFNCTLS = $(physics)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
CLASS_CP_ORTHOG   = $(physics)/include/class_defs/CP_OPERATIONS/class_cporthog.h
CLASS_CP_INTEGRATE = $(physics)/include/class_defs/CP_OPERATIONS/class_cpintegrate.h
CLASS_CP_RSPACEION = $(physics)/include/class_defs/CP_OPERATIONS/class_cprspaceion.h
CLASS_CHARM_ATOM = $(physics)/../include/Atoms.h
CLASS_ATOM_INTEGRATE = $(physics)/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h
CLASS_ATOM_OUTPUT = $(physics)/include/class_defs/ATOM_OPERATIONS/class_atomoutput.h
CLASS_CP_GEN_WAVE  = $(physics)/include/class_defs/CP/class_gen_wave.h
CLASS_VX_SMPL  = $(physics)/include/class_defs/ATOM_OPERATIONS/class_vx_smpl.h
CLASS_CHARM_CONFIG = $(physics)/../include/configure.h

#-------------------------------------------------------
#               PROTO_TYPES

COORD_ENT     = $(physics)/interface/coords/proto_coords_entry.h
COORD_LOC     = $(physics)/interface/coords/proto_coords_local.h
COORD_CP_ENT  = $(physics)/interface/coords_cp/proto_coords_cp_entry.h
COORD_CP_LOC  = $(physics)/interface/coords_cp/proto_coords_cp_local.h
CPEWALD_ENT   = $(physics)/interface/cp_ewald/proto_cp_ewald_entry.h
CPEWALD_LOC   = $(physics)/interface/cp_ewald/proto_cp_ewald_local.h
HANDLE_ENT    = $(physics)/interface/handle/proto_handle_entry.h
HANDLE_LOC    = $(physics)/interface/handle/proto_handle_local.h
INTER_ENT     = $(physics)/interface/inter_params/proto_inter_params_entry.h
INTER_LOC     = $(physics)/interface/inter_params/proto_inter_params_local.h
INTRA_LOC     = $(physics)/interface/intra_params/proto_intra_params_local.h
INTRA_ENT     = $(physics)/interface/intra_params/proto_intra_params_entry.h
LISTS_ENT     = $(physics)/interface/lists/proto_lists_entry.h
LISTS_LOC     = $(physics)/interface/lists/proto_lists_local.h
MOL_ENT       = $(physics)/interface/mol_params/proto_mol_params_entry.h 
MOL_LOC       = $(physics)/interface/mol_params/proto_mol_params_local.h
SEARCH_ENT    = $(physics)/interface/search_base/proto_search_entry.h
SEARCH_LOC    = $(physics)/interface/search_base/proto_search_local.h
SIM_LOC       = $(physics)/interface/sim_params/proto_sim_params_local.h
SIM_ENT       = $(physics)/interface/sim_params/proto_sim_params_entry.h
SURF_PRMS_ENT = $(physics)/interface/surf_params/proto_surf_params_entry.h
SURF_PRMS_LOC = $(physics)/interface/surf_params/proto_surf_params_local.h
VPS_ENT       = $(physics)/interface/vps_params/proto_vps_params_entry.h
VPS_LOC       = $(physics)/interface/vps_params/proto_vps_params_local.h
WEIGH_NODE    = $(physics)/interface/lists/weights_nodes_128.h
PATH_INIT_ENT = $(physics)/interface/path_integral/proto_path_init_entry.h
PARSE_LOC     = $(physics)/interface/parse/proto_parse_local.h
PARSE_ENT     = $(physics)/interface/parse/proto_parse_entry.h
MATH          = $(physics)/mathlib/proto_math.h
FRIEND_ENT    = $(physics)/friend_lib/proto_friend_lib_entry.h

#-------------------------------------------------------
