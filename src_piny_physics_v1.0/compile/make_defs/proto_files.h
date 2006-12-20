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

STANDARD      = $(INCLUDES)/standard_include.h
DEFINES       = $(CODE)/include/class_defs/piny_constants.h 

#-------------------------------------------------------
#              TYPEDEFS 

CLASS_GEN     = $(CODE)/include/class_defs/allclass_gen.h
CLASS_MDINT   = $(CODE)/include/class_defs/allclass_mdintegrate.h
CLASS_MDATM   = $(CODE)/include/class_defs/allclass_mdatoms.h
CLASS_MDINTER = $(CODE)/include/class_defs/allclass_mdinter.h
CLASS_MDINTRA = $(CODE)/include/class_defs/allclass_mdintra.h
CLASS_CP      = $(CODE)/include/class_defs/allclass_cp.h
TYP_STAT      = $(CODE)/include/class_defs/typedefs_stat.h
TYP_PAR       = $(CODE)/include/class_defs/typedefs_par.h
DBG_FLAGS     = $(CODE)/../include/debug_flags.h

CLASS_PhysicsAtomPosinit = $(CODE)/include/class_defs/PINY_INIT/PhysicsAtomPosInit.h
CLASS_PARAINIT  = $(CODE)/include/class_defs/PINY_INIT/PhysicsParamTrans.h
CLASS_CP_NONLOC = $(CODE)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h
CLASS_CP_LOC    = $(CODE)/include/class_defs/CP_OPERATIONS/class_cplocal.h
CLASS_CP_XCFNCTLS = $(CODE)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
CLASS_CP_ORTHOG   = $(CODE)/include/class_defs/CP_OPERATIONS/class_cporthog.h
CLASS_CP_INTEGRATE = $(CODE)/include/class_defs/CP_OPERATIONS/class_cpintegrate.h
CLASS_CP_RSPACEION = $(CODE)/include/class_defs/CP_OPERATIONS/class_cprspaceion.h
CLASS_CHARM_ATOM = $(CODE)/../include/Atoms.h
CLASS_ATOM_INTEGRATE = $(CODE)/include/class_defs/ATOM_OPERATIONS/class_atomintegrate.h
CLASS_ATOM_OUTPUT = $(CODE)/include/class_defs/ATOM_OPERATIONS/class_atomoutput.h
CLASS_CP_GEN_WAVE  = $(CODE)/include/class_defs/CP/class_gen_wave.h
CLASS_VX_SMPL  = $(CODE)/include/class_defs/ATOM_OPERATIONS/class_vx_smpl.h
CLASS_CHARM_CONFIG = $(CODE)/../include/configure.h

#-------------------------------------------------------
#               PROTO_TYPES

COORD_ENT     = $(CODE)/interface/coords/proto_coords_entry.h
COORD_LOC     = $(CODE)/interface/coords/proto_coords_local.h
COORD_CP_ENT  = $(CODE)/interface/coords_cp/proto_coords_cp_entry.h
COORD_CP_LOC  = $(CODE)/interface/coords_cp/proto_coords_cp_local.h
CPEWALD_ENT   = $(CODE)/interface/cp_ewald/proto_cp_ewald_entry.h
CPEWALD_LOC   = $(CODE)/interface/cp_ewald/proto_cp_ewald_local.h
HANDLE_ENT    = $(CODE)/interface/handle/proto_handle_entry.h
HANDLE_LOC    = $(CODE)/interface/handle/proto_handle_local.h
INTER_ENT     = $(CODE)/interface/inter_params/proto_inter_params_entry.h
INTER_LOC     = $(CODE)/interface/inter_params/proto_inter_params_local.h
INTRA_LOC     = $(CODE)/interface/intra_params/proto_intra_params_local.h
INTRA_ENT     = $(CODE)/interface/intra_params/proto_intra_params_entry.h
LISTS_ENT     = $(CODE)/interface/lists/proto_lists_entry.h
LISTS_LOC     = $(CODE)/interface/lists/proto_lists_local.h
MOL_ENT       = $(CODE)/interface/mol_params/proto_mol_params_entry.h 
MOL_LOC       = $(CODE)/interface/mol_params/proto_mol_params_local.h
SEARCH_ENT    = $(CODE)/interface/search_base/proto_search_entry.h
SEARCH_LOC    = $(CODE)/interface/search_base/proto_search_local.h
SIM_LOC       = $(CODE)/interface/sim_params/proto_sim_params_local.h
SIM_ENT       = $(CODE)/interface/sim_params/proto_sim_params_entry.h
SURF_PRMS_ENT = $(CODE)/interface/surf_params/proto_surf_params_entry.h
SURF_PRMS_LOC = $(CODE)/interface/surf_params/proto_surf_params_local.h
VPS_ENT       = $(CODE)/interface/vps_params/proto_vps_params_entry.h
VPS_LOC       = $(CODE)/interface/vps_params/proto_vps_params_local.h
WEIGH_NODE    = $(CODE)/interface/lists/weights_nodes_128.h
PATH_INIT_ENT = $(CODE)/interface/path_integral/proto_path_init_entry.h
PARSE_LOC     = $(CODE)/interface/parse/proto_parse_local.h
PARSE_ENT     = $(CODE)/interface/parse/proto_parse_entry.h
MATH          = $(CODE)/mathlib/proto_math.h
FRIEND_ENT    = $(CODE)/friend_lib/proto_friend_lib_entry.h

#-------------------------------------------------------
