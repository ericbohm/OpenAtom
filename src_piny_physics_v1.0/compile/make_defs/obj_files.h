#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================
#
#  Object and DECL files for PI_MD.
#
#===========================================================================
#ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#===========================================================================


FRIEND_FILES          = piny_malloc.o piny_pup.o friend_lib.o

INTERFACE_FILES       = parse.o \
                        interface_hand.o search_base_class.o \
                        data_base_handle.o \
                        control_sim_params.o set_sim_dict.o set_sim_params.o \
                        set_atm_NHC.o read_hmat.o \
                        control_surf_params.o set_surf_dict.o \
                        control_inter_params.o spline_fit.o \
                        set_inter_dict.o get_clong.o set_exclude.o \
                        path_integral_init.o exl_sort.o read_coord.o 

INTERFACE_CP_FILES    = set_wave_params.o set_coef_NHC.o \
                        control_set_cp_ewald.o set_cp_ewald.o \
                        search_base_cp.o \
                        set_vps_dict.o control_vps_params.o \
                        weight_node_gauss_hermite.o

INTERFACE_INTRA_FILES = close_intra_params.o control_intra_params.o \
                        control_res_params.o fetch_residue.o \
                        fetch_resbond_prm.o fetch_free_energy_index.o \
                        fetch_freeze.o init_intra_params.o \
                        manipulate_res_bonds.o replicate_mol.o residue_bond.o \
                        set_atm_mask.o set_atm_morph.o set_atm_params.o \
                        set_bend_bnd_params.o set_bend_params.o \
                        set_bond_params.o set_intra_dict.o \
                        set_intra_dict_pot.o set_intra_potent.o \
                        set_mol_name_params.o set_onfo_params.o \
                        set_res_bond_params.o set_res_def_params.o \
                        set_res_name_params.o set_res_morph_params.o \
                        set_grp_con_params.o set_tors_params.o intra_coefs.o \
                        fetch_hydrog_mass.o

INTERFACE_MOL_FILES   = control_mol_params.o control_set_mol_params.o \
                        set_base_file_params.o set_free_params.o \
                        set_mol_dict.o set_mol_params.o \
                        set_surf_params.o

MAIN_FILES            = PhysicsAtomPosInit.o vx_smpl.o ParaInfoInit.o configure.o

LIB_FILES             = PhysicsAtomPosInit.o vx_smpl.o Interface_ctrl.o \
                        ParaInfoInit.o configure.o

MATH_FILES            = fft_generic.o mathlib.o

MAIN_DECLS            = 
LIB_DECLS             = Interface_ctrl.decl.h

ABINITIO_PHYSICS_FILES = cp_eke.o cp_nl_energy_forc.o cp_nlmat.o \
                         cp_hart_ext.o cp_struct_fact.o \
                         cp_process_grad.o cp_pz_lda.o cp_min_std.o \
                         cp_min_CG.o cp_control_integrate.o  \
                         cp_lowdin.o cp_rspace_ion.o \
                         cp_grad_rho_ctrl.o cp_becke.o gen_wave.o \
                         cp_dynamics.o cp_vel_sampl.o cp_isokin.o \
                         cp_ees_nonlocal.o

CLASSICAL_PHYSICS_FILES = control_integrate.o write_output.o write_gen_header.o \
                          integration_drivers.o

#===========================================================================

DECLS = $(MAIN_DECLS) $(SPEC_DECLS)

OBJS = $(FRIEND_FILES) $(INTERFACE_FILES) \
       $(INTERFACE_CP_FILES) $(INTERFACE_INTRA_FILES) $(INTERFACE_MOL_FILES) \
       $(MATH_FILES) $(MAIN_FILES) $(SPEC_FILES) \
       $(ABINITIO_PHYSICS_FILES) $(CLASSICAL_PHYSICS_FILES)

OBJS_COMP = $(DECLS) $(OBJS) 

LIB_OBJS = $(FRIEND_FILES) $(INTERFACE_FILES) \
       $(INTERFACE_CP_FILES) $(INTERFACE_INTRA_FILES) $(INTERFACE_MOL_FILES) \
       $(MATH_FILES) $(LIB_FILES) $(SPEC_FILES) \
       $(ABINITIO_PHYSICS_FILES) $(CLASSICAL_PHYSICS_FILES)

LIB_OBJS_COMP = $(LIB_DECLS) $(LIB_OBJS)

#===========================================================================
