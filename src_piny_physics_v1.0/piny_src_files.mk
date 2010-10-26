# Source files required to build the PINY (physics) portion of OpenAtom

FRIEND_FILES          = piny_malloc.C piny_pup.C friend_lib.C

INTERFACE_FILES       = parse.C \
                        interface_hand.C search_base_class.C \
                        data_base_handle.C \
                        control_sim_params.C set_sim_dict.C set_sim_params.C \
                        set_atm_NHC.C read_hmat.C \
                        control_surf_params.C set_surf_dict.C \
                        control_inter_params.C spline_fit.C \
                        set_inter_dict.C get_clong.C set_exclude.C \
                        path_integral_init.C exl_sort.C read_coord.C 

INTERFACE_CP_FILES    = set_wave_params.C  \
                        control_set_cp_ewald.C set_cp_ewald.C set_perd_corrs.C \
                        search_base_cp.C \
                        set_vps_dict.C control_vps_params.C \
                        weight_node_gauss_hermite.C gen_wave.C

INTERFACE_INTRA_FILES = close_intra_params.C control_intra_params.C \
                        control_res_params.C fetch_residue.C \
                        fetch_resbond_prm.C fetch_free_energy_index.C \
                        fetch_freeze.C init_intra_params.C \
                        manipulate_res_bonds.C replicate_mol.C residue_bond.C \
                        set_atm_mask.C set_atm_morph.C set_atm_params.C \
                        set_bend_bnd_params.C set_bend_params.C \
                        set_bond_params.C set_intra_dict.C \
                        set_intra_dict_pot.C set_intra_potent.C \
                        set_mol_name_params.C set_onfo_params.C \
                        set_res_bond_params.C set_res_def_params.C \
                        set_res_name_params.C set_res_morph_params.C \
                        set_grp_con_params.C set_tors_params.C intra_coefs.C \
                        fetch_hydrog_mass.C

INTERFACE_MOL_FILES   = control_mol_params.C control_set_mol_params.C \
                        set_base_file_params.C set_free_params.C \
                        set_mol_dict.C set_mol_params.C \
                        set_surf_params.C

# If building standalone piny?
MAIN_FILES            = PhysicsAtomPosInit.C vx_smpl.C Parainfoinit.C configure.C
# If building as part of OpenAtom
LIB_FILES             = PhysicsAtomPosInit.C vx_smpl.C Interface_ctrl.C \
                        Parainfoinit.C configure.C

MATH_FILES            = fft_generic.f mathlib.C

ABINITIO_PHYSICS_FILES = cp_eke.C cp_nl_energy_forc.C cp_nlmat.C \
                         cp_hart_ext.C cp_struct_fact.C \
                         cp_process_grad.C cp_pz_lda.C cp_min_std.C \
                         cp_min_CG.C cp_control_integrate.C  \
                         cp_lowdin.C cp_rspace_ion.C \
                         cp_grad_rho_ctrl.C cp_becke.C gen_wave.C \
                         cp_dynamics.C cp_vel_sampl.C cp_isokin.C \
                         cp_ees_nonlocal.C cp_pw_lda.C cp_pbe.C

CLASSICAL_PHYSICS_FILES = control_integrate.C write_output.C write_gen_header.C \
                          integration_drivers.C

#===========================================================================
#----------------- Sources needed to build the physics part of OpenAtom ---------------------------
piny_src_files = $(FRIEND_FILES) $(INTERFACE_FILES) \
       $(INTERFACE_CP_FILES) $(INTERFACE_INTRA_FILES) $(INTERFACE_MOL_FILES) \
       $(MATH_FILES) $(LIB_FILES) $(SPEC_FILES) \
       $(ABINITIO_PHYSICS_FILES) $(CLASSICAL_PHYSICS_FILES)

#----------------- Fragile sources that need to be compiled carefully ---------------------------
libphysics_fragile_src = \
              set_wave_params.C set_coef_NHC.C read_coef.C mall_properties.C \
              mall_coef.C control_set_cp_ewald.C set_perd_corrs.C \
              set_cp_ewald.C search_base_cp.C proj_vel_cp.C set_vps_dict.C \
              samp_vel_cp.C control_vps_params.C weight_node_gauss_hermite.C \
              control_vc_smpl.C control_vcnhc_smpl.C control_scale_cp.C  \
              gen_wave.C parse.C interface_hand.C search_base_class.C  \
              data_base_handle.C control_sim_params.C set_sim_dict.C  \
              set_sim_params.C set_atm_NHC.C read_coord.C molecule_decomp.C  \
              read_hmat.C mall_coord.C control_surf_params.C set_surf_dict.C \
              control_inter_params.C set_inter_dict.C get_clong.C \
              control_vnhc_smpl.C control_vx_smpl.C control_scale_class.C \
              proj_vel_class.C samp_vel_class.C set_exclude.C exl_sort.C \
              path_integral_init.C close_intra_params.C \
              control_intra_params.C control_res_params.C fetch_residue.C \
              fetch_resbond_prm.C fetch_free_energy_index.C fetch_freeze.C \
              init_intra_params.C manipulate_res_bonds.C replicate_mol.C \
              residue_bond.C set_atm_mask.C set_atm_morph.C set_atm_params.C \
              set_bend_bnd_params.C set_bend_params.C set_bond_params.C \
              set_intra_dict.C set_intra_dict_pot.C set_intra_potent.C \
              intra_coefs.C set_mol_name_params.C set_onfo_params.C \
              set_res_bond_params.C set_res_def_params.C \
              set_res_name_params.C set_res_morph_params.C \
              set_grp_con_params.C set_tors_params.C fetch_hydrog_mass.C \
              control_mol_params.C control_set_mol_params.C \
              set_base_file_params.C set_surf_params.C set_free_params.C \
              set_mol_dict.C set_mol_params.C configure.C

