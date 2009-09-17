# The relevant source files for this project
include $(pinymake)/make_defs/obj_files.h
libphysics_src=
libphysics_obj= $(addsuffix .o, $(basename $(libphysics_src)) ) $(LIB_OBJS)
libphysics_intf=
CPPFLAGS += -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 

# Specify the list of directories whose contents should be stripped from prerequisite lists 
# during dependency generation
DEPSTRIPDIRS += 
# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
fileTypes     = $(sort $(suffix $(libphysics_src) $(libphysics_intf)) )
#$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(STANDARD_INC)) )
vpath standard_include.h $(STANDARD_INC)
vpath %.C $(pinysrcdirs)
vpath fft_generic.f $(physics)/mathlib

# The primary target for this module
$(libphysics): $(libphysics_obj) | $(LIB_DECLS)
	$(info-ar)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

#----------------- Temporary legacy variable definitions ---------------------------
include $(pinymake)/make_defs/proto_files.h
include $(pinymake)/targetdeps.mk

z3dfft_dec_noimsl.o : $(physics)/mathlib/z3dfft_dec_noimsl.f $(PROTO)

# All the source files that need to compiled more carefully 
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
libphysics_fragile_obj = $(addsuffix .o, $(basename $(libphysics_fragile_src)) )
# Change the optimization/debugging flags just for these source files
$(libphysics_fragile_obj): OPT = $(OPT_CARE)

#-----------------------------------------------------------
# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(libphysics_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(libphysics_intf:.ci=.di)
endif

