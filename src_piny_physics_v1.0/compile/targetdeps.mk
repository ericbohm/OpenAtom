#===================================================================
#         AB INITO PHYSICS  
#===================================================================


#=================================================================
cp_eke.o    : $(physics)/abinito_physics/cp_nonlocal/cp_eke.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)  $(DBG_FLAGS) \
               $(CLASS_CP_NONLOC) 
#------------------------------------------------------------------
cp_nl_energy_forc.o : $(physics)/abinito_physics/cp_nonlocal/cp_nl_energy_forc.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) 
#------------------------------------------------------------------
cp_nlmat.o  : $(physics)/abinito_physics/cp_nonlocal/cp_nlmat.C \
               $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) 
#------------------------------------------------------------------
cp_ees_nonlocal.o  : $(physics)/abinito_physics/cp_nonlocal/cp_ees_nonlocal.C \
               $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) 
#------------------------------------------------------------------
cp_hart_ext.o : $(physics)/abinito_physics/cp_local/cp_hart_ext.C \
               $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) \
               $(CLASS_CP) $(FFTW)   $(DBG_FLAGS) \
               $(COMPLEX)   \
               $(CLASS_CP_LOC) 
#------------------------------------------------------------------
cp_struct_fact.o : $(physics)/abinito_physics/cp_nonlocal/cp_struct_fact.C \
               $(STANDARD) $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
               $(CLASS_MDINTER) $(CLASS_MDINTRA) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_LOC) $(CLASS_CP_NONLOC)
#------------------------------------------------------------------
cp_process_grad.o : $(physics)/abinito_physics/cp_xc_fnctl/cp_process_grad.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) 
#------------------------------------------------------------------
cp_grad_rho_ctrl.o : $(physics)/abinito_physics/cp_xc_fnctl/cp_grad_rho_ctrl.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) 


#------------------------------------------------------------------
cp_becke.o : $(physics)/abinito_physics/cp_xc_fnctl/cp_becke.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS)


#------------------------------------------------------------------
cp_pz_lda.o : $(physics)/abinito_physics/cp_xc_fnctl/cp_pz_lda.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS)


#------------------------------------------------------------------
cp_control_integrate.o : $(physics)/abinito_physics/cp_integrate/cp_control_integrate.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE)

#------------------------------------------------------------------
cp_dynamics.o : $(physics)/abinito_physics/cp_integrate/cp_dynamics.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE)

#------------------------------------------------------------------
cp_vel_sampl.o : $(physics)/abinito_physics/cp_integrate/cp_vel_sampl.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE)

#------------------------------------------------------------------
cp_isokin.o : $(physics)/abinito_physics/cp_integrate/cp_isokin.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE)

#------------------------------------------------------------------
cp_rspace_ion.o : $(physics)/abinito_physics/cp_ions/cp_rspace_ion.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_RSPACEION)

#------------------------------------------------------------------
cp_min_std.o : $(physics)/abinito_physics/cp_integrate/cp_min_std.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)  $(DBG_FLAGS) \
               $(CLASS_CP_INTEGRATE)


#------------------------------------------------------------------
cp_min_CG.o : $(physics)/abinito_physics/cp_integrate/cp_min_CG.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)  \
               $(CLASS_CP_INTEGRATE)


#------------------------------------------------------------------
cp_lowdin.o : $(physics)/abinito_physics/cp_orthog/cp_lowdin.C \
               $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_ORTHOG) $(MATH)


#===================================================================
#===================================================================
#         CLASSICAL PHYSICS  
#===================================================================

#------------------------------------------------------------------
control_integrate.o : $(physics)/classical_physics/integrate/control_integrate.C \
               $(STANDARD) \
               $(CLASS_CHARM_ATOM)   $(DBG_FLAGS) \
               $(CLASS_ATOM_INTEGRATE)

#------------------------------------------------------------------
integration_drivers.o : $(physics)/classical_physics/integrate/integration_drivers.C \
               $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_INTEGRATE)   $(DBG_FLAGS)

#------------------------------------------------------------------
write_gen_header.o : $(physics)/classical_physics/output/write_gen_header.C \
               $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_OUTPUT)   $(DBG_FLAGS)

#------------------------------------------------------------------
write_output.o : $(physics)/classical_physics/output/write_output.C \
               $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_OUTPUT)   $(DBG_FLAGS)

#------------------------------------------------------------------

#===================================================================
#=================================================================
#               FRIEND_FILES
#=================================================================


#=================================================================
piny_malloc.o     : $(physics)/friend_lib/piny_malloc.C \
                   $(STANDARD) $(DEFINES) $(FRND_ENT)
#-----------------------------------------------------------------
piny_pup.o     : $(physics)/friend_lib/piny_pup.C \
                       $(STANDARD) $(DEFINES) \
                         $(FRND_ENT)


#-----------------------------------------------------------------
friend_lib.o   : $(physics)/friend_lib/friend_lib.C \
                      $(STANDARD) $(DEFINES) \
                         $(FRND_ENT)


#=================================================================
#==================================================================
#           INTERFACE_CP_FILES
#==================================================================


#==================================================================
MOL_PARMS1 = $(physics)/interface/mol_params
DMOL_PARMS1 = $(physics)/interface/mol_params
#==================================================================



#==================================================================
set_wave_params.o   : $(MOL_PARMS1)/set_params/set_wave_params.C \
                 $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_CP) \
                         $(MOL_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_coef_NHC.o    : $(physics)/interface/coords_cp/set_coef_NHC.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(FRND_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
read_coef.o     : $(physics)/interface/coords_cp/read_coef.C \
                     $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(ENR_CPCON_ENT) $(COORD_CP_ENT) $(HANDLE_ENT) \
                         $(FRND_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
mall_properties.o     : $(physics)/interface/coords_cp/mall_properties.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) \
                         $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
mall_coef.o     : $(physics)/interface/coords_cp/mall_coef.C \
                     $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_set_cp_ewald.o : $(physics)/interface/cp_ewald/control_set_cp_ewald.C \
              $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(CPEWALD_ENT) $(CPEWALD_LOC) $(ENR_CP_LOC) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_perd_corrs.o : $(physics)/interface/cp_ewald/set_perd_corrs.C \
              $(STANDARD) $(DEFINES) \
                   $(CLASS_GEN) $(CPEWALD_ENTRY) \
                   $(CPEWALD_CORR) $(MATH) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_cp_ewald.o     : $(physics)/interface/cp_ewald/set_cp_ewald.C \
                  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(CPEWALD_LOC) $(MATH) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
search_base_cp.o  : $(physics)/interface/search_base/search_base_cp.C \
                   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(SEARCH_ENT) $(VPS_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
proj_vel_cp.o     : $(physics)/interface/vel_sampl_cp/proj_vel_cp.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) $(CLASS_MDINTRA) \
                         $(SMPL_CP_LOC) $(SMPL_CP_ENT) $(ENR_CPCON_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_vps_dict.o     : $(physics)/interface/vps_params/set_vps_dict.C \
                  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
samp_vel_cp.o     : $(physics)/interface/vel_sampl_cp/samp_vel_cp.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(SMPL_CP_LOC)  $(MATH)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_vps_params.o  : $(physics)/interface/vps_params/control_vps_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(COMM_WRAP) $(MATH)
	$(COBJ_CARE)

#------------------------------------------------------------------
weight_node_gauss_hermite.o     : $(physics)/interface/vps_params/weight_node_gauss_hermite.C \
                  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_vc_smpl.o     : $(physics)/interface/vel_sampl_cp/control_vc_smpl.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_vcnhc_smpl.o  : $(physics)/interface/vel_sampl_cp/control_vcnhc_smpl.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_scale_cp.o    : $(physics)/interface/vel_sampl_cp/control_scale_cp.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP)  \
                         $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
gen_wave.o  : $(physics)/interface/coords_cp/gen_wave.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) \
                         $(CLASS_CP_GEN_WAVE)  \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(MATH)
	$(COBJ_CARE)

#==================================================================



#==================================================================
#             INTERFACE_FILES
#==================================================================


#==================================================================
parse.o     : $(physics)/interface/parse/parse.C \
                         $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(TYP_PAR) $(TYP_STAT) $(PARSE_ENT) $(PARSE_LOC) \
                         $(SIM_ENT) $(MOL_ENT) $(INTRA_ENT) $(COORD_ENT) \
                         $(COORD_LOC) $(CPEWALD_ENT) $(INTER_ENT) \
                         $(INTER_LOC) $(VPS_ENT) $(LISTS_ENT) $(SCRATCH_ENT) \
                         $(ENR_CPCON_ENT) $(REAL_LOC) $(SAMPL_CLASS_ENT) \
                         $(SAMPL_CP_ENT) $(SAMPL_CP_LOC) $(COORD_CP_ENT) \
                         $(COORD_CP_LOC) $(MATH) $(PATH_INIT_ENT) \
                         $(FRND_ENT) $(COMM_ENT) $(COMM_LOC) $(COMM_WRAP) \
                         $(INT_CPMIN_ENT)
	$(COBJ_CARE)

#==================================================================
interface_hand.o     : $(physics)/interface/handle/interface_hand.C \
                $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
search_base_class.o   : $(physics)/interface/search_base/search_base_class.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(SEARCH_LOC) \
                         $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
data_base_handle.o    : $(physics)/interface/search_base/data_base_handle.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(SEARCH_ENT) $(INTER_ENT) $(INTER_LOC) \
                         $(INTRA_LOC)  $(HANDLE_ENT) $(FRND_ENT)
	$(COBJ_CARE)

#=========================================================================
control_sim_params.o : $(physics)/interface/sim_params/control_sim_params.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_ENT) $(SIM_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_sim_dict.o     : $(physics)/interface/sim_params/set_sim_dict.C \
                  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA)  \
                         $(CLASS_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_sim_params.o  : $(physics)/interface/sim_params/set_sim_params.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) $(SIM_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) $(TYP_STAT)
	$(COBJ_CARE)

#=========================================================================
set_atm_NHC.o     : $(physics)/interface/coords/set_atm_NHC.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(FRND_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
read_coord.o     : $(physics)/interface/coords/read_coord.C \
                    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
molecule_decomp.o  : $(physics)/interface/coords/molecule_decomp.C \
                  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINT)\
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(COORD_LOC) $(FRND_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
read_hmat.o     : $(physics)/interface/coords/read_hmat.C \
                     $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(HANDLE_ENT) $(FRND_ENT) $(MATH)
	$(COBJ_CARE)

#------------------------------------------------------------------
mall_coord.o     : $(physics)/interface/coords/mall_coord.C \
                    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#=========================================================================
control_surf_params.o : $(physics)/interface/surf_params/control_surf_params.C \
              $(STANDARD) $(DEFINES) \
                        $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                        $(TYP_PAR) $(CLASS_MDINTRA) \
                        $(CLASS_CP) \
                        $(SURF_PRMS_ENT) $(SURF_PRMS_LOC) \
                        $(INTRA_LOC) $(SEARCH_ENT) \
                        $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_surf_dict.o   : $(physics)/interface/surf_params/set_surf_dict.C \
                  $(STANDARD) $(DEFINES) \
                        $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                        $(CLASS_MDINTER) \
                        $(SURF_PRMS_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#=========================================================================
control_inter_params.o : $(physics)/interface/inter_params/control_inter_params.C \
              $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(INTER_ENT) $(INTER_LOC) $(INTRA_LOC) $(SEARCH_ENT) \
                         $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_inter_dict.o   : $(physics)/interface/inter_params/set_inter_dict.C \
                  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
get_clong.o     : $(physics)/interface/inter_params/get_clong.C \
                     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
spline_fit.o     : $(physics)/interface/inter_params/spline_fit.C \
                    $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT)


#==================================================================
control_vnhc_smpl.o   : $(physics)/interface/vel_sampl_class/control_vnhc_smpl.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_vx_smpl.o     : $(physics)/interface/vel_sampl_class/control_vx_smpl.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COORD_LOC) \
                         $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_scale_class.o : $(physics)/interface/vel_sampl_class/control_scale_class.C \
              $(STANDARD) $(DEFINES) \
                        $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                        $(CLASS_MDINTRA) $(CLASS_GEN) $(TYP_PAR) \
                        $(SMPL_CLASS_ENT) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
proj_vel_class.o     : $(physics)/interface/vel_sampl_class/proj_vel_class.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP)
	$(COBJ_CARE)

#------------------------------------------------------------------
samp_vel_class.o     : $(physics)/interface/vel_sampl_class/samp_vel_class.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(MATH)
	$(COBJ_CARE)

#==================================================================
set_exclude.o     : $(physics)/interface/lists/set_exclude.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) $(WEIGH_NODE)
	$(COBJ_CARE)

#------------------------------------------------------------------
exl_sort.o     : $(physics)/interface/lists/exl_sort.C \
                      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(LISTS_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
path_integral_init.o     : $(physics)/interface/path_integral/path_integral_init.C \
                      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(PATH_INIT_ENT) $(MATH) $(FRND_ENT)
	$(COBJ_CARE)

#==================================================================





#==================================================================
#               INTERFACE_INTRA_FILES
#==================================================================


#==================================================================
SET_PARM  = $(physics)/interface/intra_params/set_params
DSET_PARM = $(physics)/interface/intra_params/set_params
#==================================================================


#==================================================================
close_intra_params.o : $(physics)/interface/intra_params/close_intra_params.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) $(CLASS_GEN) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_intra_params.o : $(physics)/interface/intra_params/control_intra_params.C \
              $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) $(PIMD_LOC) \
                         $(MATH)
	$(COBJ_CARE)

#------------------------------------------------------------------
control_res_params.o  : $(physics)/interface/intra_params/control_res_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
fetch_residue.o     : $(physics)/interface/intra_params/fetch_residue.C \
                 $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
fetch_resbond_prm.o  : $(physics)/interface/intra_params/fetch_resbond_prm.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC)
	$(COBJ_CARE)

#------------------------------------------------------------------
fetch_free_energy_index.o     : $(physics)/interface/intra_params/fetch_free_energy_index.C \
              \
                         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC)
	$(COBJ_CARE)

#------------------------------------------------------------------
fetch_freeze.o     : $(physics)/interface/intra_params/fetch_freeze.C \
                  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
init_intra_params.o  : $(physics)/interface/intra_params/init_intra_params.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
manipulate_res_bonds.o : $(physics)/interface/intra_params/manipulate_res_bonds.C \
              $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
replicate_mol.o     : $(physics)/interface/intra_params/replicate_mol.C \
                 $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)
#------------------------------------------------------------------
residue_bond.o     : $(physics)/interface/intra_params/residue_bond.C \
                  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_INTRA) $(CLASS_MDINT)\
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#==================================================================
set_atm_mask.o     : $(SET_PARM)/set_atm_mask.C \
                  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_atm_morph.o     : $(SET_PARM)/set_atm_morph.C \
                 $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_atm_params.o     : $(SET_PARM)/set_atm_params.C \
                $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_bend_bnd_params.o  : $(SET_PARM)/set_bend_bnd_params.C \
              $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_bend_params.o     : $(SET_PARM)/set_bend_params.C \
               $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_bond_params.o     : $(SET_PARM)/set_bond_params.C \
               $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_intra_dict.o     : $(SET_PARM)/set_intra_dict.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_intra_dict_pot.o  : $(SET_PARM)/set_intra_dict_pot.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_intra_potent.o  : $(physics)/interface/intra_params/set_intra_potent.C \
                 $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(FRND_ENT)
	$(COBJ_CARE)
#------------------------------------------------------------------
intra_coefs.o     : $(physics)/interface/intra_params/intra_coefs.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_mol_name_params.o : $(SET_PARM)/set_mol_name_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_onfo_params.o     : $(SET_PARM)/set_onfo_params.C \
               $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_res_bond_params.o  : $(SET_PARM)/set_res_bond_params.C \
              $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_res_def_params.o   : $(SET_PARM)/set_res_def_params.C \
              $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_res_name_params.o : $(SET_PARM)/set_res_name_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_res_morph_params.o : $(SET_PARM)/set_res_morph_params.C \
              $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_grp_con_params.o  : $(SET_PARM)/set_grp_con_params.C \
               $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
set_tors_params.o     : $(SET_PARM)/set_tors_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT)
	$(COBJ_CARE)

#------------------------------------------------------------------
fetch_hydrog_mass.o : $(physics)/interface/intra_params/fetch_hydrog_mass.C \
                 $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#==================================================================
#==================================================================
#            INTERFACE_MOL_FILES
#==================================================================


#==================================================================
MOL_PARMS = $(physics)/interface/mol_params
DMOL_PARMS = $(physics)/interface/mol_params
#==================================================================


#==================================================================
control_mol_params.o  : $(MOL_PARMS)/control_mol_params.C \
               $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(MOL_ENT) $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT)
	$(COBJ_CARE)

#----------------------------------------------------------------
control_set_mol_params.o : $(MOL_PARMS)/set_params/control_set_mol_params.C \
              \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT)
	$(COBJ_CARE)

#----------------------------------------------------------------
set_base_file_params.o : $(MOL_PARMS)/set_params/set_base_file_params.C \
              $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(MOL_LOC)
	$(COBJ_CARE)

#----------------------------------------------------------------
set_surf_params.o     : $(MOL_PARMS)/set_params/set_surf_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(MOL_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#----------------------------------------------------------------
set_free_params.o     : $(MOL_PARMS)/set_params/set_free_params.C \
               $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(MOL_LOC) $(HANDLE_ENT)
	$(COBJ_CARE)

#----------------------------------------------------------------
set_mol_dict.o     : $(MOL_PARMS)/set_params/set_mol_dict.C \
                  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) $(CLASS_CP) \
                         $(MOL_LOC) $(FRND_ENT)
	$(COBJ_CARE)

#----------------------------------------------------------------
set_mol_params.o     : $(MOL_PARMS)/set_params/set_mol_params.C \
                $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(MOL_LOC) $(FRND_ENT) $(HANDLE_ENT)
	$(COBJ_CARE)

#==================================================================

#=================================================================
#              MAIN_FILES 
#=================================================================


#==========================================================================
Interface_ctrl.o     : $(physics)/piny_to_driver/Interface_ctrl.C \
                        $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) \
                         $(physics)/piny_to_driver/Interface_ctrl.h

#	cp -f $(physics)/piny_to_driver/Interface_ctrl.h $(physics)/include/class_defs


#==========================================================================
PhysicsAtomPosInit.o     : $(physics)/piny_to_driver/PhysicsAtomPosInit.C  \
              $(STANDARD) $(DEFINES) \
                           $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                           $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                           $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) \
                           $(CLASS_PhysicsAtomPosinit)

#==========================================================================
vx_smpl.o     : $(physics)/classical_physics/vel_smpl_atms/vx_smpl.C  \
                       $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) $(CLASS_VX_SMPL)

                                                                               
#==========================================================================
Parainfoinit.o     : $(physics)/piny_to_driver/Parainfoinit.C \
                   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(DBG_FLAGS) \
                         $(CLASS_PARAINIT)


#==========================================================================
configure.o     : $(physics)/piny_to_driver/configure.C \
                   $(STANDARD) $(DEFINES) $(TYP_PAR) \
                       $(CLASS_CHARM_CONFIG) $(HANDLE_ENT) $(FRIEND_ENT)
	$(COBJ_CARE)

#==========================================================================


#===================================================================
#               MATH_FILES
#===================================================================


#=================================================================
mathlib.o     : $(physics)/mathlib/mathlib.C \
                       $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) \
                         $(MATH) $(COMM_WRAP)


#------------------------------------------------------------------
fft_generic.o       : $(physics)/mathlib/fft_generic.f $(DEFINES) $(STANDARD) $(MATH)

#===================================================================
