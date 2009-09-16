#===================================================================
#         AB INITO PHYSICS  
#===================================================================


#=================================================================
cp_eke.o    :  $(STANDARD) $(FFTW) \
               $(COMPLEX)  $(DBG_FLAGS) \
               $(CLASS_CP_NONLOC) \
               $(physics)/abinito_physics/cp_nonlocal/cp_eke.C
	$(COBJ) $(physics)/abinito_physics/cp_nonlocal/cp_eke.C
#------------------------------------------------------------------
cp_nl_energy_forc.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) \
               $(physics)/abinito_physics/cp_nonlocal/cp_nl_energy_forc.C
	$(COBJ) $(physics)/abinito_physics/cp_nonlocal/cp_nl_energy_forc.C

#------------------------------------------------------------------
cp_nlmat.o  :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) \
               $(physics)/abinito_physics/cp_nonlocal/cp_nlmat.C
	$(COBJ) $(physics)/abinito_physics/cp_nonlocal/cp_nlmat.C
#------------------------------------------------------------------
cp_ees_nonlocal.o  :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) \
               $(physics)/abinito_physics/cp_nonlocal/cp_ees_nonlocal.C
	$(COBJ) $(physics)/abinito_physics/cp_nonlocal/cp_ees_nonlocal.C
#------------------------------------------------------------------
cp_hart_ext.o :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) \
               $(CLASS_CP) $(FFTW)   $(DBG_FLAGS) \
               $(COMPLEX)   \
               $(CLASS_CP_LOC) \
               $(physics)/abinito_physics/cp_local/cp_hart_ext.C
	$(COBJ) $(physics)/abinito_physics/cp_local/cp_hart_ext.C

#------------------------------------------------------------------
cp_struct_fact.o :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
               $(CLASS_MDINTER) $(CLASS_MDINTRA) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_LOC) $(CLASS_CP_NONLOC)\
               $(physics)/abinito_physics/cp_nonlocal/cp_struct_fact.C
	$(COBJ) $(physics)/abinito_physics/cp_nonlocal/cp_struct_fact.C

#------------------------------------------------------------------
cp_process_grad.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(physics)/abinito_physics/cp_xc_fnctl/cp_process_grad.C
	$(COBJ) $(physics)/abinito_physics/cp_xc_fnctl/cp_process_grad.C

#------------------------------------------------------------------
cp_grad_rho_ctrl.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(physics)/abinito_physics/cp_xc_fnctl/cp_grad_rho_ctrl.C
	$(COBJ) $(physics)/abinito_physics/cp_xc_fnctl/cp_grad_rho_ctrl.C

#------------------------------------------------------------------
cp_becke.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(physics)/abinito_physics/cp_xc_fnctl/cp_becke.C
	$(COBJ) $(physics)/abinito_physics/cp_xc_fnctl/cp_becke.C

#------------------------------------------------------------------
cp_pz_lda.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(physics)/abinito_physics/cp_xc_fnctl/cp_pz_lda.C
	$(COBJ) $(physics)/abinito_physics/cp_xc_fnctl/cp_pz_lda.C

#------------------------------------------------------------------
cp_control_integrate.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE) \
               $(physics)/abinito_physics/cp_integrate/cp_control_integrate.C
	$(COBJ) $(physics)/abinito_physics/cp_integrate/cp_control_integrate.C
#------------------------------------------------------------------
cp_dynamics.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE) \
               $(physics)/abinito_physics/cp_integrate/cp_dynamics.C
	$(COBJ) $(physics)/abinito_physics/cp_integrate/cp_dynamics.C
#------------------------------------------------------------------
cp_vel_sampl.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE) \
               $(physics)/abinito_physics/cp_integrate/cp_vel_sampl.C
	$(COBJ) $(physics)/abinito_physics/cp_integrate/cp_vel_sampl.C
#------------------------------------------------------------------
cp_isokin.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE) \
               $(physics)/abinito_physics/cp_integrate/cp_isokin.C
	$(COBJ) $(physics)/abinito_physics/cp_integrate/cp_isokin.C
#------------------------------------------------------------------
cp_rspace_ion.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_RSPACEION) \
               $(physics)/abinito_physics/cp_ions/cp_rspace_ion.C
	$(COBJ) $(physics)/abinito_physics/cp_ions/cp_rspace_ion.C
#------------------------------------------------------------------
cp_min_std.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)  $(DBG_FLAGS) \
               $(CLASS_CP_INTEGRATE) \
               $(physics)/abinito_physics/cp_integrate/cp_min_std.C
	$(COBJ) $(physics)/abinito_physics/cp_integrate/cp_min_std.C

#------------------------------------------------------------------
cp_min_CG.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)  \
               $(CLASS_CP_INTEGRATE) \
               $(physics)/abinito_physics/cp_integrate/cp_min_CG.C
	$(COBJ) $(physics)/abinito_physics/cp_integrate/cp_min_CG.C

#------------------------------------------------------------------
cp_lowdin.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_ORTHOG) $(MATH) \
               $(physics)/abinito_physics/cp_orthog/cp_lowdin.C
	$(COBJ) $(physics)/abinito_physics/cp_orthog/cp_lowdin.C

#===================================================================
#===================================================================
#         CLASSICAL PHYSICS  
#===================================================================

#------------------------------------------------------------------
control_integrate.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM)   $(DBG_FLAGS) \
               $(CLASS_ATOM_INTEGRATE) \
               $(physics)/classical_physics/integrate/control_integrate.C
	$(COBJ) $(physics)/classical_physics/integrate/control_integrate.C
#------------------------------------------------------------------
integration_drivers.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_INTEGRATE)   $(DBG_FLAGS) \
               $(physics)/classical_physics/integrate/integration_drivers.C
	$(COBJ) $(physics)/classical_physics/integrate/integration_drivers.C
#------------------------------------------------------------------
write_gen_header.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_OUTPUT)   $(DBG_FLAGS) \
               $(physics)/classical_physics/output/write_gen_header.C
	$(COBJ) $(physics)/classical_physics/output/write_gen_header.C
#------------------------------------------------------------------
write_output.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_OUTPUT)   $(DBG_FLAGS) \
               $(physics)/classical_physics/output/write_output.C
	$(COBJ) $(physics)/classical_physics/output/write_output.C
#------------------------------------------------------------------

#===================================================================
#=================================================================
#               FRIEND_FILES
#=================================================================


#=================================================================
piny_malloc.o     :      $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(physics)/friend_lib/piny_malloc.C
	$(COBJ) $(DUAL_FFTW) -I$(FFT_HOME)/include $(physics)/friend_lib/piny_malloc.C

#-----------------------------------------------------------------
piny_pup.o     :          $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(physics)/friend_lib/piny_pup.C
	$(COBJ) $(physics)/friend_lib/piny_pup.C

#-----------------------------------------------------------------
friend_lib.o   :         $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(physics)/friend_lib/friend_lib.C
	$(COBJ) $(physics)/friend_lib/friend_lib.C

#=================================================================
#==================================================================
#           INTERFACE_CP_FILES
#==================================================================


#==================================================================
MOL_PARMS1 = $(physics)/interface/mol_params
DMOL_PARMS1 = $(physics)/interface/mol_params
#==================================================================



#==================================================================
set_wave_params.o   :    $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS1)/set_params/set_wave_params.C
	$(COBJ_CARE) $(MOL_PARMS1)/set_params/set_wave_params.C

#------------------------------------------------------------------
set_coef_NHC.o    :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(physics)/interface/coords_cp/set_coef_NHC.C
	$(COBJ_CARE) $(physics)/interface/coords_cp/set_coef_NHC.C

#------------------------------------------------------------------
read_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(ENR_CPCON_ENT) $(COORD_CP_ENT) $(HANDLE_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) \
                         $(physics)/interface/coords_cp/read_coef.C
	$(COBJ_CARE) $(physics)/interface/coords_cp/read_coef.C

#------------------------------------------------------------------
mall_properties.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) \
                         $(FRND_ENT) \
                         $(physics)/interface/coords_cp/mall_properties.C
	$(COBJ_CARE) $(physics)/interface/coords_cp/mall_properties.C

#------------------------------------------------------------------
mall_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) $(FRND_ENT) \
                         $(physics)/interface/coords_cp/mall_coef.C
	$(COBJ_CARE) $(physics)/interface/coords_cp/mall_coef.C

#------------------------------------------------------------------
control_set_cp_ewald.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(CPEWALD_ENT) $(CPEWALD_LOC) $(ENR_CP_LOC) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(physics)/interface/cp_ewald/control_set_cp_ewald.C
	$(COBJ_CARE) $(physics)/interface/cp_ewald/control_set_cp_ewald.C

#------------------------------------------------------------------
set_perd_corrs.o : $(STANDARD) $(DEFINES) \
                   $(CLASS_GEN) $(CPEWALD_ENTRY) \
                   $(CPEWALD_CORR) $(MATH) $(FRND_ENT) \
                   $(physics)/interface/cp_ewald/set_perd_corrs.C
	$(COBJ_CARE) $(physics)/interface/cp_ewald/set_perd_corrs.C

#------------------------------------------------------------------
set_cp_ewald.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(CPEWALD_LOC) $(MATH) $(FRND_ENT) \
                         $(physics)/interface/cp_ewald/set_cp_ewald.C
	$(COBJ_CARE) $(physics)/interface/cp_ewald/set_cp_ewald.C

#------------------------------------------------------------------
search_base_cp.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(SEARCH_ENT) $(VPS_LOC) $(HANDLE_ENT) \
                         $(physics)/interface/search_base/search_base_cp.C
	$(COBJ_CARE) $(physics)/interface/search_base/search_base_cp.C

#------------------------------------------------------------------
proj_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) $(CLASS_MDINTRA) \
                         $(SMPL_CP_LOC) $(SMPL_CP_ENT) $(ENR_CPCON_ENT) \
                         $(physics)/interface/vel_sampl_cp/proj_vel_cp.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_cp/proj_vel_cp.C

#------------------------------------------------------------------
set_vps_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
                         $(physics)/interface/vps_params/set_vps_dict.C
	$(COBJ_CARE) $(physics)/interface/vps_params/set_vps_dict.C

#------------------------------------------------------------------
samp_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(SMPL_CP_LOC)  $(MATH) \
                         $(physics)/interface/vel_sampl_cp/samp_vel_cp.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_cp/samp_vel_cp.C

#------------------------------------------------------------------
control_vps_params.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(COMM_WRAP) $(MATH) \
                         $(physics)/interface/vps_params/control_vps_params.C
	$(COBJ_CARE) $(physics)/interface/vps_params/control_vps_params.C

#------------------------------------------------------------------
weight_node_gauss_hermite.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
                         $(physics)/interface/vps_params/weight_node_gauss_hermite.C
	$(COBJ_CARE) $(physics)/interface/vps_params/weight_node_gauss_hermite.C

#------------------------------------------------------------------
control_vc_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
                         $(physics)/interface/vel_sampl_cp/control_vc_smpl.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_cp/control_vc_smpl.C

#------------------------------------------------------------------
control_vcnhc_smpl.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
                         $(physics)/interface/vel_sampl_cp/control_vcnhc_smpl.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_cp/control_vcnhc_smpl.C

#------------------------------------------------------------------
control_scale_cp.o    :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP)  \
                         $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP) \
                         $(physics)/interface/vel_sampl_cp/control_scale_cp.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_cp/control_scale_cp.C

#------------------------------------------------------------------
gen_wave.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) \
                         $(CLASS_CP_GEN_WAVE)  \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(MATH) \
                         $(physics)/interface/coords_cp/gen_wave.C
	$(COBJ_CARE) $(physics)/interface/coords_cp/gen_wave.C

#==================================================================



#==================================================================
#             INTERFACE_FILES
#==================================================================


#==================================================================
parse.o     :            $(STANDARD) $(DEFINES) \
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
                         $(INT_CPMIN_ENT) \
                         $(physics)/interface/parse/parse.C
	$(COBJ_CARE) $(physics)/interface/parse/parse.C

#==================================================================
interface_hand.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(HANDLE_ENT) \
                        $(physics)/interface/handle/interface_hand.C
	$(COBJ_CARE) $(physics)/interface/handle/interface_hand.C

#------------------------------------------------------------------
search_base_class.o   :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(SEARCH_LOC) \
                         $(FRND_ENT) \
                         $(physics)/interface/search_base/search_base_class.C
	$(COBJ_CARE) $(physics)/interface/search_base/search_base_class.C

#------------------------------------------------------------------
data_base_handle.o    :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(SEARCH_ENT) $(INTER_ENT) $(INTER_LOC) \
                         $(INTRA_LOC)  $(HANDLE_ENT) $(FRND_ENT) \
                         $(physics)/interface/search_base/data_base_handle.C
	$(COBJ_CARE) $(physics)/interface/search_base/data_base_handle.C

#=========================================================================
control_sim_params.o :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_ENT) $(SIM_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(physics)/interface/sim_params/control_sim_params.C
	$(COBJ_CARE) $(physics)/interface/sim_params/control_sim_params.C

#------------------------------------------------------------------
set_sim_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA)  \
                         $(CLASS_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_LOC) $(FRND_ENT) \
                         $(physics)/interface/sim_params/set_sim_dict.C
	$(COBJ_CARE) $(physics)/interface/sim_params/set_sim_dict.C

#------------------------------------------------------------------
set_sim_params.o  :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) $(SIM_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) $(TYP_STAT)\
                         $(physics)/interface/sim_params/set_sim_params.C
	$(COBJ_CARE) $(physics)/interface/sim_params/set_sim_params.C

#=========================================================================
set_atm_NHC.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(physics)/interface/coords/set_atm_NHC.C
	$(COBJ_CARE) $(physics)/interface/coords/set_atm_NHC.C

#------------------------------------------------------------------
read_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP) \
                         $(physics)/interface/coords/read_coord.C
	$(COBJ_CARE) $(physics)/interface/coords/read_coord.C

#------------------------------------------------------------------
molecule_decomp.o  :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINT)\
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(COORD_LOC) $(FRND_ENT) $(COMM_WRAP) \
                         $(physics)/interface/coords/molecule_decomp.C
	$(COBJ_CARE) $(physics)/interface/coords/molecule_decomp.C

#------------------------------------------------------------------
read_hmat.o     :        $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(HANDLE_ENT) $(FRND_ENT) $(MATH) \
                         $(physics)/interface/coords/read_hmat.C
	$(COBJ_CARE) $(physics)/interface/coords/read_hmat.C

#------------------------------------------------------------------
mall_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_LOC) $(FRND_ENT) \
                         $(physics)/interface/coords/mall_coord.C
	$(COBJ_CARE) $(physics)/interface/coords/mall_coord.C

#=========================================================================
control_surf_params.o : $(STANDARD) $(DEFINES) \
                        $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                        $(TYP_PAR) $(CLASS_MDINTRA) \
                        $(CLASS_CP) \
                        $(SURF_PRMS_ENT) $(SURF_PRMS_LOC) \
                        $(INTRA_LOC) $(SEARCH_ENT) \
                        $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
                        $(physics)/interface/surf_params/control_surf_params.C
	$(COBJ_CARE) $(physics)/interface/surf_params/control_surf_params.C

#------------------------------------------------------------------
set_surf_dict.o   :     $(STANDARD) $(DEFINES) \
                        $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                        $(CLASS_MDINTER) \
                        $(SURF_PRMS_LOC) $(FRND_ENT) \
                        $(physics)/interface/surf_params/set_surf_dict.C
	$(COBJ_CARE) $(physics)/interface/surf_params/set_surf_dict.C

#=========================================================================
control_inter_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(INTER_ENT) $(INTER_LOC) $(INTRA_LOC) $(SEARCH_ENT) \
                         $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
                         $(physics)/interface/inter_params/control_inter_params.C
	$(COBJ_CARE) $(physics)/interface/inter_params/control_inter_params.C

#------------------------------------------------------------------
set_inter_dict.o   :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(physics)/interface/inter_params/set_inter_dict.C
	$(COBJ_CARE) $(physics)/interface/inter_params/set_inter_dict.C

#------------------------------------------------------------------
get_clong.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(physics)/interface/inter_params/get_clong.C
	$(COBJ_CARE) $(physics)/interface/inter_params/get_clong.C

#------------------------------------------------------------------
spline_fit.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(physics)/interface/inter_params/spline_fit.C
	$(COBJ) $(physics)/interface/inter_params/spline_fit.C

#==================================================================
control_vnhc_smpl.o   :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
                         $(physics)/interface/vel_sampl_class/control_vnhc_smpl.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_class/control_vnhc_smpl.C

#------------------------------------------------------------------
control_vx_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COORD_LOC) \
                         $(COMM_WRAP) \
                         $(physics)/interface/vel_sampl_class/control_vx_smpl.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_class/control_vx_smpl.C

#------------------------------------------------------------------
control_scale_class.o : $(STANDARD) $(DEFINES) \
                        $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                        $(CLASS_MDINTRA) $(CLASS_GEN) $(TYP_PAR) \
                        $(SMPL_CLASS_ENT) $(COMM_WRAP) \
                        $(physics)/interface/vel_sampl_class/control_scale_class.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_class/control_scale_class.C

#------------------------------------------------------------------
proj_vel_class.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
                         $(physics)/interface/vel_sampl_class/proj_vel_class.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_class/proj_vel_class.C

#------------------------------------------------------------------
samp_vel_class.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(MATH) \
                         $(physics)/interface/vel_sampl_class/samp_vel_class.C
	$(COBJ_CARE) $(physics)/interface/vel_sampl_class/samp_vel_class.C

#==================================================================
set_exclude.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) $(WEIGH_NODE) \
                         $(physics)/interface/lists/set_exclude.C
	$(COBJ_CARE) $(physics)/interface/lists/set_exclude.C

#------------------------------------------------------------------
exl_sort.o     :         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(LISTS_LOC) $(FRND_ENT) \
                         $(physics)/interface/lists/exl_sort.C
	$(COBJ_CARE) $(physics)/interface/lists/exl_sort.C

#------------------------------------------------------------------
path_integral_init.o     :         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(PATH_INIT_ENT) $(MATH) $(FRND_ENT) \
                         $(physics)/interface/path_integral/path_integral_init.C
	$(COBJ_CARE) $(physics)/interface/path_integral/path_integral_init.C

#==================================================================





#==================================================================
#               INTERFACE_INTRA_FILES
#==================================================================


#==================================================================
SET_PARM  = $(physics)/interface/intra_params/set_params
DSET_PARM = $(physics)/interface/intra_params/set_params
#==================================================================


#==================================================================
close_intra_params.o :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) $(CLASS_GEN) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/close_intra_params.C
	$(COBJ_CARE) $(physics)/interface/intra_params/close_intra_params.C

#------------------------------------------------------------------
control_intra_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) $(PIMD_LOC) \
                         $(MATH) \
                         $(physics)/interface/intra_params/control_intra_params.C
	$(COBJ_CARE) $(physics)/interface/intra_params/control_intra_params.C

#------------------------------------------------------------------
control_res_params.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/control_res_params.C
	$(COBJ_CARE) $(physics)/interface/intra_params/control_res_params.C

#------------------------------------------------------------------
fetch_residue.o     :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(physics)/interface/intra_params/fetch_residue.C
	$(COBJ_CARE) $(physics)/interface/intra_params/fetch_residue.C

#------------------------------------------------------------------
fetch_resbond_prm.o  :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) \
                         $(physics)/interface/intra_params/fetch_resbond_prm.C
	$(COBJ_CARE) $(physics)/interface/intra_params/fetch_resbond_prm.C

#------------------------------------------------------------------
fetch_free_energy_index.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) \
                       $(physics)/interface/intra_params/fetch_free_energy_index.C
	$(COBJ_CARE) $(physics)/interface/intra_params/fetch_free_energy_index.C

#------------------------------------------------------------------
fetch_freeze.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/fetch_freeze.C
	$(COBJ_CARE) $(physics)/interface/intra_params/fetch_freeze.C

#------------------------------------------------------------------
init_intra_params.o  :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/init_intra_params.C
	$(COBJ_CARE) $(physics)/interface/intra_params/init_intra_params.C

#------------------------------------------------------------------
manipulate_res_bonds.o : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/manipulate_res_bonds.C
	$(COBJ_CARE) $(physics)/interface/intra_params/manipulate_res_bonds.C

#------------------------------------------------------------------
replicate_mol.o     :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/replicate_mol.C
	$(COBJ_CARE) $(physics)/interface/intra_params/replicate_mol.C
#------------------------------------------------------------------
residue_bond.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_INTRA) $(CLASS_MDINT)\
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/residue_bond.C
	$(COBJ_CARE) $(physics)/interface/intra_params/residue_bond.C

#==================================================================
set_atm_mask.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_mask.C
	$(COBJ_CARE) $(SET_PARM)/set_atm_mask.C

#------------------------------------------------------------------
set_atm_morph.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_morph.C
	$(COBJ_CARE) $(SET_PARM)/set_atm_morph.C

#------------------------------------------------------------------
set_atm_params.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_params.C
	$(COBJ_CARE) $(SET_PARM)/set_atm_params.C

#------------------------------------------------------------------
set_bend_bnd_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bend_bnd_params.C
	$(COBJ_CARE) $(SET_PARM)/set_bend_bnd_params.C

#------------------------------------------------------------------
set_bend_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bend_params.C
	$(COBJ_CARE) $(SET_PARM)/set_bend_params.C

#------------------------------------------------------------------
set_bond_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bond_params.C
	$(COBJ_CARE) $(SET_PARM)/set_bond_params.C

#------------------------------------------------------------------
set_intra_dict.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(SET_PARM)/set_intra_dict.C
	$(COBJ_CARE) $(SET_PARM)/set_intra_dict.C

#------------------------------------------------------------------
set_intra_dict_pot.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(SET_PARM)/set_intra_dict_pot.C
	$(COBJ_CARE) $(SET_PARM)/set_intra_dict_pot.C

#------------------------------------------------------------------
set_intra_potent.o  :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(FRND_ENT) \
                         $(physics)/interface/intra_params/set_intra_potent.C
	$(COBJ_CARE) $(physics)/interface/intra_params/set_intra_potent.C
#------------------------------------------------------------------
intra_coefs.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(physics)/interface/intra_params/intra_coefs.C
	$(COBJ_CARE) $(physics)/interface/intra_params/intra_coefs.C

#------------------------------------------------------------------
set_mol_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_mol_name_params.C
	$(COBJ_CARE) $(SET_PARM)/set_mol_name_params.C

#------------------------------------------------------------------
set_onfo_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_onfo_params.C
	$(COBJ_CARE) $(SET_PARM)/set_onfo_params.C

#------------------------------------------------------------------
set_res_bond_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_bond_params.C
	$(COBJ_CARE) $(SET_PARM)/set_res_bond_params.C

#------------------------------------------------------------------
set_res_def_params.o   : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_def_params.C
	$(COBJ_CARE) $(SET_PARM)/set_res_def_params.C

#------------------------------------------------------------------
set_res_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_name_params.C
	$(COBJ_CARE) $(SET_PARM)/set_res_name_params.C

#------------------------------------------------------------------
set_res_morph_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_morph_params.C
	$(COBJ_CARE) $(SET_PARM)/set_res_morph_params.C

#------------------------------------------------------------------
set_grp_con_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_grp_con_params.C
	$(COBJ_CARE) $(SET_PARM)/set_grp_con_params.C

#------------------------------------------------------------------
set_tors_params.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(SET_PARM)/set_tors_params.C
	$(COBJ_CARE) $(SET_PARM)/set_tors_params.C

#------------------------------------------------------------------
fetch_hydrog_mass.o :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(physics)/interface/intra_params/fetch_hydrog_mass.C
	$(COBJ_CARE) $(physics)/interface/intra_params/fetch_hydrog_mass.C

#==================================================================
#==================================================================
#            INTERFACE_MOL_FILES
#==================================================================


#==================================================================
MOL_PARMS = $(physics)/interface/mol_params
DMOL_PARMS = $(physics)/interface/mol_params
#==================================================================


#==================================================================
control_mol_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(MOL_ENT) $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MOL_PARMS)/control_mol_params.C
	$(COBJ_CARE) $(MOL_PARMS)/control_mol_params.C

#----------------------------------------------------------------
control_set_mol_params.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MOL_PARMS)/set_params/control_set_mol_params.C
	$(COBJ_CARE) $(MOL_PARMS)/set_params/control_set_mol_params.C

#----------------------------------------------------------------
set_base_file_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(MOL_LOC) \
                         $(MOL_PARMS)/set_params/set_base_file_params.C
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_base_file_params.C

#----------------------------------------------------------------
set_surf_params.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_surf_params.C
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_surf_params.C

#----------------------------------------------------------------
set_free_params.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_free_params.C
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_free_params.C

#----------------------------------------------------------------
set_mol_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) $(CLASS_CP) \
                         $(MOL_LOC) $(FRND_ENT) \
                         $(MOL_PARMS)/set_params/set_mol_dict.C
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_mol_dict.C

#----------------------------------------------------------------
set_mol_params.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(MOL_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_mol_params.C
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_mol_params.C

#==================================================================

#=================================================================
#              MAIN_FILES 
#=================================================================


#==========================================================================
Interface_ctrl.o     :           $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) \
                         $(physics)/piny_to_driver/Interface_ctrl.C \
                         $(physics)/piny_to_driver/Interface_ctrl.h
	cp -f $(physics)/piny_to_driver/Interface_ctrl.h $(physics)/include/class_defs
	$(COBJ) $(physics)/piny_to_driver/Interface_ctrl.C

#==========================================================================
PhysicsAtomPosInit.o     : $(STANDARD) $(DEFINES) \
                           $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                           $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                           $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) \
                           $(CLASS_PhysicsAtomPosinit) \
                           $(physics)/piny_to_driver/PhysicsAtomPosInit.C 
	$(COBJ) $(physics)/piny_to_driver/PhysicsAtomPosInit.C
#==========================================================================
vx_smpl.o     :          $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) $(CLASS_VX_SMPL) \
                         $(physics)/classical_physics/vel_smpl_atms/vx_smpl.C 
	$(COBJ) $(physics)/classical_physics/vel_smpl_atms/vx_smpl.C
                                                                               
#==========================================================================
ParaInfoInit.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(DBG_FLAGS) \
                         $(CLASS_PARAINIT) \
                         $(physics)/piny_to_driver/Parainfoinit.C
	$(COBJ) $(physics)/piny_to_driver/Parainfoinit.C

#==========================================================================
configure.o     :      $(STANDARD) $(DEFINES) $(TYP_PAR) \
                       $(CLASS_CHARM_CONFIG) $(HANDLE_ENT) $(FRIEND_ENT) \
                       $(physics)/piny_to_driver/configure.C
	$(COBJ_CARE) $(physics)/piny_to_driver/configure.C

#==========================================================================


#===================================================================
#               MATH_FILES
#===================================================================


#=================================================================
mathlib.o     :          $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) \
                         $(MATH) $(COMM_WRAP) \
                         $(physics)/mathlib/mathlib.C
	$(COBJ) $(physics)/mathlib/mathlib.C

#------------------------------------------------------------------
fft_generic.o       :    $(DEFINES) $(STANDARD) $(MATH) \
                         $(physics)/mathlib/fft_generic.f
	$(FOBJ) $(physics)/mathlib/fft_generic.f

#===================================================================
