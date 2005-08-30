#==================================================================
#           INTERFACE_CP_FILES
#==================================================================


#==================================================================
MOL_PARMS1 = $(CODE)/interface/mol_params
DMOL_PARMS1 = $(DCODE)/interface/mol_params
#==================================================================



#==================================================================
set_wave_params.o   :    $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS1)/set_params/set_wave_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS1)/set_params/set_wave_params.C

#------------------------------------------------------------------
set_coef_NHC.o    :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords_cp/set_coef_NHC.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/set_coef_NHC.C

#------------------------------------------------------------------
read_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(ENR_CPCON_ENT) $(COORD_CP_ENT) $(HANDLE_ENT) \
                         $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords_cp/read_coef.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/read_coef.C

#------------------------------------------------------------------
mall_properties.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) \
                         $(FRND_ENT) \
                         $(CODE)/interface/coords_cp/mall_properties.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/mall_properties.C

#------------------------------------------------------------------
mall_coef.o     :        $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_CP) \
                         $(COORD_CP_ENT) $(COORD_CP_LOC) $(FRND_ENT) \
                         $(CODE)/interface/coords_cp/mall_coef.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords_cp/mall_coef.C

#------------------------------------------------------------------
control_set_cp_ewald.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(CPEWALD_ENT) $(CPEWALD_LOC) $(ENR_CP_LOC) \
                         $(MATH) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/cp_ewald/control_set_cp_ewald.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/cp_ewald/control_set_cp_ewald.C

#------------------------------------------------------------------
set_cp_ewald.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(CPEWALD_LOC) $(MATH) $(FRND_ENT) \
                         $(CODE)/interface/cp_ewald/set_cp_ewald.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/cp_ewald/set_cp_ewald.C

#------------------------------------------------------------------
search_base_cp.o  :      $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(SEARCH_ENT) $(VPS_LOC) $(HANDLE_ENT) \
                         $(CODE)/interface/search_base/search_base_cp.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/search_base/search_base_cp.C

#------------------------------------------------------------------
proj_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) $(CLASS_MDINTRA) \
                         $(SMPL_CP_LOC) $(SMPL_CP_ENT) $(ENR_CPCON_ENT) \
                         $(CODE)/interface/vel_sampl_cp/proj_vel_cp.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/proj_vel_cp.C

#------------------------------------------------------------------
set_vps_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
                         $(CODE)/interface/vps_params/set_vps_dict.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vps_params/set_vps_dict.C

#------------------------------------------------------------------
samp_vel_cp.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(SMPL_CP_LOC)  $(MATH) \
                         $(CODE)/interface/vel_sampl_cp/samp_vel_cp.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/samp_vel_cp.C

#------------------------------------------------------------------
control_vps_params.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) \
                         $(VPS_ENT) $(SEARCH_ENT) $(INTRA_LOC) $(VPS_LOC) \
                         $(HANDLE_ENT) $(FRND_ENT) $(COMM_WRAP) $(MATH) \
                         $(CODE)/interface/vps_params/control_vps_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vps_params/control_vps_params.C

#------------------------------------------------------------------
weight_node_gauss_hermite.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) \
                         $(VPS_LOC)  $(FRND_ENT) \
                         $(CODE)/interface/vps_params/weight_node_gauss_hermite.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vps_params/weight_node_gauss_hermite.C

#------------------------------------------------------------------
control_vc_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_cp/control_vc_smpl.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/control_vc_smpl.C

#------------------------------------------------------------------
control_vcnhc_smpl.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(SMPL_CP_ENT) $(SMPL_CP_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_cp/control_vcnhc_smpl.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/control_vcnhc_smpl.C

#------------------------------------------------------------------
control_scale_cp.o    :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP)  \
                         $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_cp/control_scale_cp.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_cp/control_scale_cp.C

#------------------------------------------------------------------
#==================================================================



