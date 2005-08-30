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
                         $(CODE)/interface/parse/parse.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/parse/parse.C

#==================================================================
interface_hand.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) \
                         $(HANDLE_ENT) \
                        $(CODE)/interface/handle/interface_hand.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/handle/interface_hand.C

#------------------------------------------------------------------
search_base_class.o   :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(SEARCH_LOC) \
                         $(FRND_ENT) \
                         $(CODE)/interface/search_base/search_base_class.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/search_base/search_base_class.C

#------------------------------------------------------------------
data_base_handle.o    :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(SEARCH_ENT) $(INTER_ENT) $(INTER_LOC) \
                         $(INTRA_LOC)  $(HANDLE_ENT) $(FRND_ENT) \
                         $(CODE)/interface/search_base/data_base_handle.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/search_base/data_base_handle.C

#=========================================================================
control_sim_params.o :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_ENT) $(SIM_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(CODE)/interface/sim_params/control_sim_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/sim_params/control_sim_params.C

#------------------------------------------------------------------
set_sim_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA)  \
                         $(CLASS_CP) $(TYP_PAR) $(TYP_STAT)\
                         $(SIM_LOC) $(FRND_ENT) \
                         $(CODE)/interface/sim_params/set_sim_dict.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/sim_params/set_sim_dict.C

#------------------------------------------------------------------
set_sim_params.o  :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(TYP_PAR) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) $(SIM_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(COMM_WRAP) $(TYP_STAT)\
                         $(CODE)/interface/sim_params/set_sim_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/sim_params/set_sim_params.C

#=========================================================================
set_atm_NHC.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords/set_atm_NHC.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/set_atm_NHC.C

#------------------------------------------------------------------
read_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(COORD_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MATH) $(COMM_WRAP) \
                         $(CODE)/interface/coords/read_coord.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/read_coord.C

#------------------------------------------------------------------
molecule_decomp.o  :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINT)\
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(COORD_LOC) $(FRND_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/coords/molecule_decomp.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/molecule_decomp.C

#------------------------------------------------------------------
read_hmat.o     :        $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_ENT) $(HANDLE_ENT) $(FRND_ENT) $(MATH) \
                         $(CODE)/interface/coords/read_hmat.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/read_hmat.C

#------------------------------------------------------------------
mall_coord.o     :       $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(COORD_LOC) $(FRND_ENT) \
                         $(CODE)/interface/coords/mall_coord.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/coords/mall_coord.C

#=========================================================================
control_surf_params.o : $(STANDARD) $(DEFINES) \
                        $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                        $(TYP_PAR) $(CLASS_MDINTRA) \
                        $(CLASS_CP) \
                        $(SURF_PRMS_ENT) $(SURF_PRMS_LOC) \
                        $(INTRA_LOC) $(SEARCH_ENT) \
                        $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
                        $(CODE)/interface/surf_params/control_surf_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/surf_params/control_surf_params.C

#------------------------------------------------------------------
set_surf_dict.o   :     $(STANDARD) $(DEFINES) \
                        $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                        $(CLASS_MDINTER) \
                        $(SURF_PRMS_LOC) $(FRND_ENT) \
                        $(CODE)/interface/surf_params/set_surf_dict.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/surf_params/set_surf_dict.C

#=========================================================================
control_inter_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(INTER_ENT) $(INTER_LOC) $(INTRA_LOC) $(SEARCH_ENT) \
                         $(FRND_ENT) $(HANDLE_ENT) $(COMM_WRAP) \
                         $(CODE)/interface/inter_params/control_inter_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/inter_params/control_inter_params.C

#------------------------------------------------------------------
set_inter_dict.o   :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(CODE)/interface/inter_params/set_inter_dict.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/inter_params/set_inter_dict.C

#------------------------------------------------------------------
get_clong.o     :        $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(CODE)/interface/inter_params/get_clong.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/inter_params/get_clong.C

#------------------------------------------------------------------
spline_fit.o     :       $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(INTER_LOC) $(FRND_ENT) \
                         $(CODE)/interface/inter_params/spline_fit.C
	$(ECHO) $@
	$(COBJ) $(CODE)/interface/inter_params/spline_fit.C

#==================================================================
control_vnhc_smpl.o   :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_class/control_vnhc_smpl.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/control_vnhc_smpl.C

#------------------------------------------------------------------
control_vx_smpl.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(SMPL_CLASS_ENT) $(SMPL_CLASS_LOC) $(COORD_LOC) \
                         $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_class/control_vx_smpl.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/control_vx_smpl.C

#------------------------------------------------------------------
control_scale_class.o : $(STANDARD) $(DEFINES) \
                        $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                        $(CLASS_MDINTRA) $(CLASS_GEN) $(TYP_PAR) \
                        $(SMPL_CLASS_ENT) $(COMM_WRAP) \
                        $(CODE)/interface/vel_sampl_class/control_scale_class.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/control_scale_class.C

#------------------------------------------------------------------
proj_vel_class.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(COMM_WRAP) \
                         $(CODE)/interface/vel_sampl_class/proj_vel_class.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/proj_vel_class.C

#------------------------------------------------------------------
samp_vel_class.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(INTRA_CON_ENT) $(SMPL_CLASS_LOC) $(MATH) \
                         $(CODE)/interface/vel_sampl_class/samp_vel_class.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/vel_sampl_class/samp_vel_class.C

#==================================================================
set_exclude.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(LISTS_ENT) $(LISTS_LOC) $(FRND_ENT) $(WEIGH_NODE) \
                         $(CODE)/interface/lists/set_exclude.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/set_exclude.C

#------------------------------------------------------------------
exl_sort.o     :         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) $(TYP_PAR) \
                         $(LISTS_LOC) $(FRND_ENT) \
                         $(CODE)/interface/lists/exl_sort.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/lists/exl_sort.C

#------------------------------------------------------------------
path_integral_init.o     :         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_GEN) $(CLASS_MDINTRA) \
                         $(PATH_INIT_ENT) $(MATH) $(FRND_ENT) \
                         $(CODE)/interface/path_integral/path_integral_init.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/path_integral/path_integral_init.C

#==================================================================





