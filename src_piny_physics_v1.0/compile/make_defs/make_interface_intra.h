#==================================================================
#               INTERFACE_INTRA_FILES
#==================================================================


#==================================================================
SET_PARM  = $(CODE)/interface/intra_params/set_params
DSET_PARM = $(DCODE)/interface/intra_params/set_params
#==================================================================


#==================================================================
close_intra_params.o :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) $(CLASS_GEN) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/close_intra_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/close_intra_params.C

#------------------------------------------------------------------
control_intra_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_ENT) $(INTRA_LOC) $(FRND_ENT) $(PIMD_LOC) \
                         $(MATH) \
                         $(CODE)/interface/intra_params/control_intra_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/control_intra_params.C

#------------------------------------------------------------------
control_res_params.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/control_res_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/control_res_params.C

#------------------------------------------------------------------
fetch_residue.o     :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/fetch_residue.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_residue.C

#------------------------------------------------------------------
fetch_resbond_prm.o  :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) \
                         $(CODE)/interface/intra_params/fetch_resbond_prm.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_resbond_prm.C

#------------------------------------------------------------------
fetch_free_energy_index.o     : \
                         $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) \
                       $(CODE)/interface/intra_params/fetch_free_energy_index.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_free_energy_index.C

#------------------------------------------------------------------
fetch_freeze.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/fetch_freeze.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_freeze.C

#------------------------------------------------------------------
init_intra_params.o  :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/init_intra_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/init_intra_params.C

#------------------------------------------------------------------
manipulate_res_bonds.o : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/manipulate_res_bonds.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/manipulate_res_bonds.C

#------------------------------------------------------------------
replicate_mol.o     :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER)\
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/replicate_mol.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/replicate_mol.C
#------------------------------------------------------------------
residue_bond.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_INTRA) $(CLASS_MDINT)\
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/residue_bond.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/residue_bond.C

#==================================================================
set_atm_mask.o     :     $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_mask.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_atm_mask.C

#------------------------------------------------------------------
set_atm_morph.o     :    $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_morph.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_atm_morph.C

#------------------------------------------------------------------
set_atm_params.o     :   $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_atm_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_atm_params.C

#------------------------------------------------------------------
set_bend_bnd_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bend_bnd_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_bend_bnd_params.C

#------------------------------------------------------------------
set_bend_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bend_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_bend_params.C

#------------------------------------------------------------------
set_bond_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_bond_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_bond_params.C

#------------------------------------------------------------------
set_intra_dict.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(SET_PARM)/set_intra_dict.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_intra_dict.C

#------------------------------------------------------------------
set_intra_dict_pot.o  :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(SET_PARM)/set_intra_dict_pot.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_intra_dict_pot.C

#------------------------------------------------------------------
set_intra_potent.o  :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(SEARCH_ENT) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/set_intra_potent.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/set_intra_potent.C
#------------------------------------------------------------------
intra_coefs.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(CODE)/interface/intra_params/intra_coefs.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/intra_coefs.C

#------------------------------------------------------------------
set_mol_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_mol_name_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_mol_name_params.C

#------------------------------------------------------------------
set_onfo_params.o     :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_onfo_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_onfo_params.C

#------------------------------------------------------------------
set_res_bond_params.o  : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_bond_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_bond_params.C

#------------------------------------------------------------------
set_res_def_params.o   : $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_def_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_def_params.C

#------------------------------------------------------------------
set_res_name_params.o :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_name_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_name_params.C

#------------------------------------------------------------------
set_res_morph_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) \
                         $(SET_PARM)/set_res_morph_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_res_morph_params.C

#------------------------------------------------------------------
set_grp_con_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_GEN) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(SET_PARM)/set_grp_con_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_grp_con_params.C

#------------------------------------------------------------------
set_tors_params.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(SET_PARM)/set_tors_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(SET_PARM)/set_tors_params.C

#------------------------------------------------------------------
fetch_hydrog_mass.o :    $(STANDARD) $(DEFINES) \
                         $(CLASS_MDINT) $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(TYP_PAR) $(CLASS_INTRA) \
                         $(INTRA_LOC) $(FRND_ENT) \
                         $(CODE)/interface/intra_params/fetch_hydrog_mass.C
	$(ECHO) $@
	$(COBJ_CARE) $(CODE)/interface/intra_params/fetch_hydrog_mass.C

#==================================================================
