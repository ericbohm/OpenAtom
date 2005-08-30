#==================================================================
#            INTERFACE_MOL_FILES
#==================================================================


#==================================================================
MOL_PARMS = $(CODE)/interface/mol_params
DMOL_PARMS = $(DCODE)/interface/mol_params
#==================================================================


#==================================================================
control_mol_params.o  :  $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(MOL_ENT) $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MOL_PARMS)/control_mol_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/control_mol_params.C

#----------------------------------------------------------------
control_set_mol_params.o : \
                         $(STANDARD) $(DEFINES) \
                         $(TYP_PAR) $(CLASS_MDINT) $(CLASS_MDATM)\
                         $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(CLASS_GEN) $(CLASS_CP) \
                         $(MOL_LOC) $(HANDLE_ENT) $(FRND_ENT) \
                         $(MOL_PARMS)/set_params/control_set_mol_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/control_set_mol_params.C

#----------------------------------------------------------------
set_base_file_params.o : $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(MOL_LOC) \
                         $(MOL_PARMS)/set_params/set_base_file_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_base_file_params.C

#----------------------------------------------------------------
set_surf_params.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_surf_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_surf_params.C

#----------------------------------------------------------------
set_free_params.o     :  $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) $(CLASS_CP) \
                         $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(MOL_LOC) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_free_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_free_params.C

#----------------------------------------------------------------
set_mol_dict.o     :     $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) \
                         $(CLASS_MDINTRA) $(CLASS_CP) \
                         $(MOL_LOC) $(FRND_ENT) \
                         $(MOL_PARMS)/set_params/set_mol_dict.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_mol_dict.C

#----------------------------------------------------------------
set_mol_params.o     :   $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(TYP_PAR) $(CLASS_MDINTRA) \
                         $(CLASS_CP) $(CLASS_MDINT) $(CLASS_MDATM) \
                         $(CLASS_MDINTER) \
                         $(MOL_LOC) $(FRND_ENT) $(HANDLE_ENT) \
                         $(MOL_PARMS)/set_params/set_mol_params.C
	$(ECHO) $@
	$(COBJ_CARE) $(MOL_PARMS)/set_params/set_mol_params.C

#==================================================================

