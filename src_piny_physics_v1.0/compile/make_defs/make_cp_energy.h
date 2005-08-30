#==================================================================
#    CP ENERGY ROUTINES 
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
#=====================================================================
cp_energy_control.o   :  $(STANDARD) $(DEFINES) \
                         $(TYP_GEN) $(TYP_CLASS) $(TYP_BND) $(TYP_CP) \
                         $(ENR_CTRL_ENT) $(ENR_CTRL_LOC) \
                         $(ENR_CTRL_CP_ENT) $(ENR_CTRL_CP_LOC) \
                         $(INTRA_ENT) $(REAL_ENT) $(REC_ENT) \
                         $(ENR_CTRL_CP_ENT) $(COMM_WRAP) \
                         $(CODE)/energy/control_cp/cp_energy_control.c
        $(ECHO) $@
        $(COBJ) $(CODE)/energy/control_cp/cp_energy_control.c
#=====================================================================
