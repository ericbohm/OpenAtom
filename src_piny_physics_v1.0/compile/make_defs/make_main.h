#=================================================================
#              MAIN_FILES 
#=================================================================


#==========================================================================
Interface_ctrl.o     :           $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) \
                         $(LIB_DECLS) \
                         $(CODE)/piny_to_driver/Interface_ctrl.C \
                         $(CODE)/piny_to_driver/Interface_ctrl.h
	$(ECHO) $@
	cp -f $(CODE)/piny_to_driver/Interface_ctrl.h $(CODE)/include/class_defs
	$(COBJ) $(CODE)/piny_to_driver/Interface_ctrl.C

#==========================================================================
Interface_ctrl.decl.h  : $(CODE)/piny_to_driver/Interface_ctrl.ci
	$(ECHO) $@
	$(COBJ_DECL) $(CODE)/piny_to_driver/Interface_ctrl.ci
	test ! -d  $(CODE)/include/charm_defs && mkdir $(CODE)/include/charm_defs || true
	/bin/cp -f Interface_ctrl.decl.h $(CODE)/include/charm_defs
	/bin/cp -f Interface_ctrl.def.h $(CODE)/include/charm_defs

#==========================================================================
PhysicsAtomPosInit.o     : $(STANDARD) $(DEFINES) \
                           $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                           $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                           $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) \
                           $(CLASS_PhysicsAtomPosinit) \
                           $(CODE)/piny_to_driver/PhysicsAtomPosInit.C 
	$(ECHO) $@
	$(COBJ) $(CODE)/piny_to_driver/PhysicsAtomPosInit.C
#==========================================================================
vx_smpl.o     :          $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDINT) \
                         $(CLASS_MDATM) $(CLASS_MDINTER) $(CLASS_MDINTRA) \
                         $(MAIN_ENT) $(MAIN_LOC) $(PARSE_ENT) $(CLASS_VX_SMPL) \
                         $(CODE)/classical_physics/vel_smpl_atms/vx_smpl.C 
	$(ECHO) $@
	$(COBJ) $(CODE)/classical_physics/vel_smpl_atms/vx_smpl.C
                                                                               
#==========================================================================
ParaInfoInit.o     :      $(STANDARD) $(DEFINES) \
                         $(CLASS_GEN) $(CLASS_CP) $(DBG_FLAGS) \
                         $(CLASS_PARAINIT) \
                         $(CODE)/piny_to_driver/Parainfoinit.C
	$(ECHO) $@
	$(COBJ) $(CODE)/piny_to_driver/Parainfoinit.C
#==========================================================================




