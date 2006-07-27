#===================================================================
#         CLASSICAL PHYSICS  
#===================================================================

#------------------------------------------------------------------
control_integrate.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM)   $(DBG_FLAGS) \
               $(CLASS_ATOM_INTEGRATE) \
               $(CODE)/classical_physics/integrate/control_integrate.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/classical_physics/integrate/control_integrate.C
#------------------------------------------------------------------
integration_drivers.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_INTEGRATE)   $(DBG_FLAGS) \
               $(CODE)/classical_physics/integrate/integration_drivers.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/classical_physics/integrate/integration_drivers.C
#------------------------------------------------------------------
write_gen_header.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_OUTPUT)   $(DBG_FLAGS) \
               $(CODE)/classical_physics/output/write_gen_header.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/classical_physics/output/write_gen_header.C
#------------------------------------------------------------------
write_output.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_OUTPUT)   $(DBG_FLAGS) \
               $(CODE)/classical_physics/output/write_output.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/classical_physics/output/write_output.C
#------------------------------------------------------------------

#===================================================================
