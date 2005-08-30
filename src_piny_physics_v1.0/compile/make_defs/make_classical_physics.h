#===================================================================
#         CLASSICAL PHYSICS  
#===================================================================

#------------------------------------------------------------------
control_integrate.o :  $(STANDARD) \
               $(CLASS_CHARM_ATOM) \
               $(CLASS_ATOM_INTEGRATE) \
               $(CODE)/classical_physics/integrate/control_integrate.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/classical_physics/integrate/control_integrate.C
#------------------------------------------------------------------

#===================================================================
