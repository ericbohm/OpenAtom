#===================================================================
#         AB INITO PHYSICS  
#===================================================================


#=================================================================
cp_eke.o    :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_NONLOC) \
               $(CODE)/abinito_physics/cp_nonlocal/cp_eke.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_nonlocal/cp_eke.C
#------------------------------------------------------------------
cp_nl_energy_forc.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) \
               $(CODE)/abinito_physics/cp_nonlocal/cp_nl_energy_forc.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_nonlocal/cp_nl_energy_forc.C

#------------------------------------------------------------------
cp_nlmat.o  :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CP_OPERATIONS) $(CLASS_CP_NONLOC) \
               $(CODE)/abinito_physics/cp_nonlocal/cp_nlmat.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_nonlocal/cp_nlmat.C
#------------------------------------------------------------------
cp_hart_ext.o :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDATM) \
               $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_LOC) \
               $(CODE)/abinito_physics/cp_local/cp_hart_ext.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_local/cp_hart_ext.C

#------------------------------------------------------------------
cp_struct_fact.o :  $(STANDARD) $(CLASS_GEN) $(CLASS_MDINT) $(CLASS_MDATM) \
               $(CLASS_MDINTER) $(CLASS_MDINTRA) $(CLASS_CP) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_LOC) $(CLASS_CP_NONLOC)\
               $(CODE)/abinito_physics/cp_nonlocal/cp_struct_fact.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_nonlocal/cp_struct_fact.C

#------------------------------------------------------------------
cp_process_grad.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(CODE)/abinito_physics/cp_xc_fnctl/cp_process_grad.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_xc_fnctl/cp_process_grad.C

#------------------------------------------------------------------
cp_grad_rho_ctrl.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(CODE)/abinito_physics/cp_xc_fnctl/cp_grad_rho_ctrl.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_xc_fnctl/cp_grad_rho_ctrl.C

#------------------------------------------------------------------
cp_becke.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(CODE)/abinito_physics/cp_xc_fnctl/cp_becke.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_xc_fnctl/cp_becke.C

#------------------------------------------------------------------
cp_pz_lda.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CPXC_FNCTLS) \
               $(CODE)/abinito_physics/cp_xc_fnctl/cp_pz_lda.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_xc_fnctl/cp_pz_lda.C

#------------------------------------------------------------------
cp_control_integrate.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE) \
               $(CODE)/abinito_physics/cp_integrate/cp_control_integrate.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_integrate/cp_control_integrate.C
#------------------------------------------------------------------
cp_dynamics.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_INTEGRATE) \
               $(CODE)/abinito_physics/cp_integrate/cp_dynamics.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_integrate/cp_dynamics.C
#------------------------------------------------------------------
cp_rspace_ion.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_RSPACEION) \
               $(CODE)/abinito_physics/cp_ions/cp_rspace_ion.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_ions/cp_rspace_ion.C
#------------------------------------------------------------------
cp_min_std.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)  \
               $(CLASS_CP_INTEGRATE) \
               $(CODE)/abinito_physics/cp_integrate/cp_min_std.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_integrate/cp_min_std.C

#------------------------------------------------------------------
cp_min_CG.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)  \
               $(CLASS_CP_INTEGRATE) \
               $(CODE)/abinito_physics/cp_integrate/cp_min_CG.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_integrate/cp_min_CG.C

#------------------------------------------------------------------
cp_lowdin.o :  $(STANDARD) $(FFTW) \
               $(COMPLEX)   \
               $(CLASS_CP_ORTHOG) $(MATH) \
               $(CODE)/abinito_physics/cp_orthog/cp_lowdin.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/abinito_physics/cp_orthog/cp_lowdin.C

#------------------------------------------------------------------
gen_wave.o    :  $(STANDARD) $(COMPLEX) $(TYP_PAR) \
               $(CLASS_GEN) $(CLASS_CP) $(CLASS_MDATM) \
               $(MATH) $(HANDLE_ENT)\
               $(GEN_WAVE) \
               $(CODE)/interface/coords_cp/gen_wave.C
	$(ECHO) $@
	$(COBJ_TEST) $(CODE)/interface/coords_cp/gen_wave.C

#===================================================================
