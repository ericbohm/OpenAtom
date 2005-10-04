#=====================================================================
include ../../compile/Makefile.config.gen
include ../../compile/Makefile.config.charm
include ../../compile/Makefile.config.piny
#=====================================================================


#===========================================================================
RM    = /bin/rm -f
ECHO  = /bin/echo "making"
HEAD  = ..
MAIN  = ../main
TOPINC  = ../../include
STATE = ../cp_state_ctrl
RHO   = ../cp_density_ctrl
STRUCTFACT   = ../structure_factor
ORTHO = ../orthog_ctrl
UTIL  = ../utility
PAIRCALC  = ../paircalc
MATH  = ../mathlib
PHYSICS = ../../src_piny_physics_v1.0/
LOADB = ../load_balance
#===========================================================================


#===========================================================================
FLAGS  = -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib \
         -I$(HEAD) -I$(MAIN) -I$(STATE) -I$(RHO) \
         -I$(ORTHO) -I$(UTIL) -I$(MATH) -I$(LOADB) \
	 -I../misc -I$(STRUCTFACT) -I$(TOPINC) \
         -I../$(STANDARD_INC) #-DGROUP_REDUCTION
LIB    = -L$(FFT_HOME)/lib -module commlib -lfftw \
         -lrfftw -lm -module fftlib                      #-memory paranoid
#===========================================================================
PCHEADERS= $(PAIRCALC)/pairCalculator.h\
	$(PAIRCALC)/pairutil.h \
	$(PAIRCALC)/ckPairCalculator.h\
	$(MAIN)/ckPairCalculator.decl.h\
	$(MAIN)/ckPairCalculator.def.h

#===========================================================================


OBJS   = ckPairCalculator.o pairCalculator.o matmul.o cpaimd.o StructFactorCache.o StructureFactor.o \
	 CP_State_GSpacePlane.o CP_State_ParticlePlane.o \
         CP_State_RealSpacePlane.o CP_Rho_RealSpacePlane.o \
	 CP_Rho_GSpacePlane.o \
         util.o ortho.o sim_subroutines.o  \
         groups.o map.o para_grp_parse.o

OBJS_DECL = ../main/cpaimd.decl.h
OBJS_DEF = ../main/cpaimd.decl.h
TARGET_EXE = ../../binary/CharmDriver.x
TARGET_LIB = ../../binary/libCharmDriver.a
#===========================================================================

#===========================================================================
CHARMLIBS:
	(cd $(CHARMBASE)/tmp/libs/ck-libs/sparseContiguousReducer;\
	make OPTS="$(OPTS)";)
#	(cd $(CHARMBASE)/tmp/libs/ck-libs/paircalc;\
#	make OPTS="$(OPTS) $(OPT_FORT)" FFTW_HOME="$(FFT_HOME)";)
	(cd $(CHARMBASE)/tmp/libs/ck-libs/fftlib;\
	make OPTS="$(OPTS)" FFTW_HOME="$(FFT_HOME)";)
#===========================================================================


#===========================================================================
summary:$(OBJS)
	$(CHARMC) -o $(TARGET_EXE).summary $(OBJS) $(LIB) \
        -language charm++ -module CkMulticast -module CkSparseContiguousReducer \
        -tracemode summary $(OPTS) $(FLAGS) $(OPT_FORT)
#===========================================================================


#===========================================================================
proj:	$(OBJS_DECL) $(OBJS)
	$(CHARMC) -o $(TARGET_EXE).proj $(OBJS) $(LIB) -language charm++ \
        -module CkMulticast -module CkSparseContiguousReducer \
        -tracemode projections $(OPTS) $(FLAGS) $(OPT_FORT)
#===========================================================================


#===========================================================================
exe:	$(OBJS_DECL) $(OBJS)
	$(CHARMC) -o $(TARGET_EXE) $(OBJS) $(LIB) -language charm++ \
        -module CkMulticast -module CkSparseContiguousReducer \
        $(OPTS) $(FLAGS) $(OPT_FORT)
#===========================================================================

#===========================================================================
../main/matmul.decl.h : $(MAIN)/matmul.ci
	$(ECHO) $@
	$(CHARMC) $(MAIN)/matmul.ci $(OPTS) $(FLAGS) $(OPT_FORT)
	/bin/cp -f matmul.decl.h $(MAIN)/matmul.decl.h
	/bin/cp -f matmul.def.h $(MAIN)/matmul.def.h
#---------------------------------------------------------------------------
matmul.o : $(MAIN)/matmul.decl.h $(MAIN)/matmul.C $(MAIN)/matmul.h \
          $(ORTHO)/ortho.h 
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/matmul.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
ckPairCalculator.decl.h: $(PAIRCALC)/ckPairCalculator.ci
	$(CHARMC) $(PAIRCALC)/ckPairCalculator.ci 

ckPairCalculator.o: $(PAIRCALC)/ckPairCalculator.C $(PAIRCALC)/ckPairCalculator.h ckPairCalculator.decl.h 
	$(CHARMC)  -c $(PAIRCALC)/ckPairCalculator.C $(OPTS) $(FLAGS)

pairCalculator.o: $(PAIRCALC)/pairCalculator.C $(PAIRCALC)/pairCalculator.h $(PAIRCALC)/ckPairCalculator.h ckPairCalculator.decl.h
	$(CHARMC) -c  $(PAIRCALC)/pairCalculator.C $(OPTS) $(FLAGS)

#---------------------------------------------------------------------------
#===========================================================================
driverlib:	CHARMLIBS $(OBJS_DECL) $(OBJS) 
	$(CHARMC) -o $(TARGET_LIB) $(OBJS) $(LIB) -language charm++ \
        -module CkMulticast -module CkSparseContiguousReducer \
        $(OPTS) $(FLAGS) $(OPT_FORT)

cpaimd:    $(OBJS_DECL) $(OBJS)
	$(CHARMC) -o cpaimd $(OBJS) $(LIB) -language charm++ \
	 -module CkMulticast -module CkSparseContiguousReducer \
	 $(OPTS) $(FLAGS) $(OPT_FORT)

#===========================================================================


#===========================================================================
../main/cpaimd.decl.h : $(MAIN)/cpaimd.ci
	$(ECHO) $@
	$(CHARMC) $(MAIN)/cpaimd.ci $(OPTS) $(FLAGS) $(OPT_FORT)
	/bin/cp -f cpaimd.decl.h $(MAIN)/cpaimd.decl.h
	/bin/cp -f cpaimd.def.h $(MAIN)/cpaimd.def.h
#---------------------------------------------------------------------------
cpaimd.o : $(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.C $(MAIN)/cpaimd.h \
	$(MAIN)/ReductionClients.h \
	$(ORTHO)/ortho.h $(HEAD)/misc/sim_subroutines.h \
	$(MAIN)/groups.h $(STATE)/CP_State_Plane.h \
	$(TOPINC)/debug_flags.h \
	$(TOPINC)/RunDescriptor.h \
	$(STRUCTFACT)/StructureFactor.h $(STRUCTFACT)/StructFactorCache.h
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/cpaimd.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
CP_State_GSpacePlane.o: $(STATE)/CP_State_GSpacePlane.C \
                        $(STATE)/CP_State_Plane.h \
                        $(MAIN)/cpaimd.decl.h \
                        $(MAIN)/cpaimd.def.h \
			$(TOPINC)/debug_flags.h \
			$(TOPINC)/RunDescriptor.h \
			$(HEAD)/misc/sim_subroutines.h \
			$(STRUCTFACT)/StructureFactor.h \
			$(STRUCTFACT)/StructFactorCache.h \
			$(MAIN)/groups.h \
                        $(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h  \
			$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpintegrate.h 

	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_GSpacePlane.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
CP_State_ParticlePlane.o: $(STATE)/CP_State_ParticlePlane.C \
                          $(STATE)/CP_State_Plane.h \
                          $(STRUCTFACT)/StructureFactor.h \
                          $(STRUCTFACT)/StructFactorCache.h \
                          $(MAIN)/cpaimd.decl.h \
                          $(MAIN)/cpaimd.def.h \
				$(TOPINC)/debug_flags.h \
				$(TOPINC)/RunDescriptor.h \
				$(HEAD)/misc/sim_subroutines.h \
				$(MAIN)/groups.h \
				$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
         $(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cplocal.h
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_ParticlePlane.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
CP_State_RealSpacePlane.o: $(STATE)/CP_State_RealSpacePlane.C \
                           $(STATE)/CP_State_Plane.h \
                           $(MAIN)/cpaimd.decl.h \
			   $(TOPINC)/debug_flags.h \
			   $(TOPINC)/RunDescriptor.h \
			   $(MAIN)/cpaimd.def.h \
			   $(HEAD)/misc/sim_subroutines.h \
			   $(MAIN)/groups.h 
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_RealSpacePlane.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
CP_Rho_RealSpacePlane.o: $(RHO)/CP_Rho_RealSpacePlane.C \
                         $(STATE)/CP_State_Plane.h \
                         $(MAIN)/cpaimd.decl.h \
	 		 $(TOPINC)/debug_flags.h \
	 		 $(TOPINC)/RunDescriptor.h \
                         $(MAIN)/cpaimd.def.h \
                         $(HEAD)/misc/sim_subroutines.h \
                         $(MAIN)/groups.h \
                         $(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
	$(ECHO) $@
	$(CHARMC) -c $(RHO)/CP_Rho_RealSpacePlane.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
CP_Rho_GSpacePlane.o: $(RHO)/CP_Rho_GSpacePlane.C \
		      $(STATE)/CP_State_Plane.h \
                      $(MAIN)/cpaimd.decl.h \
                      $(MAIN)/cpaimd.def.h \
		      $(TOPINC)/debug_flags.h \
		      $(TOPINC)/RunDescriptor.h \
		      $(HEAD)/misc/sim_subroutines.h \
		      $(MAIN)/groups.h \
		      $(STATE)/CP_State_Plane.h \
		      $(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
		      $(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h

	$(ECHO) $@
	$(CHARMC) -c $(RHO)/CP_Rho_GSpacePlane.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
util.o: $(UTIL)/util.C $(UTIL)/util.h $(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h
	$(ECHO) $@
	$(CHARMC) -c $(UTIL)/util.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
para_grp_parse.o: $(UTIL)/para_grp_parse.C $(UTIL)/para_grp_parse.h \
                  $(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h 
	$(ECHO) $@
	$(CHARMC) -c $(UTIL)/para_grp_parse.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
map.o:  $(LOADB)/map.C 	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h 
	$(ECHO) $@
	$(CHARMC) -c $(LOADB)/map.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
ortho.o: $(ORTHO)/ortho.C $(ORTHO)/ortho.h $(MAIN)/cpaimd.decl.h \
	 $(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
         $(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cporthog.h  
	$(ECHO) $@
	$(CHARMC) -c $(ORTHO)/ortho.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
StructureFactor.o: $(STRUCTFACT)/StructureFactor.h \
	$(STRUCTFACT)/StructureFactor.C \
	$(TOPINC)/debug_flags.h \
        $(TOPINC)/RunDescriptor.h \
	$(STRUCTFACT)/StructFactorCache.h \
	$(MAIN)/cpaimd.decl.h \
	$(MAIN)/cpaimd.def.h \
	$(HEAD)/misc/sim_subroutines.h \
	$(MAIN)/groups.h \
	$(MAIN)/cpaimd.h \
	$(UTIL)/util.h 
	$(ECHO) $@
	$(CHARMC) -c $(STRUCTFACT)/StructureFactor.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
StructFactorCache.o: $(HEAD)/structure_factor/StructFactorCache.h $(HEAD)/structure_factor/StructFactorCache.C \
	$(STRUCTFACT)/StructureFactor.h \
	$(MAIN)/cpaimd.decl.h \
	$(TOPINC)/debug_flags.h \
        $(TOPINC)/RunDescriptor.h \
	$(MAIN)/cpaimd.def.h \
	$(HEAD)/misc/sim_subroutines.h \
	$(MAIN)/groups.h \
	$(MAIN)/cpaimd.h \
	$(UTIL)/util.h 
	$(ECHO) $@
	$(CHARMC) -c $(STRUCTFACT)/StructFactorCache.C $(OPTS) $(FLAGS) $(OPT_FORT)

#---------------------------------------------------------------------------
sim_subroutines.o: $(HEAD)/misc/sim_subroutines.h $(HEAD)/misc/sim_subroutines.C\
	$(MAIN)/cpaimd.decl.h \
	$(TOPINC)/debug_flags.h \
        $(TOPINC)/RunDescriptor.h \
	$(MAIN)/cpaimd.decl.h 
	$(ECHO) $@
	$(CHARMC) -c $(HEAD)/misc/sim_subroutines.C $(OPTS) $(FLAGS) $(OPT_FORT)
#---------------------------------------------------------------------------
groups.o: $(MAIN)/groups.C $(MAIN)/groups.h \
          $(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h 
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/groups.C $(OPTS) $(FLAGS) $(OPT_FORT)
#===========================================================================


#===========================================================================
clean:	
	$(RM)  ../conv-host *.o ../charmrun
	$(RM) $(HEAD)/*.def.h $(HEAD)/*.decl.h   $(HEAD)/*~
	$(RM) $(MAIN)/*.def.h $(MAIN)/*.decl.h   $(MAIN)/*~
	$(RM) $(STATE)/*.def.h $(STATE)/*.decl.h $(STATE)/*~
	$(RM) $(RHO)/*.def.h $(RHO)/*.decl.h     $(RHO)/*~
	$(RM) $(ORTHO)/*.def.h $(ORTHO)/*.decl.h $(ORTHO)/*~
	$(RM) $(UTIL)/*.def.h $(UTIL)/*.decl.h   $(UTIL)/*~
	$(RM) $(MATH)/*.def.h $(MATH)/*.decl.h   $(MATH)/*~
	$(RM) $(LOADB)/*.def.h $(LOADB)/*.decl.h $(LOADB)/*~
	$(RM) $(STRUCTFACT)/*.def.h $(STRUCTFACT)/*.decl.h $(STRUCTFACT)/*~
	$(RM) $(TARGET_EXE) $(TARGET_LIB) ../core *.def.h *.decl.h
#===========================================================================
