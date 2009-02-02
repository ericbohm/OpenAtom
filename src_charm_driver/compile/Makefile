#===========================================================================
include ../../compile/Makefile.config
#===========================================================================


#===========================================================================
RM		= /bin/rm -f
ECHO		= /bin/echo "making"
HEAD		= ..
MAIN		= ../main
TOPINC		= ../../include
UBER		= ../uber
STATE		= ../cp_state_ctrl
RHO		= ../cp_density_ctrl
STRUCTFACT	= ../structure_factor
ORTHO		= ../orthog_ctrl
UTIL		= ../utility
FFTSLAB		= ../fft_slab_ctrl
PAIRCALC	= ../paircalc
MATH		= ../mathlib
PHYSICS		= ../../src_piny_physics_v1.0/
LOADB		= ../load_balance
SRC_COMPILE	= ../compile
#===========================================================================


#===========================================================================
FLAGS	= \
	-I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib -I$(HEAD) \
	-I$(MAIN) -I$(STATE) -I$(RHO) -I$(ORTHO) -I$(UTIL) -I$(MATH) \
	-I$(LOADB) -I$(FFTSLAB) -I$(STRUCTFACT) -I$(TOPINC) -I$(UBER) \
        -I$(PAIRCALC) -I$(SRC_COMPILE) -I../$(STANDARD_INC) 

LIB	= -L$(LOCAL_LIB_PATH) -L$(FFT_HOME)/lib -module commlib \
	-lm $(LOCAL_LIB) 
#===========================================================================


#===========================================================================
PCHEADERS = \
	$(PAIRCALC)/pairCalculator.h \
	$(PAIRCALC)/pairutil.h \
	$(PAIRCALC)/paircalcMessages.h \
	$(PAIRCALC)/MessageDataCollator.h \
	$(PAIRCALC)/InputDataHandler.h \
	$(PAIRCALC)/ForwardPathPostCollationTriggers.h \
	$(PAIRCALC)/ckPairCalculator.h \
	$(MAIN)/paircalcMessages.decl.h \
	$(MAIN)/paircalcMessages.def.h \
	$(MAIN)/RDMAMessages.decl.h \
	$(MAIN)/RDMAMessages.def.h \
	$(MAIN)/inputDataHandler.decl.h \
	$(MAIN)/inputDataHandler.def.h \
	$(MAIN)/ckPairCalculator.decl.h \
	$(MAIN)/ckPairCalculator.def.h
#===========================================================================


#===========================================================================
OBJS	= ckPairCalculator.o pairCalculator.o CLA_Matrix.o cpaimd.o \
	StructFactorCache.o  StructureFactor.o CP_State_GSpacePlane.o GSpaceDriver.o \
	CP_State_ParticlePlane.o CP_State_RealParticlePlane.o \
	CP_State_RealSpacePlane.o CP_Rho_RealSpacePlane.o \
	CP_Rho_GSpacePlane.o CP_Rho_GHartExt.o CP_Rho_RHartExt.o util.o \
	ortho.o fftCache.o stateSlab.o rhoSlab.o eesCache.o \
	groups.o MapTable.o PeList.o para_grp_parse.o MapFile.o \
	InstanceController.o

OBJS_DECL	= $(MAIN)/gSpaceDriver.decl.h $(MAIN)/ortho.decl.h $(MAIN)/structureFactor.decl.h \
		  $(MAIN)/cpaimd.decl.h $(MAIN)/groups.decl.h $(MAIN)/RDMAMessages.decl.h $(MAIN)/paircalcMessages.decl.h \
		  $(MAIN)/inputDataHandler.decl.h $(MAIN)/ckPairCalculator.decl.h 
OBJS_DEF 	= $(MAIN)/cpaimd.decl.h
TARGET_EXE	= ../../binary/CharmDriver.x
TARGET_LIB 	= ../../binary/libCharmDriver.a
#===========================================================================


#===========================================================================
CHARMLIBS:
#	(cd $(CHARMBASE)/tmp/libs/ck-libs/sparseContiguousReducer;\
#	make OPTS="$(OPTS)";)
#	(cd $(CHARMBASE)/tmp/libs/ck-libs/paircalc;\
#	make OPTS="$(OPTS) " FFTW_HOME="$(FFT_HOME)";)
#	(cd $(CHARMBASE)/tmp/libs/ck-libs/fftlib;\
#	make OPTS="$(OPTS) -nobs " FFTW_HOME="$(FFT_HOME)";)
#===========================================================================


#===========================================================================
summary:$(OBJS)
	$(CHARMC) -o $(TARGET_EXE).summary $(OBJS) $(LIB) \
        -language charm++ -module CkMulticast \
        -tracemode summary $(OPTS) $(FLAGS)
#===========================================================================


#===========================================================================
proj:	$(OBJS_DECL) $(OBJS)
	$(CHARMC) -o $(TARGET_EXE).proj $(OBJS) $(LIB) -language charm++ \
        -module CkMulticast \
        -tracemode projections $(OPTS) $(FLAGS)
#===========================================================================


#===========================================================================
exe:	$(OBJS_DECL) $(OBJS)
	$(CHARMC) -o $(TARGET_EXE) $(OBJS) $(LIB) -language charm++ \
        -module CkMulticast \
        $(OPTS) $(FLAGS)
#===========================================================================


#===========================================================================
$(MAIN)/cpaimd.decl.h $(MAIN)/groups.decl.h $(MAIN)/CPcharmParaInfo.decl.h: $(MAIN)/cpaimd.ci
	$(ECHO) $@
	$(CHARMC) $(MAIN)/cpaimd.ci $(OPTS) $(FLAGS)
	cp cpaimd.decl.h cpaimd.def.h $(MAIN)
	cp groups.decl.h groups.def.h $(MAIN)
	cp CPcharmParaInfo.decl.h CPcharmParaInfo.def.h $(MAIN)
#---------------------------------------------------------------------------
$(MAIN)/CLA_Matrix.decl.h : $(MAIN)/CLA_Matrix.ci
	$(ECHO) $@
	$(CHARMC) $(MAIN)/CLA_Matrix.ci $(OPTS) $(FLAGS)
	/bin/cp -f CLA_Matrix.decl.h $(MAIN)/CLA_Matrix.decl.h
	/bin/cp -f CLA_Matrix.def.h $(MAIN)/CLA_Matrix.def.h
#---------------------------------------------------------------------------
$(MAIN)/paircalcMessages.decl.h $(MAIN)/RDMAMessages.decl.h $(MAIN)/inputDataHandler.decl.h $(MAIN)/ckPairCalculator.decl.h : $(PAIRCALC)/ckPairCalculator.ci
	$(CHARMC) $(PAIRCALC)/ckPairCalculator.ci 
	cp paircalcMessages.decl.h paircalcMessages.def.h $(MAIN)
	cp RDMAMessages.decl.h RDMAMessages.def.h $(MAIN)
	cp inputDataHandler.decl.h inputDataHandler.def.h $(MAIN)
	cp ckPairCalculator.decl.h ckPairCalculator.def.h $(MAIN)
#---------------------------------------------------------------------------
$(MAIN)/gspace.decl.h $(MAIN)/gSpaceDriver.decl.h: $(STATE)/gspace.ci
	$(ECHO) $@
	$(CHARMC) $(STATE)/gspace.ci
	cp gSpaceDriver.decl.h gSpaceDriver.def.h $(MAIN)
	cp gStatePlane.decl.h gStatePlane.def.h $(MAIN)
	cp gParticlePlane.decl.h gParticlePlane.def.h $(MAIN)
#---------------------------------------------------------------------------
$(MAIN)/ortho.decl.h: $(ORTHO)/ortho.ci
	$(ECHO) $@
	$(CHARMC) $^
	cp ortho.decl.h ortho.def.h $(MAIN)
#---------------------------------------------------------------------------
$(MAIN)/structureFactor.decl.h: $(STRUCTFACT)/structureFactor.ci 
	$(ECHO) $@
	$(CHARMC) $(STRUCTFACT)/structureFactor.ci
	cp structureFactor.decl.h structureFactor.def.h $(MAIN)
	cp structureFactorCache.decl.h structureFactorCache.def.h $(MAIN)
#---------------------------------------------------------------------------
CLA_Matrix.o : $(MAIN)/CLA_Matrix.decl.h $(MAIN)/CLA_Matrix.C \
	$(MAIN)/CLA_Matrix.h $(ORTHO)/ortho.h 
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/CLA_Matrix.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
ckPairCalculator.o: $(PAIRCALC)/ckPairCalculator.C \
	$(PAIRCALC)/ckPairCalculator.h $(MAIN)/paircalcMessages.decl.h \
	$(MAIN)/inputDataHandler.decl.h $(MAIN)/ckPairCalculator.decl.h 
	$(CHARMC)  -c $(PAIRCALC)/ckPairCalculator.C $(OPTS) $(FLAGS)

pairCalculator.o: $(PAIRCALC)/pairCalculator.C \
	$(PAIRCALC)/pairCalculator.h $(PAIRCALC)/ckPairCalculator.h \
	$(MAIN)/paircalcMessages.decl.h $(MAIN)/RDMAMessages.decl.h \
	$(MAIN)/inputDataHandler.decl.h $(MAIN)/ckPairCalculator.decl.h
	$(CHARMC) -c  $(PAIRCALC)/pairCalculator.C $(OPTS) $(FLAGS)
#===========================================================================


#===========================================================================
driverlib: CHARMLIBS $(OBJS_DECL) $(OBJS) 
	$(CHARMC) -o $(TARGET_LIB) $(OBJS) $(LIB) -language charm++ \
        -module CkMulticast $(OPTS) $(FLAGS)

cpaimd: $(OBJS_DECL) $(OBJS)
	$(CHARMC) -o cpaimd $(OBJS) $(LIB) -language charm++ \
	-module CkMulticast $(OPTS) $(FLAGS)
#===========================================================================


#===========================================================================
cpaimd.o: $(MAIN)/gSpaceDriver.decl.h $(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.C $(MAIN)/cpaimd.h \
	$(MAIN)/ReductionClients.h $(MAIN)/TimeKeeper.h \
	$(MAIN)/groups.h $(STATE)/CP_State_Plane.h \
	$(ORTHO)/ortho.h $(FFTSLAB)/fftCacheSlab.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/eesDataClass.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(LOADB)/PeList.h $(LOADB)/MapTable.h $(LOADB)/IntMap.h \
	$(UTIL)/MapFile.h \
	$(STRUCTFACT)/StructureFactor.h $(STRUCTFACT)/StructFactorCache.h
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/cpaimd.C $(FLAGS)
#---------------------------------------------------------------------------
CP_State_GSpacePlane.o: $(STATE)/CP_State_GSpacePlane.C \
	$(STATE)/CP_State_Plane.h $(MAIN)/RDMAMessages.decl.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h $(MAIN)/groups.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/eesDataClass.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h \
	$(STRUCTFACT)/StructureFactor.h $(STRUCTFACT)/StructFactorCache.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h  \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpintegrate.h 
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_GSpacePlane.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
GSpaceDriver.o: $(STATE)/GSpaceDriver.C $(STATE)/GSpaceDriver.h $(STATE)/GSpaceRTH.h \
	$(MAIN)/gSpaceDriver.decl.h $(TOPINC)/debug_flags.h 
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/GSpaceDriver.C $(OPTS) $(FLAGS)

#---------------------------------------------------------------------------
CP_State_ParticlePlane.o: $(STATE)/CP_State_ParticlePlane.C \
	$(STATE)/CP_State_Plane.h \
	$(STRUCTFACT)/StructureFactor.h $(STRUCTFACT)/StructFactorCache.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h $(MAIN)/groups.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/eesDataClass.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cplocal.h
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_ParticlePlane.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
CP_State_RealParticlePlane.o: $(STATE)/CP_State_RealParticlePlane.C \
	$(STATE)/CP_State_Plane.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h $(MAIN)/groups.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/eesDataClass.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cplocal.h
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_RealParticlePlane.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
CP_State_RealSpacePlane.o: $(STATE)/CP_State_RealSpacePlane.C \
	$(STATE)/CP_State_Plane.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h $(MAIN)/groups.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
	$(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h
	$(ECHO) $@
	$(CHARMC) -c $(STATE)/CP_State_RealSpacePlane.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
CP_Rho_RealSpacePlane.o: $(RHO)/CP_Rho_RealSpacePlane.C \
	$(STATE)/CP_State_Plane.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h $(MAIN)/groups.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
	$(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
	$(ECHO) $@
	$(CHARMC) -c $(RHO)/CP_Rho_RealSpacePlane.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
CP_Rho_GSpacePlane.o: $(RHO)/CP_Rho_GSpacePlane.C \
	$(STATE)/CP_State_Plane.h $(STATE)/CP_State_Plane.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
	$(MAIN)/groups.h $(MAIN)/cpaimd.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
	$(ECHO) $@
	$(CHARMC) -c $(RHO)/CP_Rho_GSpacePlane.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
CP_Rho_GHartExt.o: $(RHO)/CP_Rho_GHartExt.C \
	$(STATE)/CP_State_Plane.h $(STATE)/CP_State_Plane.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h \
	$(MAIN)/groups.h $(MAIN)/cpaimd.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/eesDataClass.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
	$(ECHO) $@
	$(CHARMC) -c $(RHO)/CP_Rho_GHartExt.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
CP_Rho_RHartExt.o: $(RHO)/CP_Rho_RHartExt.C \
	$(STATE)/CP_State_Plane.h $(STATE)/CP_State_Plane.h \
	$(MAIN)/cpaimd.decl.h $(MAIN)/cpaimd.def.h \
	$(MAIN)/groups.h $(MAIN)/cpaimd.h \
	$(TOPINC)/CPcharmParaInfo.h $(TOPINC)/eesDataClass.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(FFTSLAB)/fftCacheSlab.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpnonlocal.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cpxcfnctls.h
	$(ECHO) $@
	$(CHARMC) -c $(RHO)/CP_Rho_RHartExt.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
util.o: $(UTIL)/util.C $(UTIL)/util.h $(TOPINC)/debug_flags.h \
	$(TOPINC)/RunDescriptor.h $(TOPINC)/CPcharmParaInfo.h   
	$(ECHO) $@
	$(CHARMC) -c $(UTIL)/util.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
MapFile.o: $(UTIL)/MapFile.C $(UTIL)/MapFile.h $(TOPINC)/debug_flags.h \
	$(TOPINC)/RunDescriptor.h $(TOPINC)/CPcharmParaInfo.h \
	$(LOADB)/MapTable.h
	$(ECHO) $@
	$(CHARMC) -c $(UTIL)/MapFile.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
para_grp_parse.o: $(UTIL)/para_grp_parse.C $(UTIL)/para_grp_parse.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(TOPINC)/CPcharmParaInfo.h   
	$(ECHO) $@
	$(CHARMC) -c $(UTIL)/para_grp_parse.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
PeList.o: $(LOADB)/PeList.C $(LOADB)/PeList.h $(TOPINC)/debug_flags.h \
	$(TOPINC)/RunDescriptor.h $(TOPINC)/CPcharmParaInfo.h
	$(ECHO) $@
	$(CHARMC) -c $(LOADB)/PeList.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
MapTable.o: $(LOADB)/MapTable.C $(LOADB)/MapTable.h $(TOPINC)/debug_flags.h\
	$(TOPINC)/RunDescriptor.h $(TOPINC)/CPcharmParaInfo.h \
	$(LOADB)/IntMap.h
	$(ECHO) $@
	$(CHARMC) -c $(LOADB)/MapTable.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
ortho.o: $(ORTHO)/ortho.C $(ORTHO)/ortho.h $(MAIN)/ortho.decl.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(PHYSICS)/include/class_defs/CP_OPERATIONS/class_cporthog.h \
	$(TOPINC)/CPcharmParaInfo.h 
	$(ECHO) $@
	$(CHARMC) -c $(ORTHO)/ortho.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
StructureFactor.o: $(STRUCTFACT)/StructureFactor.h \
	$(STRUCTFACT)/StructureFactor.C $(STRUCTFACT)/StructFactorCache.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h $(MAIN)/structureFactor.decl.h \
	$(MAIN)/groups.h $(FFTSLAB)/fftCacheSlab.h $(UTIL)/util.h 
	$(ECHO) $@
	$(CHARMC) -c $(STRUCTFACT)/StructureFactor.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
StructFactorCache.o: $(STRUCTFACT)/StructFactorCache.C $(STRUCTFACT)/StructFactorCache.h $(STRUCTFACT)/StructureFactor.h \
	  $(MAIN)/structureFactor.decl.h $(MAIN)/groups.h $(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
      $(TOPINC)/RunDescriptor.h $(FFTSLAB)/fftCacheSlab.h $(UTIL)/util.h 
	$(ECHO) $@
	$(CHARMC) -c $(STRUCTFACT)/StructFactorCache.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
fftCache.o: $(FFTSLAB)/fftCacheSlab.h $(FFTSLAB)/fftCache.C \
        $(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
        $(TOPINC)/RunDescriptor.h \
	$(MAIN)/cpaimd.decl.h
	$(ECHO) $@
	$(CHARMC) -c $(FFTSLAB)/fftCache.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
eesCache.o: $(MAIN)/eesCache.h $(MAIN)/eesCache.C \
        $(TOPINC)/eesDataClass.h $(TOPINC)/CPcharmParaInfo.h \
	$(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
	$(MAIN)/cpaimd.decl.h 
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/eesCache.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
stateSlab.o: $(FFTSLAB)/fftCacheSlab.h $(FFTSLAB)/stateSlab.C\
        $(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
        $(TOPINC)/RunDescriptor.h \
	$(MAIN)/cpaimd.decl.h 
	$(ECHO) $@
	$(CHARMC) -c $(FFTSLAB)/stateSlab.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
rhoSlab.o: $(FFTSLAB)/fftCacheSlab.h $(FFTSLAB)/rhoSlab.C\
        $(TOPINC)/CPcharmParaInfo.h $(TOPINC)/debug_flags.h \
        $(TOPINC)/RunDescriptor.h \
	$(MAIN)/cpaimd.decl.h 
	$(ECHO) $@
	$(CHARMC) -c $(FFTSLAB)/rhoSlab.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
groups.o: $(MAIN)/groups.C $(MAIN)/groups.h \
        $(TOPINC)/debug_flags.h $(TOPINC)/RunDescriptor.h \
        $(TOPINC)/CPcharmParaInfo.h
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/groups.C $(OPTS) $(FLAGS)
#---------------------------------------------------------------------------
InstanceController.o: $(MAIN)/InstanceController.h  $(MAIN)/InstanceController.C \
	$(MAIN)/cpaimd.decl.h 
	$(ECHO) $@
	$(CHARMC) -c $(MAIN)/InstanceController.C $(OPTS) $(FLAGS)

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
	$(RM) $(FFTSLAB)/*.def.h $(FFTSLAB)/*.decl.h   $(FFTSLAB)/*~
	$(RM) $(MATH)/*.def.h $(MATH)/*.decl.h   $(MATH)/*~
	$(RM) $(LOADB)/*.def.h $(LOADB)/*.decl.h $(LOADB)/*~
	$(RM) $(STRUCTFACT)/*.def.h $(STRUCTFACT)/*.decl.h $(STRUCTFACT)/*~
	$(RM) $(TARGET_EXE) $(TARGET_LIB) ../core *.def.h *.decl.h
#===========================================================================
