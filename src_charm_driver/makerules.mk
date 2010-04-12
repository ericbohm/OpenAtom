# The relevant source files for this module
libdriver_src = \
                cpaimd.C InstanceController.C pcCreationManager.C groups.C eesCache.C CLA_Matrix.C \
                CP_State_GSpacePlane.C CP_State_RealSpacePlane.C GSpaceDriver.C pcCommManager.C \
                CP_LargeSP_RhoGSpacePlane.C  CP_LargeSP_RhoRealSpacePlane.C \
                CP_State_ParticlePlane.C CP_State_RealParticlePlane.C \
                CP_Rho_RealSpacePlane.C CP_Rho_GSpacePlane.C CP_Rho_GHartExt.C CP_Rho_RHartExt.C \
                CP_VanderWaalsR.C CP_VanderWaalsG.C \
                ckPairCalculator.C pcBuilder.C \
                ortho.C pcSectionManager.C orthoBuilder.C \
                StructureFactor.C StructFactorCache.C \
                fftCache.C stateSlab.C rhoSlab.C \
                MapTable.C PeList.C \
                util.C MapFile.C para_grp_parse.C matrix2file.C
libdriver_obj = $(addsuffix .o, $(basename $(libdriver_src)) )
libdriver_intf= ckPairCalculator.ci ortho.ci gspace.ci CLA_Matrix.ci structureFactor.ci startupMessages.ci cpaimd.ci 


# Specify the list of directories whose contents should be stripped from prerequisite lists 
# during dependency generation
DEPSTRIPDIRS += 
# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
fileTypes     = $(sort $(suffix $(libdriver_src) $(libdriver_intf)) )
$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(alldriverdirs) $(STANDARD_INC)) )
# Explicitly add the driver dir to the vpath for headers so that decl files including such headers
# can have their dependencies located by make
vpath %.h $(driver)

# The primary target for this module
$(libdriver): $(libdriver_obj)
	$(info-ar)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(libdriver_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
$(depFiles): CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include -I$(STANDARD_INC)
-include $(depFiles)
-include $(libdriver_intf:.ci=.di)
endif

