# Define the base directory of the openatom tree
base := ../..

# Include basic make rules
include $(base)/compile/srcdirs.mk
include $(base)/compile/basicrules.mk
include $(base)/compile/Makefile.config

# Directory structure and vpath related stuff
# Overriding values set in srcdirs.mk for a shorter, cleaner paths
driver        = ..
build_drv     = .
CPAIMD        = $(base)
PINY_STD_INC  = $(STANDARD_INC)
DEPSTRIPDIRS += 
VPATH         = $(alldriverdirs) $(PINY_STD_INC) 

# Specify the compilation target
TARGET    = libCharmDriver.a

# The relevant source files for this project
SRC       = \
            cpaimd.C InstanceController.C groups.C eesCache.C CLA_Matrix.C \
			CP_State_GSpacePlane.C CP_State_RealSpacePlane.C GSpaceDriver.C \
			CP_LargeSP_RhoGSpacePlane.C  CP_LargeSP_RhoRealSpacePlane.C \
			CP_State_ParticlePlane.C CP_State_RealParticlePlane.C \
			CP_Rho_RealSpacePlane.C CP_Rho_GSpacePlane.C CP_Rho_GHartExt.C CP_Rho_RHartExt.C \
            ckPairCalculator.C pairCalculator.C \
			ortho.C pcSectionManager.C \
            StructureFactor.C StructFactorCache.C \
			fftCache.C stateSlab.C rhoSlab.C \
            MapTable.C PeList.C \
			util.C MapFile.C para_grp_parse.C
            
HDR       = 
OBJ       = $(addsuffix .o, $(basename $(SRC)) )
INTF      = ckPairCalculator.ci ortho.ci gspace.ci CLA_Matrix.ci structureFactor.ci cpaimd.ci 

########### Specify the flags, include paths, libraries etc.###########
# Preprocessor flags (Used for all preprocessing, compilation and linking)
CPPFLAGS += -I$(driver) -I$(build_drv) -I$(base) -I$(base)/include -I$(PINY_STD_INC)
# Compiler flags (Used for compilation and linking)
CFLAGS   += 
CXXFLAGS += 
FFLAGS   += 
# Linker flags (Used for all linking)
LDFLAGS  += 
# Additional libraries
LDLIBS   += 

.PHONY: all clean again 

all: $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.a) $(TARGET) 

again: 
	$(MAKE) clean; $(MAKE)

# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(SRC)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(INTF:.ci=.di)
endif

