# Define the base directory of the openatom tree
base := ../..
makedir := $(base)/makefiles
cfg  := $(base)/config.mk

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(cfg)

# The relevant source files for this test suite
TARGET          = pairCalcTest
libpaircalc     = libpaircalc.a
libpaircalc_src = ckPairCalculator.C pcBuilder.C \
				  ortho.C CLA_Matrix.o pcSectionManager.C orthoBuilder.C \
				  pcCreationManager.C pcCommManager.C
libpaircalc_intf= ckPairCalculator.ci ortho.ci CLA_Matrix.ci startupMessages.ci
libpaircalc_obj = $(addsuffix .o, $(basename $(libpaircalc_src)) )


other_reqd_src  = configure.C interface_hand.C piny_malloc.C friend_lib.C \
				  InstanceController.C PeList.C MapTable.C MapFile.C
pctest_src      = CP_State_GSpacePlane.C pairCalcTest.C $(other_reqd_src)
pctest_obj      = $(addsuffix .o, $(basename $(pctest_src)) )
pctest_intf     = gspace.ci pairCalcTest.ci


# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
fileTypes     = $(sort $(suffix $(pctest_src) $(pctest_intf)) )
$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(alldriverdirs) $(STANDARD_INC)) )
# Explicitly add the driver dir to the vpath for headers so that decl files including such headers
# can have their dependencies located by make
vpath %.h $(driver)
# Explicitly add piny dirs to vpaths because the test suite seems to need source from there
vpath standard_include.h $(STANDARD_INC)
vpath %.C $(pinysrcdirs)
vpath %.h $(pinyinc)/class_defs $(pinyinc)/proto_defs $(classical) $(interface) $(physics)/friend_lib


CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include -I$(STANDARD_INC) 

$(TARGET): $(pctest_obj) $(libpaircalc) $(libmath)
	$(info-ld)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	@echo "=========== Target produced from commit hash: $(REVNUM)"
	
$(libpaircalc): $(libpaircalc_obj)
	$(info-ar)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(pctest_src) $(libpaircalc_src) ) ) )
ifneq ($(MAKECMDGOALS),clean)
$(depFiles): CPPFLAGS += -I. -I$(driver) -I$(base) -I$(base)/include -I$(STANDARD_INC)
-include $(depFiles)
-include $(pctest_intf:.ci=.di)
-include $(libpaircalc_intf:.ci=.di)
include $(mathlib)/makerules.mk
endif

