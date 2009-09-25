# The relevant source files for this project
include $(physics)/piny_src_files.mk

libphysics_src= $(piny_src_files)
libphysics_obj= $(addsuffix .o, $(basename $(libphysics_src)) )
libphysics_intf=

# Identify all the object files that need to be generated carefully
libphysics_fragile_obj = $(addsuffix .o, $(basename $(libphysics_fragile_src)) )
# Change the optimization/debugging flags just for these source files
$(libphysics_fragile_obj): OPT = $(OPT_CARE)

# Specify the list of directories whose contents should be stripped from prerequisite lists 
# during dependency generation
DEPSTRIPDIRS += 
# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
vpath standard_include.h $(STANDARD_INC)
vpath %.C $(pinysrcdirs)
vpath %.h $(pinyinc)/class_defs $(pinyinc)/proto_defs $(classical) $(interface) $(physics)/friend_lib
vpath fft_generic.f $(physics)/mathlib
vpath z3dfft_dec_noimsl.f $(physics)/mathlib


# The primary target for this module
$(libphysics): $(libphysics_obj)
	$(info-ar)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

#-----------------------------------------------------------
# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(libphysics_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(libphysics_intf:.ci=.di)
endif

