# The relevant source files for this project
include $(pinymake)/make_defs/obj_files.h
libphysics_src=
libphysics_obj= $(addsuffix .o, $(basename $(libphysics_src)) ) $(LIB_OBJS)
libphysics_intf=
CPPFLAGS += -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 

# Specify the list of directories whose contents should be stripped from prerequisite lists 
# during dependency generation
DEPSTRIPDIRS += 
# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
fileTypes     = $(sort $(suffix $(libphysics_src) $(libphysics_intf)) )
#$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(STANDARD_INC)) )
vpath standard_include.h $(STANDARD_INC)
vpath %.C $(pinysrcdirs)
vpath fft_generic.f $(physics)/mathlib

# The primary target for this module
$(libphysics): $(libphysics_obj) | $(LIB_DECLS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

#----------------- Temporary legacy variable definitions ---------------------------
include $(pinymake)/make_defs/proto_files.h
include $(pinymake)/targetdeps.mk

define COBJ_CARE
$(info-cpp)
$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $(OPT_CARE) $<
endef

z3dfft_dec_noimsl.o : $(physics)/mathlib/z3dfft_dec_noimsl.f $(PROTO)

#-----------------------------------------------------------
# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(libphysics_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(libphysics_intf:.ci=.di)
endif

