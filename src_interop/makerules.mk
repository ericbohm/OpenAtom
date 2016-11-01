# The relevant source files for this module
libinterop_src   = pdsyev.cc
libinterop_obj   = $(addsuffix .o, $(basename $(libinterop_src)) )
libinterop_intf  =

# Specify the list of directories whose contents should be stripped from prerequisite lists
# during dependency generation
DEPSTRIPDIRS +=
# Add appropriate directory(ies) to the vpaths for the source file types present in this module
# so that they can be located from the build directory. This is a small effort to avoid swamping
# VPATH with a long list of directories hurting the build times that we hope to improve
fileTypes     = $(sort $(suffix $(libinterop_src) $(libinterop_intf) ) )
$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(interoplib)) )

# The primary target for this module
$(libinterop): $(libinterop_obj)
	$(info-ar)
	ar q $@ $^

#Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.cc %.c, $(libinterop_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
$(depFiles): CPPFLAGS += -I. -I$(main) -I$(driver) -I$(base) -I$(base)/include -I$(STANDARD_INC) -I$(interoplib) -I$(CHARMINC) -I$(topinclude)
-include $(depFiles)
-include $(libinterop_intf:.ci=.di)
endif

