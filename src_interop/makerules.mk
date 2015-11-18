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

#vpath %.cpp $(interoplib)
#
#VPATH += $(interoplib)

#$(foreach suf, $(fileTypes), $(eval vpath %$(suf) $(alldriverdirs) $(STANDARD_INC)) )
#vpath %.h $(ortho)

#pdsyev.o: pdsyev.cc
#	echo "Success! This is what should have happened!"
#	echo $(fileTypes)
#	echo $(interoplib)
#	mpicxx -c $< -o $@ -I. -I$(driver) -I$(base) -I$(STANDARD_INC) -I$(CHARMINC) -I$(topinclude) #-I$(ortho)
#	#mpicxx -c $(interoplib)/pdsyev.cpp -o pdsyev.o -I$(CHARMINC)
#	#mpicxx -c pdsyev.cpp -o pdsyev.o

# The primary target for this module
$(libinterop): $(libinterop_obj) 
	$(info-ar)
	echo "Success 2! This is what should have happened!"
	ar q $@ $^
	#$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

#Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.cc %.c, $(libinterop_src)) ) )
ifneq ($(MAKECMDGOALS),clean)
$(depFiles): CPPFLAGS += -I. -I$(main) -I$(driver) -I$(base) -I$(base)/include -I$(STANDARD_INC) -I$(interoplib) -I$(CHARMINC) -I$(topinclude)
-include $(depFiles)
-include $(libinterop_intf:.ci=.di)
endif

