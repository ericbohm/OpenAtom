# Define the base directory of the openatom tree
base := ../..

# Include basic make rules
include $(base)/compile/srcdirs.mk
include $(base)/compile/basicrules.mk
include $(base)/compile/Makefile.config

# Directory structure and vpath related stuff
# Overriding values set in srcdirs.mk for a shorter, cleaner paths
physics       = ..
build_phy     = .
CPAIMD        = $(base)
DEPSTRIPDIRS += 
VPATH        += 

# Specify the compilation target
TARGET    = libPinyInterface.a

# The relevant source files for this project
include $(build_phy)/make_defs/obj_files.h
SRC       = 
HDR       = 
OBJ       = $(addsuffix .o, $(basename $(SRC)) ) $(LIB_OBJS)
INTF      = 

########### Specify the flags, include paths, libraries etc.###########
# Preprocessor flags (Used for all preprocessing, compilation and linking)
CPPFLAGS += -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 
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

$(TARGET): $(OBJ) | $(LIB_DECLS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.a) $(TARGET) 

again: 
	$(MAKE) clean; $(MAKE)

#----------------- Temporary legacy variable definitions ---------------------------
include $(MAKE_PHYS_INC)/make_head.h
include $(build_phy)/make_defs/proto_files.h
include $(build_phy)/targetdeps.mk

FOBJ          = $(FC) $(FFLAGS) $(FOPTS) $(CPPFLAGS) -c -o $@
COBJ          = $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $(OPT_FULL)
COBJ_CARE     = $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $(OPT_CARE)
COBJ_DECL     = $(CC) $(CFLAGS) 

z3dfft_dec_noimsl.o : $(PROTO) $(CODE)/mathlib/z3dfft_dec_noimsl.f
	$(FOBJ) $(CODE)/mathlib/z3dfft_dec_noimsl.f

#-----------------------------------------------------------
# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(SRC)) ) )
ifneq ($(MAKECMDGOALS),clean)
#-include $(depFiles)
#-include $(INTF:.ci=.di)
endif

