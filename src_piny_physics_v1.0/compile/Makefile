# Define the base directory of the openatom tree
base := ../..
makedir := $(base)/compile

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(makedir)/Makefile.config
include $(makedir)/compileflags.mk

# Directory structure and vpath related stuff
# Overriding values set in srcdirs.mk for a shorter, cleaner paths
DEPSTRIPDIRS += 
VPATH        += $(STANDARD_INC) 

# Specify the compilation target
TARGET    = libPinyInterface.a

# The relevant source files for this project
include $(pinymake)/make_defs/obj_files.h
SRC       = 
HDR       = 
OBJ       = $(addsuffix .o, $(basename $(SRC)) ) $(LIB_OBJS)
INTF      = 
CPPFLAGS += -I$(FFT_HOME)/include -I$(CHARMBASE)/include/fftlib 


.PHONY: all clean again 

all: $(TARGET)

$(TARGET): $(OBJ) | $(LIB_DECLS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

clean:
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o *.a) $(TARGET) 

again: 
	$(MAKE) clean; $(MAKE)

#----------------- Temporary legacy variable definitions ---------------------------
include $(pinymake)/make_defs/proto_files.h
include $(pinymake)/targetdeps.mk

FOBJ          = $(FC) $(FFLAGS) $(CPPFLAGS) -c -o $@
COBJ          = $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@
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

