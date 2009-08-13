# Define the base directory of the openatom tree
base := ..
makedir := $(base)/compile

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(makedir)/Makefile.config

# Directory structure and vpath related stuff
CPAIMD        = $(base)
PINY_STD_INC  = $(STANDARD_INC)
DEPSTRIPDIRS += 
VPATH        += $(build_drv) $(build_phy) $(build_math) 

# Specify the compilation target
TARGET    = OpenAtom

# The relevant source files for this project
SRC       = 
HDR       = 
OBJ       = $(addsuffix .o, $(basename $(SRC)) )
INTF      = 
LIBS      = CharmDriver PinyInterface MyMathLib
LIBNAMES  = $(LIBS:%=lib%.a)

########### Specify the flags, include paths, libraries etc.###########
# Preprocessor flags (Used for all preprocessing, compilation and linking)
CPPFLAGS += -I$(base) -I$(base)/include 
# Compiler flags (Used for compilation and linking)
CFLAGS   += 
CXXFLAGS += 
FFLAGS   += 
# Linker flags (Used for all linking)
LDFLAGS  += -L$(FFT_HOME)/lib -memory gnu 
# Additional libraries
LDLIBS   += -module CkMulticast -module comlib
ifeq ($(DUAL_FFTW), -DDUAL_FFTW_OFF)
  LDLIBS += -lrfftw -lfftw
else
  LDLIBS += -ldrfftw -ldfftw
endif
LDLIBS   += -lz -lm -lconv-util


.PHONY: all driver physics libs clean clean_driver clean_physics clean_libs again test translateInterface

all: $(TARGET)

$(TARGET): $(LIBS:%=-l%) $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(build_drv)/libCharmDriver.a: driver

$(build_phy)/libPinyInterface.a: physics

$(build_math)/libMyMathLib.a: libs

driver:
	@cd $(build_drv) && $(MAKE) CXXFLAGS+="$(OPTS)"

physics:
	@cd $(build_phy) && $(MAKE)

libs:
	@cd $(build_math) && $(MAKE)

# Clean only the build artifacts. Leave the binary etc around
clean: clean_driver clean_physics clean_libs
	@$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o) 

realclean: clean
	@$(RM) $(LIBNAMES) $(TARGET) charmrun 

clean_driver: 
	@cd $(build_drv) && $(MAKE) clean

clean_physics: 
	@cd $(build_phy) && $(MAKE) clean

clean_libs:
	@cd $(build_math) && $(MAKE) clean

again: 
	@$(MAKE) clean; $(MAKE)


# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(SRC)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(INTF:.ci=.di)
endif

