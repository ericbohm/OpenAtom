# Define the base directory of the openatom tree
base := ..

# Include basic make rules
makedir   = $(base)/compile
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(makedir)/Makefile.config

# Directory structure and vpath related stuff
# Overriding values set in srcdirs.mk for a shorter, cleaner paths
build         = $(base)/binary
CPAIMD        = $(base)
PINY_STD_INC  = $(STANDARD_INC)
DEPSTRIPDIRS += 
VPATH        += 

# Specify the compilation target
TARGET    = OpenAtom

# The relevant source files for this project
SRC       = 
HDR       = 
OBJ       = $(SRC:.C=.o)
INTF      = 
LIBS      = CharmDriver PinyInterface MyMathLib
LIBNAMES  = $(LIBS:%=lib%.a)

########### Specify the flags, include paths, libraries etc.###########
# Preprocessor flags (Used for all preprocessing, compilation and linking)
CPPFLAGS += -I$(driver) -I$(build) -I$(base) -I$(base)/include -I$(PINY_STD_INC)
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

$(TARGET): $(LIBNAMES) $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

libCharmDriver.a:
	cd $(driver)/compile && $(MAKE) CXXFLAGS+="$(OPTS)"

libPinyInterface.a:
	cd $(MAKE_PHYS_INC) && $(MAKE) PINYLIB

libMyMathlib.a:
	cd $(mathlib) && $(MAKE) OPTS="$(OPTS)"

driver: libCharmDriver.a

physics: libPinyPhysics.a

libs: libMyMathlib.a

# Clean only the build artifacts. Leave the binary etc around
clean: clean_driver clean_physics clean_libs
	$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o) 

clean_driver: 
	cd $(driver)/compile && $(MAKE) clean

clean_physics: 
	cd $(MAKE_PHYS_INC) && $(MAKE) clean

clean_libs:
	cd $(mathlib) && $(MAKE) clean

again: 
	$(MAKE) clean; $(MAKE)


# Include the generated files containing dependency info
-include $(SRC:.C=.d)
-include $(INTF:.ci=.di)

