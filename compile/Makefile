# Define the base directory of the openatom tree
base := ..
makedir := $(base)/compile

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(makedir)/Makefile.config
include $(makedir)/compileflags.mk

# Directory structure and vpath related stuff
DEPSTRIPDIRS += 
VPATH        +=

# Specify the compilation target
TARGET    = OpenAtom

# The relevant source files for this project
SRC       = 
HDR       = 
OBJ       = $(addsuffix .o, $(basename $(SRC)) )
INTF      = 
LIBS      = $(moddriver) $(modphysics) $(modmath)


.PHONY: all driver physics libs mathlib clean again test translateInterface

all: $(TARGET)

$(TARGET): $(LIBS:%=lib%.a) $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

driver: $(libdriver)

physics: $(libphysics)

mathlib: $(libmath)

libs: mathlib

# Clean only the build artifacts. Leave the binary etc around
clean: 
	@$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o)

again: 
	@$(MAKE) clean; $(MAKE)


ifneq ($(MAKECMDGOALS),clean)
# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(SRC)) ) )
-include $(depFiles)
-include $(INTF:.ci=.di)
# Include the source lists and dependency information for each sub-module
include $(driver)/makerules.mk
include $(physics)/makerules.mk
include $(mathlib)/makerules.mk
endif

