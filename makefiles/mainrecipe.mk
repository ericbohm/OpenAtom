# Define the base directory of the openatom tree
base := ..
makedir := $(base)/makefiles
cfg  := $(base)/config.mk

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(cfg)

# Directory structure and vpath related stuff
DEPSTRIPDIRS +=
VPATH        +=

# Specify the compilation target
TARGET    = $(executable)

# The relevant source files for this project
SRC       =
HDR       =
OBJ       = $(addsuffix .o, $(basename $(SRC)) )
INTF      =
LIBS      = $(moddriver) $(modphysics) $(modmath)


.PHONY: compile driver physics libs mathlib clean again test translateInterface

compile: $(TARGET)

$(TARGET): $(LIBS:%=lib%.a) $(OBJ)
	$(info-ld)
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

driver: $(libdriver)

physics: $(libphysics)

mathlib: $(libmath)

libs: mathlib

# Clean only the build artifacts. Leave the binary etc around
clean:
	@$(RM) $(wildcard *.decl.h *.def.h *.o)

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

