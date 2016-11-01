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
ifeq ($(DIAGONALIZER), YES)
  LIBS     += $(modinterop)
endif

# Compute the revision number (hash) of the build and feed it to the code
ifeq ($(GIT),)
  $(warning Cannot find the git binary. Will not compute the revision number)
else
  REVNUM  = $(shell $(GIT) --git-dir=$(base)/.git rev-parse HEAD)
endif
CPPFLAGS += -DOPENATOM_REVISION=$(REVNUM)

.PHONY: compile driver physics libs mathlib interoplib clean again test translateInterface

compile: $(TARGET)

$(TARGET): $(LIBS:%=lib%.a) $(OBJ)
	$(info-ld)
ifeq ($(DIAGONALIZER), YES)
	ar q libmodulecpaimd.a *.o
	$q$(CXX) -mpi -nomain-module $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(LDLIBS)
else
	$q$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
endif
	@echo "=========== Target produced from commit hash: $(REVNUM)"

driver: $(libdriver)

physics: $(libphysics)

mathlib: $(libmath)

interoplib: $(libinterop)

libs: mathlib

# Clean only the build artifacts. Leave the binary etc around
clean:
	@$(RM) $(wildcard *.decl.h *.def.h *.o *.ci.stamp)

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
include $(interoplib)/makerules.mk
endif

