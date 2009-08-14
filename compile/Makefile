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
VPATH        += $(build_driver) $(build_physics) $(build_mathlib) 

# Specify the compilation target
TARGET    = OpenAtom

# The relevant source files for this project
SRC       = 
HDR       = 
OBJ       = $(addsuffix .o, $(basename $(SRC)) )
INTF      = 
LIBS      = $(moddriver) $(modphysics) $(modmath)

# Tiny function to help generate the correct arguments for a module's sub-make invocation
# @note: Some string manipulation to automatically generate correct relative paths will be good
genMakeArgs  = -C $(build_$1) -f ../$($1)/compile/Makefile base=../..


.PHONY: all driver physics libs mathlib clean again test translateInterface

all: $(TARGET)

$(TARGET): $(LIBS:%=-l%) $(OBJ) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(build_driver)/libCharmDriver.a: driver

$(build_physics)/libPinyInterface.a: physics

$(build_mathlib)/libMyMathLib.a: libs

driver physics mathlib:
	@$(MAKE) $(call genMakeArgs,$@)

libs: mathlib

# Clean only the build artifacts. Leave the binary etc around
clean: 
	@$(RM) -r $(build_driver) $(build_physics) $(build_mathlib)
	@$(RM) $(wildcard *.decl.h *.def.h *.d *.di *.o)

again: 
	@$(MAKE) clean; $(MAKE)


# Include the generated files containing dependency info
depFiles := $(addsuffix .d, $(basename $(filter %.C %.cpp %.cxx %.c, $(SRC)) ) )
ifneq ($(MAKECMDGOALS),clean)
-include $(depFiles)
-include $(INTF:.ci=.di)
endif

