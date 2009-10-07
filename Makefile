# Define the base directory and get the openatom tree
base       := .
makedir     = $(base)/makefiles

# Include project directory hierarchy information
include $(makedir)/srcdirs.mk

# Determine if we'll need info from a configuration file
infoTargets:= help docs doxygen
shouldInclude := no
ifeq (,$(MAKECMDGOALS))
  shouldInclude := yes
endif
ifneq (,$(filter-out $(infoTargets),$(MAKECMDGOALS)))
  shouldInclude := yes
endif

# Include the configuration file if we're going to need it
ifeq (yes,$(shouldInclude))
  isExist  := $(shell test -r $(base)/config.mk || echo "no")
  ifeq (no,$(isExist))
    $(error $(base)/config.mk does not exist! Please copy $(makedir)/config.MACHINE.mk to $(base)/config.mk and customize it for your build)
  endif
  include $(base)/config.mk
endif

################## Build directory related stuff ####################
# Location in which the build directory should be created
where         = $(base)
# The name prefix of the build directory
builddir      = build
# The computed suffix of the build directory
buildsuffix   = $(shell echo $(sort $(OPT)) | sed "s| ||g")
# The actual build directory name and location
build         = $(strip $(where))/$(strip $(builddir))$(strip $(buildsuffix))

################## Test directory related stuff ####################
# The directory within which any tests are run
testop = $(build)/test-output
testop_regr = $(testop)/regression
testop_unit = $(testop)/unit

################## Arguments to the sub-makes ####################
# Define the command line args to sub-make
BUILDARGS     =-C "$(call realpath,$(build))" \
               -f "$(call realpath,$(makedir)/mainrecipe.mk)" \
               -r \
               base="$(call abs2rel,$(base),$(call realpath,$(build)))"

TESTARGS      =-C "$(call realpath,$1)" \
               -f "$(call realpath,$(makedir)/$2)" \
               -r \
               base="$(call abs2rel,$(base),$(call realpath,$1))" \
			   build="$(call abs2rel,$(build),$(call realpath,$1))"


.PHONY: all \
        compile driver physics libs \
		again clean clean_driver clean_physics clean_libs \
		test test-regr test-unit clean-test retest \
		docs doxygen

all: compile

################## Build-related targets ####################
compile driver physics libs again: $(build)
	@$(MAKE) $(BUILDARGS) $@
	@echo "=========== Build results are in the build directory: $(build)"

clean clean_driver clean_physics clean_libs: $(build)
	@$(MAKE) $(MAKEARGS) $@

realclean:
	@test ! -d $(testop) || $(RM) -r $(testop)
	@test ! -d $(build)  || $(RM) -r $(build)

################## Testing-related targets ####################
test: test-unit test-regr

test-unit:

test-regr: compile $(testop_regr)
	@$(MAKE) $(call TESTARGS,$(testop_regr),testrecipe.mk) $@
	@echo "=========== Regression test output is in the test directory: $(testop_regr)"

retest: clean-test test

clean-test:
	@test ! -d $(testop) || $(RM) -r $(testop)

################## Other targets ####################
docs:
	@cd $(docs) && $(MAKE)

doxygen:
	@cd $(docs) && doxygen $(call abs2rel,$(base),$(docs))/Doxyfile

################## Utility targets ####################
$(build) $(testop) $(testop_regr) $(testop_unit):
	@echo "=========== Creating directory: $@"
	@mkdir -p $@

help:
	@less $(makedir)/makehelp.txt

