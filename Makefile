# Define the base directory and default config file to read in
base       := .
makedir     = $(base)/makefiles
cfg        := $(base)/config.mk

# Include project directory hierarchy information
include $(makedir)/srcdirs.mk

################## Build directory related variables that can reset in the config file or cmd line ####################
# Location in which the build directory should be created
where         = $(base)
# The name prefix of the build directory
builddir     ?= build

################## Determine if we'll need info from a configuration file ##################
infoTargets:= help docs doxygen
shouldInclude := no
ifeq (,$(MAKECMDGOALS))
  shouldInclude := yes
endif
ifneq (,$(filter-out $(infoTargets),$(MAKECMDGOALS)))
  shouldInclude := yes
endif

################## Include the configuration file if we're going to need it ##################
ifeq (yes,$(shouldInclude))
  isExist  := $(shell test -r $(cfg) || echo "no")
  ifeq (no,$(isExist))
    $(error $(cfg) does not exist! Please copy $(makedir)/config.MACHINE.mk to $(cfg) and customize it for your build)
  endif
  include $(cfg)
endif

################## Compute the build directory name ####################
# The computed suffix of the build directory
buildsuffix   = $(shell echo $(sort $(OPT)) | sed "s| ||g")
# The actual build directory name and location
build         = $(strip $(where))/$(strip $(builddir))$(strip $(buildsuffix))

################## Test directory related stuff ####################
# The directory within which any tests are run
testop = $(build)/test-output
testop_regr = $(testop)/regression
testop_unit = $(testop)/unit
testop_perf = $(testop)/scaling

################## Arguments to the sub-makes ####################
# Define the command line args to sub-make
BUILDARGS     =-C "$(call realpath,$(build))" \
               -f "$(call realpath,$(makedir)/mainrecipe.mk)" \
               -r \
               base="$(call abs2rel,$(base),$(call realpath,$(build)))" \
               cfg="$(call abs2rel,$(cfg),$(call realpath,$(build)))"

TESTARGS      =-C "$(call realpath,$1)" \
               -f "$(call realpath,$(makedir)/$2)" \
               -r \
               base="$(call abs2rel,$(base),$(call realpath,$1))" \
               cfg="$(call abs2rel,$(cfg),$(call realpath,$1))" \
               build="$(call abs2rel,$(build),$(call realpath,$1))"


.PHONY: all \
        compile driver physics libs \
		again clean clean_driver clean_physics clean_libs \
		test test-regr test-unit clean_test retest \
		perf perf-setup perf-chk clean_perf
		docs doxygen

all: compile

################## Build-related targets ####################
compile driver physics libs again: $(build)
	@$(MAKE) $(BUILDARGS) $@
	@echo "=========== Build results are in the build directory: $(build)"

clean clean_driver clean_physics clean_libs: $(build)
	@$(MAKE) $(BUILDARGS) $@

realclean:
	@test ! -d $(testop) || $(RM) -r $(testop)
	@test ! -d $(build)  || $(RM) -r $(build)

################## Testing-related targets ####################
test: test-unit test-regr

test-unit:

test-regr: compile
	cd python; python test_runner.py tests.config

retest: clean_test test

clean_test:
	@test ! -d $(testop) || $(RM) -r $(testop)

################## Scaling / Performance measurement related targets ####################
perf perf-setup perf-chk: compile $(testop_perf)
	@$(MAKE) $(call TESTARGS,$(testop_perf),perfrecipe.mk) $@
	@echo "=========== Performance measurement output is in the scaling directory: $(testop_perf)"

clean_perf:
	@test ! -d $(testop) || $(RM) -r $(testop_perf)

################## Other targets ####################
docs:
	@cd $(docs) && $(MAKE)

doxygen:
	@cd $(docs) && doxygen Doxyfile

clean_docs:
	@cd $(docs) && $(MAKE) clean
	@$(RM) -r $(docs)/html

################## Utility targets ####################
$(build) $(testop) $(testop_regr) $(testop_unit) $(testop_perf):
	@echo "=========== Creating directory: $@"
	@mkdir -p $@

help:
	@less $(makedir)/makehelp.txt

