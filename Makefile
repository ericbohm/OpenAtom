# Define the base directory and get the openatom tree
base   := .
makedir = $(base)/compile
include $(makedir)/srcdirs.mk
include $(makedir)/Makefile.config

################## Build directory related stuff ####################
# Location in which the build directory should be created
where         = $(base)
# The name prefix of the build directory
builddir      = build
# The computed suffix of the build directory
buildsuffix   = $(shell echo $(sort $(OPT)) | sed "s| ||g")
# The actual build directory name and location
build         = $(strip $(where))/$(strip $(builddir))$(strip $(buildsuffix))

# Define the command line args to sub-make
MAKEARGS      =-C $(call realpath,$(build)) \
               -f $(call realpath,$(makedir)/Makefile) \
               -r \
               base=$(call abs2rel,$(base),$(call realpath,$(build)))

.PHONY: all driver physics libs clean clean_driver clean_physics clean_libs again test docs doxygen

all driver physics libs again: $(build)
	@$(MAKE) $(MAKEARGS) $@
	@echo "=========== Build results are in the build directory: $(build)"

clean clean_driver clean_physics clean_libs: $(build)
	@$(MAKE) $(MAKEARGS) $@

realclean:
	@test ! -d $(build) || $(RM) -r $(build)

docs:
	@cd $(docs) && $(MAKE)

doxygen:
	@cd $(docs) && $(DOXYGEN) $(docs)/Doxyfile

$(build):
	@echo "=========== Creating build directory: $@"
	@mkdir -p $@

help:
	@less $(makedir)/makehelp.txt

