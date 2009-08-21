# Define the base directory and get the openatom tree
base   := .
makedir = $(base)/compile
include $(makedir)/srcdirs.mk

################## Build directory related stuff ####################
# Location in which the build directory should be created
where         = $(base)
# The name prefix of the build directory
builddir      = build
# The actual build directory name and location
build         = $(strip $(where)/$(builddir)$(buildsuffix))

# Define the command line args to sub-makes
MAKEARGS      = -C $(build) -f $(abspath $(makedir)/Makefile) 
# No available, assured portable utility to compute the relative path between any two absolute paths
# Hence use absolute path to specify the project root if build directory is elsewhere
ifneq ($(abspath $(where)),$(abspath $(base)))
MAKEARGS     += base=$(abspath $(CURDIR))
endif

.PHONY: all driver physics libs clean clean_driver clean_physics clean_libs again test docs doxygen

all driver physics libs clean clean_driver clean_physics clean_libs again: $(build)
	@$(MAKE) $(MAKEARGS) $@

realclean:
	@test ! -d $(build) || $(RM) -r $(build)

docs:
	@cd $(docs) && $(MAKE)

doxygen:
	@cd $(docs) && $(DOXYGEN) $(docs)/Doxyfile

$(build):
	@echo "Creating build directory $@"
	@mkdir -p $@

