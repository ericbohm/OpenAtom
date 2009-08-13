# Define the base, build and include directories of the openatom tree
base   := .
makedir = $(base)/compile
include $(makedir)/srcdirs.mk

# Define the command line args to sub-makes
MAKEARGS = -C $(build) -f ../$(makedir)/Makefile


.PHONY: all driver physics libs clean clean_driver clean_physics clean_libs again test docs doxygen

all driver physics libs clean clean_driver clean_physics clean_libs again: $(build) 
	@$(MAKE) $(MAKEARGS) $@

realclean:
	@test ! -d $(build) || { $(MAKE) $(MAKEARGS) $@ && rmdir $(build); }

docs:
	@cd $(docs) && $(MAKE)

doxygen:
	@cd $(docs) && $(DOXYGEN) $(docs)/Doxyfile

$(build): 
	@mkdir -p $(build)

