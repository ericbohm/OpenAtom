# Define the base, build and include directories of the openatom tree
base   := .
makedir = $(base)/compile
include $(makedir)/srcdirs.mk
build   = $(base)/binary

.PHONY: all driver physics libs clean clean_driver clean_physics clean_libs again test docs doxygen

all: $(build)
	cd $(build) && $(MAKE) -f ../$(makedir)/Makefile $@

driver physics libs clean clean_driver clean_physics clean_libs again: $(build) 
	cd $(build) && $(MAKE) -f ../$(makedir)/Makefile $@

realclean:
	cd $(build) && $(MAKE) -f ../$(makedir)/Makefile $@
	rmdir $(build)

docs:
	cd $(docs) && $(MAKE)

doxygen:
	cd $(docs) && $(DOXYGEN) $(docs)/Doxyfile

$(build): 
	mkdir -p $(build)

