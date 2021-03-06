# Define the base directory of the openatom tree
base := ..
makedir := $(base)/makefiles
cfg  := $(base)/config.mk

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(cfg)

# Basic test command settings
runCmd     = $(build)/charmrun
binary     = $(build)/$(executable)
opFile     = op-$*.log
refFile    = ref-$*.log

# Define all the test configurations
ees-nl0l0-p1-desc = EES: off for nonlocals; off for locals;
ees-nl0l0-p1-args = $(regrInpDir)/cpaimd_config.p1 $(regrInpDir)/water.input.min.ees-nl0l0 +logsize 10000000
ees-nl0l1-p1-desc = EES: off for nonlocals; on for locals;
ees-nl0l1-p1-args = $(regrInpDir)/cpaimd_config.p1 $(regrInpDir)/water.input.min.ees-nl0l1 +logsize 10000000 
ees-nl1l0-p1-desc = EES: on for nonlocals; off for locals;
ees-nl1l0-p1-args = $(regrInpDir)/cpaimd_config.p1 $(regrInpDir)/water.input.min.ees-nl1l0 +logsize 10000000
ees-nl1l1-p1-desc = EES: on for nonlocal; on for locals;
ees-nl1l1-p1-args = $(regrInpDir)/cpaimd_config.p1 $(regrInpDir)/water.input.min.ees-nl1l1 +logsize 10000000

# Define the list of tests to run
testList   = ees-nl0l0-p1 ees-nl0l1-p1 ees-nl1l0-p1 ees-nl1l1-p1


.PHONY: test-regr

test-regr: $(testList:%=op-%.log)

op-%-p1.log: setup $(binary)
	@printf "Running regression test $*: $($*-p1-desc) ..."
	@./tidy water >/dev/null 2>&1 || true
	@$(binary) $($*-p1-args) 2>&1 > $@
	@grep "Iter .1." $@ | cut -d' ' -f3-10| grep "Psi\[" | sort  > snip-$@
	@cat $(regrInpDir)/ref-$*-p1.log | grep "Psi\[" | sort > ref-$@
	@$(w3210)/sigcmp.pl 9 ref-$@ snip-$@ 5 && echo "\t\t Passed" || echo "\t\t\t TEST FAILED!"

setup:
	@$(LN) $(wildcard $(w3210)/*) .
	@$(LN) $(realpath $(molDbase)) ..
	@../../../utils/setup
	@touch $@
