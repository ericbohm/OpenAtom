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
ees-nl0l0-p1-args = $(testregr)/cpaimd_config.p1 $(testregr)/water.input.min.ees-nl0l0
ees-nl0l1-p1-desc = EES: off for nonlocals; on for locals;
ees-nl0l1-p1-args = $(testregr)/cpaimd_config.p1 $(testregr)/water.input.min.ees-nl0l1
ees-nl1l0-p1-desc = EES: on for nonlocals; off for locals;
ees-nl1l0-p1-args = $(testregr)/cpaimd_config.p1 $(testregr)/water.input.min.ees-nl1l0
ees-nl1l1-p1-desc = EES: on for nonlocal; on for locals;
ees-nl1l1-p1-args = $(testregr)/cpaimd_config.p1 $(testregr)/water.input.min.ees-nl1l1

# Define the list of tests to run
testList   = ees-nl0l0-p1 ees-nl0l1-p1 ees-nl1l0-p1 ees-nl1l1-p1


.PHONY: test-regr

test-regr: $(testList:%=op-%.log)

op-%-p1.log: setup $(binary)
	@printf "Running regression test $*: $($*-p1-desc) ..."
	@$(call cleanRunDebris,water)
	@$(binary) $($*-p1-args) 2>&1 > $@
	@sed -ne "/Iteration 1 done/,/Iteration 2 done/p" $@ | grep "\WPsi" | sort > snip-$@
	@cat $(testregr)/ref-$*-p1.log | grep "\WPsi" | sort > ref-$@
	@$(testregr)/sigcmp.pl 9 ref-$@ snip-$@ && echo "\t\t Passed" || echo "\t\t\t TEST FAILED!"

setup:
	@$(LN) $(wildcard $(data)/water_32M_10Ry/*) .
	@$(LN) $(realpath $(data)/DATABASE) ..
	@mkdir -p STATES_OUT
	@touch $@

cleanRunDebris  = $(RM) *.out $1.confv $1.iavg $1.confp $1.params $1.confc $1.coords_out $1

