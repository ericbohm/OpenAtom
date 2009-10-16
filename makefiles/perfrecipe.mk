# Define the base directory of the openatom tree
base := ../../..
makedir := $(base)/makefiles
cfg  := $(base)/config.mk

# Include basic make rules
include $(makedir)/srcdirs.mk
include $(makedir)/basicrules.mk
include $(cfg)

# Compute the list of run sizes (no of cores)
procList = $(shell myvar=$(procStart); while [ $$myvar -le $(procEnd) ];do echo $$myvar; myvar=$$(( $$myvar*2 ));done )

# Generate the list of targets that will be used to setup the jobs
setupTargets = $(procList:%=setup_%)
# Generate the list of targets that will be used to submit the jobs
submitTargets= $(procList:%=submit_%)

# This is the unprocessed command line that will be used to submit the job
cmd_raw = $(submitLine) $(perfInpDir)/$(perfConfig) $(perfInpDir)/$(perfInput) $(perfOtherArgs)
# An innocent placeholder-replacement mechanism to process the command line
parseIt = $(subst @X,$(build)/$(executable),$(subst @T,$(walltime),$(subst @N,$3,$(subst @C,$2,$1))))
# Run some arithmetic in bash to compute the number of nodes
numNodes = $(shell echo $$(( $1/$(ppn) )))

$(submitTargets): override base  := $(addprefix ../,$(base))
$(submitTargets): override build := $(addprefix ../,$(build))

.PHONY: perf perf-setup

perf: perf-setup perf-submit

perf-setup: $(setupTargets)

perf-submit: $(submitTargets)

submit_%: setup_% | proc_%
	@echo "Submitting the job using command ..."
	cd $| && $(call parseIt,$(cmd_raw),$*,$(call numNodes,$*))

setup_%: | proc_%
	@echo "Setting up a job using $* cores on $(call numNodes,$*) nodes in directory: $|"
	@$(LN) $(call realpath,$(perfDatabase)) $|/..
	@$(LN) $(wildcard $(call realpath,$(perfDataset))/*) $|
	@touch $@

proc_%:
	@echo "=========== Creating directory: $@"
	@mkdir -p $@/STATES_OUT

# Prevent make from attempting to delete the directories which it considers to be intermediate cruft
.PRECIOUS: proc_%

cleanRunDebris  = $(RM) *.out $1.confv $1.iavg $1.confp $1.params $1.confc $1.coords_out $1

