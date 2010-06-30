# Makefile scripts to build and run pair calculator tests

BINARIES = pairCalcTest parseTestConfig

.PHONY: build clean run runtests setup

build: ${BINARIES}

parseTestConfig: parseTestConfig.C
	g++ -o $@ $^

clean:
	@$(RM) $(wildcard *.decl.h *.def.h *.o *.a)
	@$(RM) $(wildcard *.d *.di *.ci.stamp)
	@$(RM) ${BINARIES} charmrun a.out core setup
	@$(RM) -r tmp/ STATES_OUT/ test_output/

run: build runtests

runtests: setup pairCalcTest
	mkdir test_output
	/bin/rm -f *.out
	$(MAKE) -f tests.mk

setup: testConfig.txt parseTestConfig
	mkdir -p STATES_OUT
	mkdir -p tmp
	/bin/rm -f tmp/*
	./parseTestConfig testConfig.txt
	touch $@

# Include build script for the test code
include build.mk

