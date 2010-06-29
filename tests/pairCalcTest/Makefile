# Makefile scripts to build and run pair calculator tests

CHARMDIR = $(HOME)/charm
CHARMC = $(CHARMDIR)/bin/charmc $(OPTS)
CC = g++
FCC = gfortran
FFT_HOME = $(HOME)/fftw

#Step one, build all necessary files
include build.mk

clean:
	/bin/rm -f ${BINARIES} a.out core *.o setup
	/bin/rm -f *.decl.h *.def.h charmrun
	/bin/rm -rf tmp/ STATES_OUT/
	/bin/rm -rf test_output

run: build runtests

.PHONY: runtests

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

