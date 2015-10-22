#CHARMDIR = /Users/minjung/packages/Charm
CHARMDIR = /scratch/mikida2/charm
CHARMC = $(CHARMDIR)/bin/charmc
FFT3 = /dcsdata/home/bohm/fftw3

INCLUDES=-I$(FFT3)/include
LIBS=-L$(FFT3)/lib -lfftw3 -lm -lz
CFLAGS=-g $(INCLUDE) #-D_CP_DEBUG_UTIL_VERBOSE_ -DDEBUG_FFT
LFLAGS=$(LIBS)

C_FILES=$(wildcard *.C)
CI_FILES=$(wildcard *.ci)
MODULES=$(CI_FILES:.ci=.decl.h)
OBJS=$(C_FILES:.C=.o)

default: all
all: gw_bse

gw_bse: $(MODULES) $(OBJS)
	$(CHARMC) $(LFLAGS) $(OBJS) -o $@

%.decl.h: %.ci
	$(CHARMC) $<

%.o: %.C %.h $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.decl.h *.def.h
	rm -rf charmrun gw_bse
