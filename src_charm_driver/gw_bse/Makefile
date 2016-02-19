include config.mk
CHARMC = $(CHARMDIR)/bin/charmc

INCLUDES=-I$(FFTW3)/include
LIBS=-L$(FFTW3)/lib -lfftw3 -lm -lz -module CkLoop

CFLAGS+=-g $(INCLUDES)
LFLAGS+=$(LIBS) $(LAPACKLIB)

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

%.o: %.C $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o *.decl.h *.def.h
	rm -rf charmrun gw_bse
