CHARMC=~/Programming/ppl/charm/bin/charmc

OBJS=gw_bse.o psi.o fcalculator.o

all: gw_bse.decl.h $(OBJS)
	$(CHARMC) -o gw_bse $(OBJS) -module CkMulticast -language charm++

gw_bse.decl.h: gw_bse.ci config.h
	$(CHARMC) gw_bse.ci

%.o: %.C %.h
	$(CHARMC) $<

clean:
	rm -rf *.decl.h *.def.h *.o charmrun gw_bse
