CHARMC=~/Programming/ppl/charm/bin/charmc

all: gw_bse.C gw_bse.ci gw_bse.h
	$(CHARMC) gw_bse.ci
	$(CHARMC) gw_bse.C
	$(CHARMC) gw_bse.o -o gw_bse

clean:
	rm -rf *.decl.h *.def.h *.o charmrun gw_bse
