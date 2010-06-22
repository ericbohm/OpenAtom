#get list of test directories
include tmp/*.name

.PHONY: tests $(testlist)

tests: $(testlist)
	@printf "ALL TESTS PASSED\n"

$(testlist):
	./pairCalcTest tmp/$@.config $(shell cat tmp/$@.dir)
	/bin/rm -f *.out
