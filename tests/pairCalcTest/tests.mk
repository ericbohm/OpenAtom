#get list of test directories
include tmp/*.name

.PHONY: tests $(testlist)

tests: $(testlist)
	@printf "TESTS FINISHED - detailed results in the tests_output directory\n"
	@printf "===============================================================\n"
	@printf "|                   Results Summary:                          |\n"
	@printf "===============================================================\n"
	@tail -n 1 test_output/*

$(testlist):
	./pairCalcTest tmp/$@.config $(shell cat tmp/$@.dir) 2>test_output/$@.out
