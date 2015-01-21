def summary_generator(config_dict):
	import yaml
	import checkresult
	import string

	passed = 0
	failed = 0
	test_strings = []

	tests = config_dict['tests']

	header = string.ljust('TEST NAME',16) + '\t'
	header = header + string.ljust('VARIANT',16) + '\t'
	header = header + string.ljust('PES',4) + '\t' + 'RESULT'

	prev_name = ''

	for test_dict in tests:
		test_string = ''
		if test_dict['test_name'] == prev_name:
			test_string = test_string + string.ljust('',16) + '\t'
		else:
			prev_name = test_dict['test_name']
			test_string = test_string + string.ljust(test_dict['test_name'],16) + '\t'
			if len(test_strings) > 0:
				test_strings.append('')

		test_string = test_string + string.ljust(test_dict['variant'],16) + '\t'
		test_string = test_string + string.ljust(test_dict['numpe'],4)

		output_file = test_dict['output_file']
		out_ref = test_dict['out_ref']

		test_result = checkresult.checkresult(output_file, out_ref, 1, 5)
		if test_result == True:
			test_strings.append(test_string + '\tPASSED')
			passed = passed + 1
		if test_result == False:
			test_strings.append(test_string + '\tFAILED')
			failed = failed + 1

	print '\n'
	print string.ljust('',64,'=')
	print header
	print string.ljust('',64,'-')
	for value in test_strings:
		print value
	print string.ljust('',64,'=')
	print `passed` + ' PASSED\t' + `failed` + ' FAILED'
	print '\n'
