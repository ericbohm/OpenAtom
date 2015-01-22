def summary_generator(config_dict):
	import yaml
	import checkresult
	import string

	passed = 0
	failed = 0
	test_strings = []

	tests = config_dict['tests']

	header = string.ljust('TEST NAME',20)
	header = header + string.ljust('VARIANT',20)
	header = header + string.ljust('PES',4) + 'RESULT'

	prev_name = ''

	iteration = config_dict['iteration']
	sig_figs = config_dict['sig_figs']

	for test_dict in tests:
		test_string = ''
		if test_dict['test_name'] == prev_name:
			test_string = test_string + string.ljust('',20)
		else:
			prev_name = test_dict['test_name']
			test_string = test_string + string.ljust(test_dict['test_name'],20)
			if len(test_strings) > 0:
				test_strings.append('')

		test_string = test_string + string.ljust(test_dict['variant'],20)
		test_string = test_string + string.ljust(test_dict['numpe'],4)

		output_file = test_dict['output_file']
		out_ref = test_dict['out_ref']

		test_result = checkresult.checkresult(output_file, out_ref, iteration, sig_figs)
		if test_result[0] == True:
			test_strings.append(test_string + 'PASSED')
			passed = passed + 1
		if test_result[0] == False:
			test_strings.append(test_string + 'FAILED: ' + test_result[1])
			failed = failed + 1

	print '\n'
	print string.ljust('',80,'=')
	print header
	print string.ljust('',80,'-')
	for value in test_strings:
		print value
	print string.ljust('',80,'=')
	print `passed` + ' PASSED\t' + `failed` + ' FAILED'
	print '\n'
