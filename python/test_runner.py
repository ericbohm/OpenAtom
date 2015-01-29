def execute_tests(tests):
	import os
	import shutil
	import subprocess

	import test_checker
	import test_output
	import test_util

	print '\nExecuting ' + `len(tests)` + ' tests...\n'

	# First just go through and make sure all output dirs are created fresh
	for test_dict in tests:
		if os.path.exists(test_dict['output_path']):
			shutil.rmtree(test_dict['output_path'])
		os.mkdir(test_dict['output_path'])

	test_output.print_table_header()

	passed = 0
	failed = 0
	# Then go through and execute each test one at a time
	for i in range(len(tests)):
		test_dict = tests[i]

		os.chdir(test_dict['data_path'])
		test_util.setup_dataset()
		test_util.clean_dataset()

		command = generate_command(test_dict)

		test_output.print_summary_string(i+1, len(tests), test_dict)

		f = open(test_dict['output_file'],'w')
		subprocess.call(command, stdout = f, stderr = f)
		f.close()

		test_dict['result'] = test_checker.check_result(test_dict)

		test_output.print_result_string(test_dict)

		if test_dict['result'][0]:
			passed = passed + 1
		else:
			failed = failed + 1

	test_output.print_table_footer(passed, failed)

	return failed

def generate_command(test_dict):
	import os

	data_path = test_dict['data_path']
	charmrun = os.path.relpath(test_dict['charmrun'],data_path)
	exe = os.path.relpath(test_dict['exe'],data_path)

	numpe = test_dict['numpe']
	par_cfg = os.path.relpath(test_dict['par_cfg'],data_path)
	phy_cfg = os.path.relpath(test_dict['phy_cfg'],data_path)

	command = [charmrun, '++local', '+p'+numpe, exe, par_cfg, phy_cfg]

	return command
