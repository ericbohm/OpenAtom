def execute_command(config_dict):
	import os
	import shutil
	import subprocess

	os.chdir(config_dict['data_path'])

	if os.path.exists(config_dict['output_path']):
		shutil.rmtree(config_dict['output_path'])
	os.mkdir(config_dict['output_path'])

	tests = config_dict['tests']
	for i in range(len(tests)):
		test_dict = tests[i]
		setup_dataset()
		clean_dataset()

		command = test_dict['command']
		outfile = open(test_dict['output_file'],'w')
		subprocess.call(command, stdout = outfile, stderr = outfile)
		outfile.close()
		print 'Executed ' + str(i+1) + ' out of ' + str(len(tests)) + ' tests'

def clean_dataset():
	import os
	import subprocess
	devnull = open(os.devnull, 'w')
	subprocess.call('../../utils/tidy', stdout=devnull, stderr=devnull)

def setup_dataset():
	import os
	import subprocess
	devnull = open(os.devnull, 'w')
	subprocess.call('../../utils/setup', stdout=devnull, stderr=devnull)
