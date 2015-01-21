def execute_command(config_dict):
	import os
	import shutil
	import subprocess

	os.chdir(config_dict['data_path'])
	shutil.rmtree(config_dict['output_path'])
	os.mkdir(config_dict['output_path'])
	tests = config_dict['tests']

	for test_dict in tests:
		setup_dataset()
		clean_dataset()

		command = test_dict['command']
		print 'Running command: ' + `command`
		outfile = open(test_dict['output_file'],'w')
		subprocess.call(command, stdout = outfile, stderr = outfile)
		outfile.close()

def clean_dataset():
	import subprocess
	print 'Cleaning dataset'
	subprocess.call('../../utils/tidy')

def setup_dataset():
	import subprocess
	print 'Setting up dataset'
	subprocess.call('../../utils/setup')
