def config_parser(data, exe, test):
	import config_reader
	import command_generator
	import execute_command
	import summary_generator
	import os

	config_dict = {}

	data_path = os.path.abspath(data)
	exe_path = os.path.abspath(exe)
	test_path = os.path.join(data_path,'tests',test)

	config_dict['data_path'] = data_path

	config_dict['exe_path'] = exe_path
	config_dict['charmrun'] = os.path.join(exe_path,'charmrun')
	config_dict['exe'] = os.path.join(exe_path,'OpenAtom')

	config_dict['test_path'] = test_path
	config_dict['output_path'] = os.path.join(test_path,'output')
	config_dict['config_file'] = os.path.join(test_path,'testConfig.yml')

	config_dict['tests'] = []

	config_dict = config_reader.config_reader(config_dict)
	print 'Read ' + str(len(config_dict['tests'])) + ' tests from the config'
	config_dict = command_generator.command_generator(config_dict)
	print 'Executing tests...'
	execute_command.execute_command(config_dict)
	print 'Testing complete'
	summary = summary_generator.summary_generator(config_dict)
	return summary
