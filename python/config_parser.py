#!/usr/bin/env python
def test(config_dict):
	import config_reader
	import command_generator
	import execute_command
	import summary_generator

	config_dict = config_reader.config_reader(config_dict)
	print 'Read test config......................................................'
	config_dict = command_generator.command_generator(config_dict)
	print 'Generated commands....................................................'
	execute_command.execute_command(config_dict)
	print 'Executed commands.....................................................'
	summary_generator.summary_generator(config_dict)
	print 'Testing complete......................................................'

import sys
import os
config_dict = {}

data_path = os.path.abspath(sys.argv[1])
exe_path = os.path.abspath(sys.argv[2])
test_path = os.path.join(data_path,'tests',sys.argv[3])

config_dict['data_path'] = data_path

config_dict['exe_path'] = exe_path
config_dict['charmrun'] = os.path.join(exe_path,'charmrun')
config_dict['exe'] = os.path.join(exe_path,'OpenAtom')

config_dict['test_path'] = test_path
config_dict['output_path'] = os.path.join(test_path,'output')
config_dict['config_file'] = os.path.join(test_path,'testConfig.yml')

config_dict['tests'] = []

print config_dict

test(config_dict)
