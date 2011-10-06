#!/usr/bin/env python
def test(filename, executable_path):
	import config_reader
	import command_generator
	import folder_generator
	import execute_command
	import summary_generator
	import os
	current_working_dir = os.getcwd()
	mydict = config_reader.config_reader(filename, executable_path)
	mydict = command_generator.command_generator(mydict)
	folder_generator.folder_generator(mydict)
	execute_command.execute_command(mydict, current_working_dir)
	os.chdir(current_working_dir)
	mydict_string = repr(mydict)
	dict_file = open('info_dict','w+')
	dict_file.write(mydict_string)
	dict_file.close()
	summary_generator.summary_generator('info_dict')
import sys
configfilename = sys.argv[1]
executable_path = sys.argv[2]
test(configfilename, executable_path)
