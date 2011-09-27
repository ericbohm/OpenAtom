def test(filename, set_name):
	import config_reader
	import command_generator
	import folder_generator
	import execute_command
	import summary_generator
	mydict = config_reader.config_reader(filename, set_name)
	mydict = command_generator.command_generator(mydict)
	folder_generator.folder_generator(mydict)
	execute_command.execute_command(mydict)
	mydict_string = repr(mydict)
	dict_file = open('info_dict','w+')
	dict_file.write(mydict_string)
	dict_file.close()
	summary_generator.summary_generator('info_dict')
import sys
configfilename = sys.argv[1]
set_name = sys.argv[2]
test(configfilename, set_name)
