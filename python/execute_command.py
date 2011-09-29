def execute_command(config_dict, original_working_dir):
	import os
	import stat
	current_working_dir = os.getcwd()
	os.chdir(original_working_dir)
	testscript = open('testscript', 'a')
	testscript.close()
	os.chmod('testscript', 511)
	testscript = open('testscript', 'a')
	os.chdir(current_working_dir)
	folder_list = config_dict.keys()
	folder_list.remove('general_info')
	for folder_name in folder_list:
		os.chdir(folder_name)
		testscript.write('cd '+folder_name+'\n')

		os.popen(config_dict[folder_name]['command'])
		testscript.write(config_dict[folder_name]['command']+'\n')

		os.chdir('..')
		testscript.write('cd ..\n')
	testscript.close()

