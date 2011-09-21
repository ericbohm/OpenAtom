def execute_command(config_dict):
	import os
	folder_list = config_dict.keys()
	folder_list.remove('general_info')
	for folder_name in folder_list:
		os.chdir('../build-O3')
		os.chdir(config_dict[folder_name][13])
		os.chdir(folder_name)
		os.popen(config_dict[folder_name][15])
		os.chdir('../../../python')
