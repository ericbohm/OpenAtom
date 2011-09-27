def folder_generator(config_dict):
	import os
	folder_list = config_dict.keys()
	folder_list.remove('general_info')
# [0] = key, [1] = testname, [2] = machinename, [3] = inputDir, [4] = inputPar, [5] = inputPhy, [6] = outRef, [7] = testoutput, [8] = numpe, [9] = prefix_to execute, [10] = prefix_to_data, [11] = prefix_to_top, [12] = set_name, [13] = output_folder_name, [14] = prefix_to_output_folder, [15] = command, 
	for folder_name in folder_list:
		to_execute = config_dict[folder_name][9]
		os.chdir('../build-O3')
		command = 'mkdir ' + config_dict[folder_name][13]
		os.popen(command)
		os.chdir(config_dict[folder_name][13])
		command = 'mkdir ' + folder_name
		print command
		os.popen(command)
		os.chdir(folder_name)
		command = 'ln -s ../../../' + config_dict[folder_name][3] + '* .'
		os.popen(command)
		command = 'ln -s ../../data/DATABASE ..'
		os.popen(command)
		command = 'sh ../../../utils/setup'
		os.popen(command)
		command = 'sh ./tidy '+config_dict[folder_name][12]
		os.popen(command)
		os.chdir('..')
		os.chdir('..')
		os.chdir('..')
		os.chdir('python')
