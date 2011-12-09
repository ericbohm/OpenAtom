def folder_generator(config_dict):
	import os
	testscript = open('testscript', 'a')
	folder_list = config_dict.keys()
	folder_list.remove('general_info')
	executable_path = config_dict[folder_list[0]]['executable_path']
	os.chdir(executable_path)
	
	output_folder = config_dict[folder_list[0]]['output_folder']
	command = 'mkdir ' + output_folder
	os.popen(command)
	testscript.write(command + '\n')

	os.chdir(output_folder)
	testscript.write('cd ' + output_folder + '\n')
# [0] = key, [1] = testname, [2] = machinename, [3] = inputDir, [4] = inputPar, [5] = inputPhy, [6] = outRef, [7] = testoutput, [8] = numpe, [9] = prefix_to execute, [10] = prefix_to_data, [11] = prefix_to_top, [12] = set_name, [13] = output_folder_name, [14] = prefix_to_output_folder, [15] = command, 
	for folder_name in folder_list:
		info_dict = config_dict[folder_name]

		command = 'mkdir ' + folder_name
		os.popen(command)
		testscript.write(command + '\n')

		os.chdir(folder_name)
		testscript.write('cd ' + folder_name + '\n')

		command = 'cp -ar ../../../' + info_dict['inputDir'] + '* .'
		os.popen(command)
		testscript.write(command + '\n')

		command = 'ln -s ../../data/DATABASE ..'
		os.popen(command)
		testscript.write(command + '\n')

		command = 'sh ../../../utils/setup'
		os.popen(command)
		testscript.write(command + '\n')

		command = 'sh ./tidy '+info_dict['set_name']
		os.popen(command)
		testscript.write(command + '\n')

		os.chdir('..')
		testscript.write('cd ..\n')
	testscript.close()
