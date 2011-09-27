def command_generator(config_dict):
	folder_list = config_dict.keys()
	folder_list.remove('general_info')
# [0] = key, [1] = testname, [2] = machinename, [3] = inputDir, [4] = inputPar, [5] = inputPhy, [6] = outRef, [7] = testoutput, [8] = numpe, [9] = prefix_to execute, [10] = prefix_to_data, [11] = prefix_to_top, [12] = set_name, [13] = output_folder_name, [14] = prefix_to_output_folder, [15] = command, 
	for folder_name in folder_list:
		prefix_to_execute = config_dict[folder_name][9]
		numpe = config_dict[folder_name][8]
		prefix_to_data = config_dict[folder_name][10]
		Par_value = config_dict[folder_name][4]
		Phy_value = config_dict[folder_name][5]
		outputfile = config_dict[folder_name][7]
		command = prefix_to_execute + 'charmrun +p'+numpe+' ++local '+prefix_to_execute+'OpenAtom '+prefix_to_data+Par_value+' '+prefix_to_data+Phy_value+' > '+outputfile
		config_dict[folder_name].append(command)
	return config_dict



