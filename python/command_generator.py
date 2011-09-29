def command_generator(config_dict):
	folder_list = config_dict.keys()
	folder_list.remove('general_info')
# [0] = key, [1] = testname, [2] = machinename, [3] = inputDir, [4] = inputPar, [5] = inputPhy, [6] = outRef, [7] = test_output, [8] = numpe, [9] = prefix_to execute, [10] = prefix_to_data, [11] = prefix_to_top, [12] = set_name, [13] = output_folder_name, [14] = prefix_to_output_folder, [15] = command, 
	for folder_name in folder_list:
		info_dict = config_dict[folder_name]
		numpe = info_dict['numpe']
		data_path = info_dict['data_path']
		Par_value = info_dict['inputPar']
		Phy_value = info_dict['inputPhy']
		outputfile = info_dict['test_output']
		command = '../../charmrun +p'+numpe+' ++local ../../OpenAtom regression/'+Par_value+' regression/'+Phy_value+' > '+outputfile
		config_dict[folder_name]['command'] = command
	return config_dict



