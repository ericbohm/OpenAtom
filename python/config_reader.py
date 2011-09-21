def config_reader(filename, set_name):
	import yaml
	stream = file(filename)
	configfile = yaml.load(stream)
	stream.close()
	config_dict={}
	config_dict['general_info'] = configfile[0]
	testnum = len(configfile)
	testcounter = 1
	prefix_to_execute = '../../'
	prefix_to_data = 'regression/'
	prefix_to_top = '../'
	output_folder_name = 'test-output'
	prefix_to_output_folder = '../build-O3/'+output_folder_name+'/'
# [0] = key, [1] = testname, [2] = machinename, [3] = inputDir, [4] = inputPar, [5] = inputPhy, [6] = outRef, [7] = testoutput, [8] = numpe, [9] = prefix_to execute, [10] = prefix_to_data, [11] = prefix_to_top, [12] = set_name, [13] = output_folder_name, [14] = prefix_to_output_folder, [15] = command, 
	while testcounter < testnum:
		keylist = []
		numpelist = []
		numpelist = configfile[testcounter]['numpe']
		if configfile[testcounter].has_key('keys'):
			keylist = configfile[testcounter]['keys']
		for key in keylist:
			for numpe in numpelist:
				info_list = []
				test_name = configfile[testcounter]['name'];
				folder_name = test_name + '_' + key + '_' + str(numpe)
				info_list.append(key);
				info_list.append(test_name);
				info_list.append(configfile[testcounter]['machine names'])
				info_list.append(configfile[testcounter]['inputDir'])
				info_list.append(configfile[testcounter]['inputPar'].replace('$K', key).replace('$P', str(numpe)))
				info_list.append(configfile[testcounter]['inputPhy'].replace('$K', key).replace('$P', str(numpe)))
				info_list.append(configfile[testcounter]['outRef'].replace('$K', key).replace('$P', str(numpe)))
				testoutput = folder_name + '.result'
				info_list.append(testoutput)
				info_list.append(str(numpe))
				info_list.append(prefix_to_execute)
				info_list.append(prefix_to_data)
				info_list.append(prefix_to_top)
				info_list.append(set_name)
				info_list.append(output_folder_name)
				info_list.append(prefix_to_output_folder)				
				config_dict[folder_name] = info_list;
				info_list = []
		testcounter = testcounter + 1
	return config_dict					

