def config_reader(filename, executable_path):
	import yaml
	stream = file(filename)
	configfile = yaml.load(stream)
	stream.close()
	config_dict={}
	config_dict['general_info'] = configfile[0]
	testnum = len(configfile)
	testcounter = 1
	prefix_to_execute = '../../'
	data_path = 'regression/'
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
		if configfile[testcounter].has_key('keys') == False:
			keylist.append('default_key')
		for key in keylist:
			for numpe in numpelist:
				info_dict = {}
				test_name = configfile[testcounter]['name'];
				folder_name = test_name + '_' + key + '_' + str(numpe)
				info_dict['key'] = key
				info_dict['test_name'] = test_name
				info_dict['machine_name'] = configfile[testcounter]['machine names'];
				info_dict['inputDir'] = configfile[testcounter]['inputDir']
				info_dict['inputPar'] = configfile[testcounter]['inputPar'].replace('$K', key).replace('$P', str(numpe))
				info_dict['inputPhy'] = configfile[testcounter]['inputPhy'].replace('$K', key).replace('$P', str(numpe))
				info_dict['outRef'] = configfile[testcounter]['outRef'].replace('$K', key).replace('$P', str(numpe))
				testoutput = folder_name + '.result'
				info_dict['test_output'] = testoutput
				info_dict['numpe'] = str(numpe)
				info_dict['executable_path'] = executable_path
				info_dict['set_name'] = configfile[testcounter]['set_name']
				info_dict['output_folder'] = output_folder_name
				info_dict['data_path'] = data_path		
				config_dict[folder_name] = info_dict
		testcounter = testcounter + 1
	return config_dict					

