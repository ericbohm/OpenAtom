def summary_generator(filename):
	import yaml
	import checkresult
	dict_file = file(filename)
	config_dict = yaml.load(dict_file)
	dict_file.close()
	keylist = config_dict.keys()
	keylist.remove('general_info')
	summary_True = []
	summary_False = []
	for folder_name in keylist:
		info_dict = config_dict[folder_name]
		outputfile = info_dict['executable_path']+'/'+info_dict['output_folder']+'/'+folder_name+'/'+info_dict['test_output']
		reffile = info_dict['executable_path']+'/'+info_dict['output_folder']+'/'+folder_name+'/regression/'+info_dict['outRef']
		print 'comparing test :' + folder_name
		iternum = config_dict['general_info']['iternum']
		numSigDigits = config_dict['general_info']['numSigDigits']
		test_result = checkresult.checkresult(outputfile, reffile, iternum, numSigDigits)
		if test_result == True:
			summary_True.append(folder_name)
		if test_result == False:
			summary_False.append(folder_name)
	print '\n'
	print '################################################################'
	print '				SUMMARY:'
	print 'successful tests:\n'
	for value in summary_True:
		print value
	print '\n'
	print 'failed tests:\n'
	for value in summary_False:
		print value
	print '################################################################'
		
