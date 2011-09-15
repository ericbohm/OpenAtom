def compare_iteration(testoutput, reffile, magnum):
	import compare_number
	import stripNumbers
	testfile = open(testoutput)
	ref_file = open(reffile)
	test_content = testfile.readlines()
	ref_content = ref_file.readlines()
	testfile.close()
	ref_file.close()
	key_words = []
	option_key_words = []
	option_key_words.append('ENL')
	key_words.append('EHART')
	key_words.append('EExt')
	key_words.append('EWALD_recip')
	key_words.append('EEXC        ')
	key_words.append('EGGA        ')
	key_words.append('EEXC+EGGA')
	key_words.append('EKE')
	key_words.append('EWALD_REAL')
	key_words.append('atm fmag')
	key_words.append('MagForPsi')
	ref_dict = {}
	test_dict = {}
	ref_Psi_key_list = []
	test_Psi_key_list = []
	for ref_line in ref_content:
		if ref_line.find('Psi[is') != -1:
			ref_number_list = []
			ref_number_list = stripNumbers.stripNumbers(ref_line, 0)
			dict_key = str(ref_number_list[0])+' '+str(ref_number_list[1])+' '+str(ref_number_list[2])+' '+str(ref_number_list[3])
			new_list = []
			new_list.append(ref_number_list[4])
			new_list.append(ref_number_list[5])
			ref_Psi_key_list.append(dict_key)
			ref_dict[dict_key] = new_list
	for test_line in test_content:
		if test_line.find('Psi[is') != -1:
			test_number_list = []
			test_number_list = stripNumbers.stripNumbers(test_line, 0)
			dict_key = str(test_number_list[0])+' '+str(test_number_list[1])+' '+str(test_number_list[2])+' '+str(test_number_list[3])
			new_list = []
			new_list.append(test_number_list[4])
			new_list.append(test_number_list[5])
			test_Psi_key_list.append(dict_key)
			test_dict[dict_key] = new_list
	if len(ref_Psi_key_list) != len(test_Psi_key_list):
		print 'test output file is incomplete'
		return False
	for value in ref_Psi_key_list:
		if len(ref_dict[value]) != len(test_dict[value]):
			print 'test output file is incomplete'
			return False
		counter = 0
		while counter < len(ref_dict[value]):
			if compare_number.compare_number(ref_dict[value][counter], test_dict[value][counter], magnum) == False:
				return False
			counter = counter + 1
	for ref_line in ref_content:
		for keys in option_key_words:
			ref_number_list = []
			test_number_list = []
			found_op_key = 0
			if ref_line.find(keys) != -1:
				ref_number_list = stripNumbers.stripNumbers(ref_line, 0)
				for test_line in test_content:
					if test_line.find(keys) != -1:
						found_op_key = 1
						test_number_list = stripNumbers.stripNumbers(test_line, 0)
						break
			if found_op_key == 0:
				continue
			if len(ref_number_list) != len(test_number_list):
				return False
			counter = 0
			while counter < len(ref_number_list):
				if compare_number.compare_number(test_number_list[counter], ref_number_list[counter], magnum) == False:
					return False
				counter = counter + 1
		for keys in key_words:
			ref_number_list = []
			test_number_list = []
			if ref_line.find(keys) != -1:
				found_key = 0
				ref_number_list = stripNumbers.stripNumbers(ref_line, 0)
				for test_line in test_content:
					if test_line.find(keys) != -1:
						found_key = 1
						test_number_list = stripNumbers.stripNumbers(test_line, 0)
						break
				if found_key == 0:
					print 'output file is incomplete \n'
					return False
				if keys != 'MagForPsi':
					if len(ref_number_list) != len(test_number_list):
						return False
					counter = 0
					while counter < len(ref_number_list):
						if compare_number.compare_number(test_number_list[counter], ref_number_list[counter], magnum) == False:
							return False
						counter = counter + 1
				if keys == 'MagForPsi':
					if compare_number.compare_number(test_number_list[0], ref_number_list[0], magnum) == False:
						return False			
	return True
			
			
		
