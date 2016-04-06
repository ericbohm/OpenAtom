def check_result(test_dict):
	import os

	output_file = test_dict['output_file']
	out_ref = test_dict['out_ref']
	iteration = test_dict['iteration']
	sig_figs = test_dict['sig_figs']

	# Snip the iteration being checked out into a temporary file
	output_temp = output_file+'.temp'
	snip_iteration(iteration, output_file, output_temp)

	# Compare the output reference to the temporary output
	result = compare_iteration(output_temp, out_ref, sig_figs)

	return result

def snip_iteration(num, inputfile, outputfile):
	in_file = open(inputfile)
	out_file = open(outputfile, 'w')

	# Iterations marked in the output file start with 'Iter [x]'
	iter_num = 'Iter ['+str(num)+']'

	# Bad keys are things in the output file that we don't want to check, ie time
	bad_key = []
	bad_key.append('allDoneCPForces')
	bad_key.append('EWALD_SELF')
	bad_key.append('EWALD_BGR')
	bad_key.append('Iteration time')

	# Pull out all good lines for the desired iteration and append them to file
	lines = in_file.readlines();
	in_file.close()
	for line in lines:
		if line.find(iter_num,0,len(line)-1) > -1:
			bad_key_flag= 0
			for value in bad_key:
				bad_key_flag = bad_key_flag + line.find(value,0,len(line)-1)
			if bad_key_flag != 0-(len(bad_key)):
				continue

			# Strip off the 'Bead and Iter prefixes ' before each line
			counter = line.find(iter_num,0,len(line)-1) + len(iter_num) + 1
			line_new = ''
			while counter < len(line) - 1:
				line_new = line_new + line[counter]
				counter = counter + 1
			line_new = line_new + '\n'
			out_file.write(line_new)
	out_file.close()

def compare_iteration(testoutput, reffile, magnum):
	testfile = open(testoutput)
	ref_file = open(reffile)
	test_content = testfile.readlines()
	ref_content = ref_file.readlines()
	testfile.close()
	ref_file.close()
	key_words = []

	# We are looking to match all of the following keywords in both files
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

	# First compare all Psi values found in both files
	for ref_line in ref_content:
		# First get all Psis from the reference file
		if ref_line.find('Psi[is') != -1:
			ref_number_list = []
			ref_number_list = strip_numbers(ref_line, 0)
			dict_key = str(ref_number_list[0])+' '+str(ref_number_list[1])+' '+str(ref_number_list[2])+' '+str(ref_number_list[3])
			new_list = []
			new_list.append(ref_number_list[4])
			new_list.append(ref_number_list[5])
			ref_Psi_key_list.append(dict_key)
			ref_dict[dict_key] = new_list
	# Then get all Psis from the test output
	for test_line in test_content:
		if test_line.find('Psi[is') != -1:
			test_number_list = []
			test_number_list = strip_numbers(test_line, 0)
			dict_key = str(test_number_list[0])+' '+str(test_number_list[1])+' '+str(test_number_list[2])+' '+str(test_number_list[3])
			new_list = []
			new_list.append(test_number_list[4])
			new_list.append(test_number_list[5])
			test_Psi_key_list.append(dict_key)
			test_dict[dict_key] = new_list

	# If Psi values are don't match 1-1, then return a failure
	if len(ref_Psi_key_list) != len(test_Psi_key_list):
		return [False, 'Unequal number of Psi values in output and ref file']
	for value in ref_Psi_key_list:
		if len(ref_dict[value]) != len(test_dict[value]):
			return [False, 'Psi values don\'t match']
		counter = 0
		while counter < len(ref_dict[value]):
			if compare_number(ref_dict[value][counter], test_dict[value][counter], magnum) == False:
				return [False, 'Psi values do not match']
			counter = counter + 1

	# Now we check the values of optional keywords
	for ref_line in ref_content:
		for keys in option_key_words:
			ref_number_list = []
			test_number_list = []
			found_op_key = 0
			# If we found a key in the ref file, find it in the test output
			if ref_line.find(keys) != -1:
				ref_number_list = strip_numbers(ref_line, 0)
				for test_line in test_content:
					if test_line.find(keys) != -1:
						found_op_key = 1
						test_number_list = strip_numbers(test_line, 0)
						break

			# If you didn't find the key at all, move on. Otherwise check results.
			if found_op_key == 0:
				continue

			if len(ref_number_list) != len(test_number_list):
				return [False, 'Key ' + keys + ' doesn\'t match']

			counter = 0
			while counter < len(ref_number_list):
				if compare_number(test_number_list[counter], ref_number_list[counter], magnum) == False:
					return [False, 'Key ' + keys + ' doesn\'t match']
				counter = counter + 1

		for keys in key_words:
			ref_number_list = []
			test_number_list = []
			# If we found the key in the ref file, find it in the test output
			if ref_line.find(keys) != -1:
				found_key = 0
				ref_number_list = strip_numbers(ref_line, 0)
				for test_line in test_content:
					if test_line.find(keys) != -1:
						found_key = 1
						test_number_list = strip_numbers(test_line, 0)
						break

				# If we couldn't find the key, or they don't match, the test failed
				if found_key == 0:
					return [False, 'Missing key in output: ' + keys]
				if keys != 'MagForPsi':
					if len(ref_number_list) != len(test_number_list):
						return [False, 'Key ' + keys + ' doesn\'t match']
					counter = 0
					while counter < len(ref_number_list):
						if compare_number(test_number_list[counter], ref_number_list[counter], magnum) == False:
							return [False, 'Key ' + keys + ' doesn\'t match']
						counter = counter + 1
				if keys == 'MagForPsi':
					if compare_number(test_number_list[0], ref_number_list[0], magnum) == False:
						return [False, 'Key ' + keys + ' doesn\'t match']

	return [True, '']

def compare_number(num1, num2, numsigdig):
	tempnum1 = abs(num1)
	tempnum2 = abs(num2)
	c = 0
	index = 0
	reachdigit = 0
	returnval = False
	if num1 < 0 and num2 > 0:
		return False
	if num1 > 0 and num2 < 0:
		return False
	if num1 == num2:
		return True
	if tempnum1 < 1 and tempnum2 > 1:
		return False
	if tempnum1 > 1 and tempnum2 < 1:
		return False
	tempnum1 = str(tempnum1)
	tempnum2 = str(tempnum2)
	len1 = len(tempnum1)
	len2 = len(tempnum2)
	while c < numsigdig:
		if index < len1 and index < len2:
			if tempnum1[index] == '.' and tempnum2[index] == '.':
				index = index + 1
				continue
			if tempnum1[index] == tempnum2[index] and reachdigit == 0 and (tempnum1[index] == '0' or tempnum1[index] == '.'):
				index = index + 1
				continue
			if tempnum1[index] == tempnum2[index] and reachdigit == 0 and tempnum1[index] != '0' and tempnum1[index] != '.':
				reachdigit = 1
				index = index + 1
				c = c + 1
				continue
			if tempnum1[index] == tempnum2[index] and reachdigit == 1:
				index = index + 1
				c = c + 1
				continue
			if tempnum1[index] != tempnum2[index]:
				return False
		else:
			print "Warning : the precision of input data is not enough"
			return False
	return True

def strip_numbers(inputstr, method):
	length = len(inputstr)
	strcounter = 0
	tempstr = ""
	returnlist = []
	while strcounter < length:
		if (inputstr[strcounter] < '0' or inputstr[strcounter] > '9') and inputstr[strcounter] != '.':
			tempstr = ""
		else:
			tempstr = tempstr + inputstr[strcounter]
			if strcounter + 1 < length:
				if (inputstr[strcounter + 1] < '0' or inputstr[strcounter + 1] > '9') and inputstr[strcounter + 1] != '.':
						if method == 0:
							returnlist.append(float(tempstr))
						if method == 1:
							returnlist.append(int(tempstr))
						if method == 2:
							returnlist.append(tempstr)
			if strcounter + 1 == length:
					if method == 0:
						returnlist.append(float(tempstr))
					if method == 1:
						returnlist.append(int(tempstr))
					if method == 2:
						returnlist.append(tempstr)
		strcounter = strcounter + 1;	
	return returnlist
