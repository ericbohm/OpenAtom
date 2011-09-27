def snippet(num, inputfile, outputfile):
	in_file = open(inputfile)
	out_file = open(outputfile, 'w')
# prepare the key word to identify iteration num, for example, Iter [0]
	iter_num = 'Iter ['+str(num)+']'
# bad_key list contains all the key that will appear in the output file but won'a appear in reference file.
	bad_key = []
	bad_key.append('allDoneCPForces')
	bad_key.append('EWALD_SELF')
	bad_key.append('EWALD_BGR')
	bad_key.append('Iteration time')
# list lines contains all the contents in the inputfile
	lines = in_file.readlines();
	in_file.close()
# to check every line
	for line in lines:
# if there is a iteration identifier in the line and there is no bad_key in the line then append the line in to the outputfile.
		if line.find(iter_num,0,len(line)-1) > -1:
			bad_key_flag= 0
			for value in bad_key:
				bad_key_flag = bad_key_flag + line.find(value,0,len(line)-1)
			if bad_key_flag != 0-(len(bad_key)):
				continue
# this part is needed because we don't want 'Iter [0] ' words appear in the output file.
			counter = 9
			line_new = ''
			while counter < len(line) - 1:
				line_new = line_new + line[counter]
				counter = counter + 1
			line_new = line_new + '\n'
			out_file.write(line_new)
	out_file.close()

		
			
	
