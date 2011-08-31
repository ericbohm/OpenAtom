def compare_iteration(testoutput, reffile, magnum):
	import compare_number
	import stripNumbers
	testfile = open(testoutput)
	ref_file = open(reffile)
	counter = 0
	testcontent = []
	refcontent = []
	#read file contents into a list
	while 2 > 1:
		line = testfile.readline()
		if line == '':
			break
		testcontent.append(line)
	testfile.close()
	if len(testcontent) == 0:
		print 'Error, empty test output file.'
		return False
	counter = 0
	while 2 > 1:
		line = ref_file.readline()
		if line == '':
			break
		refcontent.append(line)
	ref_file.close()
	if len(refcontent) == 0:
		print 'Error, empty ref file.'
		return False
	#start compare
	while 2 > 1:
		testlen = len(testcontent)
		reflen = len(refcontent)
		if testlen == 0 or reflen == 0:
			return True
		if testcontent[0].find('Psi[',0,len(testcontent[0])-1) < 0:
			#find the keyword
			keyword = ''
			counter = 4
			while 2 > 1:
				if testcontent[0][counter] == ' ':
					break
				keyword = keyword + testcontent[0][counter]
				counter = counter + 1
			#find the same keyword in the ref file and compare the number with ref file
			counter = 0
			while 2 > 1:
				if counter == reflen:
					del testcontent[0]
					break
				# if found the result
				if refcontent[counter].find(keyword,0,len(refcontent[counter])-1) >= 0:
					testlist = stripNumbers.stripNumbers(testcontent[0], 0)
					reflist = stripNumbers.stripNumbers(refcontent[counter], 0)
					#return false if the output is different from the ref file
					if compare_number.compare_number(testlist[1], reflist[0], magnum) == False:
						return False
					#if the result if the same then delete the compared part
					del refcontent[counter]
					del testcontent[0]
					break
				if refcontent[counter].find(keyword,0,len(refcontent[counter])-1) < 0:
					counter = counter + 1
			continue
		if testcontent[0].find('Psi[',0,len(testcontent[0])-1) >= 0:
			counter = 4
			keyword = ''
			# find the key word
			while 2 > 1:
				if testcontent[0][counter] == ']':
					break
				if testcontent[0][counter] != ' ':
					keyword = keyword + testcontent[0][counter]
				counter = counter + 1
			#locate the result in the ref file
			counter = 0
			while 2 > 1:
				if counter == reflen:
					del testcontent[0]
					break
				if refcontent[counter].find('Psi[',0,len(refcontent[counter])-1) >= 0:
					linelen = len(refcontent[counter])
					newline = ''
					counter2 = 0
					while 2 > 1:
						if counter2 == linelen:
							break
						if refcontent[counter][counter2] != ' ':
							newline = newline + refcontent[counter][counter2]
						counter2 = counter2 + 1
					if newline.find(keyword,0,len(newline))>= 0:
						break
					counter = counter + 1
					continue
				if refcontent[counter].find('Psi[',0,len(refcontent[counter])-1) < 0:
					counter = counter + 1
			testlist = stripNumbers.stripNumbers(testcontent[0], 0)
			del testlist[0]
			reflist = stripNumbers.stripNumbers(refcontent[counter], 0)
			tempcounter = counter
			counter = 4
			testlistlen = len(testlist)
			reflistlen = len(reflist)
			while 2 > 1:
				if counter == reflistlen or counter == testlistlen:
					del refcontent[tempcounter]
					del testcontent[0]
					break
				if compare_number.compare_number(testlist[counter], reflist[counter], magnum) == False:
					return False
				counter = counter + 1


				
					

				
						
						

								
			

