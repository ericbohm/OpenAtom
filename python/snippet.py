def snippet(num, inputfile, outputfile):
	in_file = open(inputfile)
	out_file = open(outputfile, 'w')
	if num > 1:
		while 2 > 1:
			line = in_file.readline()
			if line == '':
				print 'reach the end of the file'
				in_file.close()
				out_file.close()
				return
			index1 = line.find('Iteration',0,len(line)-1)
			index2 = line.find('done',0,len(line)-1)
			if index1 == -1 or index2 == -1:
 				continue
			index1 = index1 + 10
			index2 = index2 - 1
			while index1 < index2:
				iternum = ''
				iternum = iternum + line[index1]
				index1 = index1 + 1
			iternum = int(iternum)
			if iternum < num - 1:
				continue
			if iternum == num - 1:
				line = in_file.readline()
				line = in_file.readline()
				line = in_file.readline()
				while 2 > 1:
					line = in_file.readline()
					memflag = line.find('Memory',0,len(line)-1)
					equflag = line.find('=',0,len(line)-1)
					Psiflag = line.find('Psi[is=',0,len(line)-1)
					alldoneflag = line.find('allDoneCPForces',0,len(line)-1)
					if memflag == -1 and (equflag > -1 or Psiflag > -1) and alldoneflag == -1:
						out_file.write(line)
						continue
					if memflag > -1:
						in_file.close()
						out_file.close()
						return
			if iternum > num - 1:
				in_file.close()
				out_file.close()
				print 'Iteration not found!'
				return
	if num == 1:
		reachIter = 0
		while 2 > 1:
			line = in_file.readline()
			if line == '':
				in_file.close()
				out_file.close()
				return
			Psiflag = line.find('Psi[is=',0,len(line)-1)
			if reachIter == 0 and Psiflag > -1:
				out_file.write(line)
				while 2 > 1:
					line = in_file.readline()
					Psiflag = line.find('Psi[is=',0,len(line)-1)
					equflag = line.find('=',0,len(line)-1)
					memflag = line.find('Memory',0,len(line)-1)
					cpuflag = line.find('CpuTime',0,len(line)-1)
					alldoneflag = line.find('allDoneCPForces',0,len(line)-1)
					if memflag == -1 and (equflag > -1 or Psiflag > -1) and cpuflag == -1 and alldoneflag == -1:
						out_file.write(line)
						continue
					if memflag > -1:
						in_file.close()
						out_file.close()
						return

		
			
			

