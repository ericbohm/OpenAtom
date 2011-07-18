def config_reader(filename):
	import yaml
	import os
	import checkresult
	stream = file(filename)
	configfile = yaml.load(stream)
	stream.close()
	testnum = len(configfile) - 1
	testcounter = 1
	while testcounter <= testnum:
		machine_name_list = configfile[testcounter]['machine names']
		keylist = []
		if configfile[testcounter].has_key('keys'):
			keylist = configfile[testcounter]['keys']
		if configfile[testcounter].has_key('numpe'):
			numpelist = []
			for value in configfile[testcounter]['numpe']:
				numpelist.append(str(value))
		inputParlist = []
		inputPhylist = []
		outReflist = []
		commandlist = []
		finalcommandlist = []
		outputlist = []
		finaloutputlist = []
		mykey = '$K'
		if len(keylist) == 0:
			mykey = ''
		refdict = {}
		inputParlist.append(configfile[testcounter]['inputPar'])
		inputPhylist.append(configfile[testcounter]['inputPhy'])
		inputParlist = list(set(inputParlist)) #delete duplicates
		inputPhylist = list(set(inputPhylist))
		outReflist = list(set(outReflist))
		if len(keylist) > 0:
			for P_value in numpelist:
				for Par_value in inputParlist:
					for Phy_value in inputPhylist:
						if (Phy_value.find('$T', 0, len(Phy_value)) >= 0):
							mykey = '$T'
						command = '../../charmrun +p$P ++local ../../OpenAtom regression/'+Par_value+' '+'regression/'+Phy_value+' > '+'temp.'+configfile[testcounter]['name']+'.'+mykey+'.$P.result'
						outputval = 'temp.'+configfile[testcounter]['name']+'.'+mykey+'.$P.result'
						outputlist.append(outputval)
						commandlist.append(command)
			if len(machine_name_list) == 0:
				for kvalue in keylist:
					for cvalue in commandlist:
						for pvalue in numpelist:
							finalcommandlist.append(cvalue.replace('$K', kvalue).replace('$T', configfile[testcounter]['name']).replace('$P', pvalue).replace('$M', ''))
				for kvalue in keylist:
					for cvalue in outputlist:
						for pvalue in numpelist:
							tempoutval = cvalue.replace('$K', kvalue).replace('$T', configfile[testcounter]['name']).replace('$P', pvalue)
							finaloutputlist.append(tempoutval)
							temprefval = configfile[testcounter]['outRef'].replace('$K', kvalue).replace('$T', configfile[testcounter]['name']).replace('$P', pvalue)
							refdict[tempoutval] = temprefval		
			if len(machine_name_list) > 0:
				for kvalue in keylist:
					for cvalue in commandlist:
						for pvalue in numpelist:
							for mvalue in machine_name_list:
								finalcommandlist.append(cvalue.replace('$K', kvalue).replace('$T', configfile[testcounter]['name']).replace('$P', pvalue).replace('$M', '/'+mvalue))
				for kvalue in keylist:
					for cvalue in outputlist:
						for pvalue in numpelist:
							for mvalue in machine_name_list:
								finaloutputlist.append(cvalue.replace('$K', kvalue).replace('$T', configfile[testcounter]['name']).replace('$P', pvalue).replace('$M', '/'+mvalue))


		if len(keylist) == 0:
			for P_value in numpelist:
				for Par_value in inputParlist:
					for Phy_value in inputPhylist:
						if (Phy_value.find('$T', 0, len(Phy_value)) >= 0):
							mykey = '$T'
						command = '../../charmrun +p$P ++local ../../OpenAtom regression/'+Par_value+' '+'regression/'+Phy_value+' > '+'temp.'+configfile[testcounter]['name']+'.'+mykey+'.$P.result'
						outputval = 'temp.'+configfile[testcounter]['name']+'.'+mykey+'.$P.result'
						outputlist.append(outputval)
						commandlist.append(command)
			if len(machine_name_list) == 0:
				for cvalue in commandlist:
					for pvalue in numpelist:
						finalcommandlist.append(cvalue.replace('$T', configfile[testcounter]['name']).replace('$P', pvalue).replace('$M', ''))
				for cvalue in outputlist:
					for pvalue in numpelist:
						tempoutval = cvalue.replace('$T', configfile[testcounter]['name']).replace('$P', pvalue)
						finaloutputlist.append(tempoutval)
						temprefval = configfile[testcounter]['outRef'].replace('$T', configfile[testcounter]['name']).replace('$P', pvalue)
						refdict[tempoutval] = temprefval		
			if len(machine_name_list) > 0:
				for cvalue in commandlist:
					for pvalue in numpelist:
						for mvalue in machine_name_list:
							finalcommandlist.append(cvalue.replace('$T', configfile[testcounter]['name']).replace('$P', pvalue).replace('$M', '/'+mvalue))
				for cvalue in outputlist:
					for pvalue in numpelist:
						for mvalue in machine_name_list:
							finaloutputlist.append(cvalue.replace('$T', configfile[testcounter]['name']).replace('$P', pvalue).replace('$M', '/'+mvalue))
		commandlist = list(set(finalcommandlist))
		outReflist.sort()
		outputlist = list(set(finaloutputlist))
		outputlist.sort()
		for value in commandlist:
			os.chdir('..')
			os.chdir('build-O3')
			os.system('mkdir test-output')
			os.system('mkdir test-output/regression')
			os.system('mkdir test-output/regression/STATES_OUT')
			os.chdir('test-output/regression')
			os.system('ln -s ../../../' + configfile[testcounter]['inputDir'] + '* .')
			os.system('ln -s ../../data/DATABASE ..')
			for sd in os.listdir('STATES'):
				sdpath=os.path.join('STATES_OUT',sd)
				os.system('rm -r '+ sdpath)
				os.system('mkdir '+ sdpath)
			for ad in os.listdir('ATOM_COORDS_IN'):
				adpath=os.path.join('ATOM_COORDS_OUT',ad)
				os.system('rm -r '+ adpath)
				os.system('mkdir '+ adpath)
			os.system(value)
			#print value
			os.system('./tidy water')
			os.chdir('..')
			os.chdir('..')
			os.chdir('..')
			os.chdir('python')
		lengh = len(outputlist)
		counter = 0
		while counter < lengh:
			testoutput = '../build-O3/test-output/regression/' + outputlist[counter]
			reffile = '../build-O3/test-output/regression/regression/' + refdict[outputlist[counter]]
			print 'comparing file :' + outputlist[counter]
			checkresult.checkresult(testoutput, reffile, configfile[0]['iternum'], configfile[0]['numSigDigits'])
			counter = counter + 1
		testcounter=testcounter + 1
		refdict.clear()
import sys
configfilename = sys.argv[1]
config_reader(configfilename)
				


					
