def config_reader(filename, set_name):
	import yaml
	import os
	import checkresult
# this section reads the config file and store it into a dictionary
	stream = file(filename)
	configfile = yaml.load(stream)
	stream.close()
#############################################################
# since the first key is the general config information, the real test number is the length of the dictionary - 1
	testnum = len(configfile) - 1
#############################################################
# the first item in the dictionary is the general config information so the testcounter starts from 1
	testcounter = 1
# summary_True list stores the successful tests' name
	summary_True = []
# summary_False list stores the failed tests' name
	summary_False = []
#start to compare
	while testcounter <= testnum:
		testname = configfile[testcounter]['name']
		machine_name_list = configfile[testcounter]['machine names']
		keylist = []
#key is an optional attribute, we have to decide wheather the attribute is existed or not, the same to numpe
		if configfile[testcounter].has_key('keys'):
			keylist = configfile[testcounter]['keys']
		if configfile[testcounter].has_key('numpe'):
			numpelist = []
			for value in configfile[testcounter]['numpe']:
				numpelist.append(str(value))
# inputParlist stores the inputPar values
		inputParlist = []
# inputPhylist stores the inputPhy values
		inputPhylist = []
# outReflist stores Reference file names
		outReflist = []
# commandlist stores commands that launch the program with only key words inside (like %K, %T...)
		commandlist = []
# finalcommandlist stores commands that launch the program, with key words removed
		finalcommandlist = []
# outputlist stores the tempoutput files' names with only key workds inside (like %K, %T...)
		outputlist = []
# finaloutputlist stores tempoutput files' names, with key words removed
		finaloutputlist = []
# suppose there is a key
		mykey = '$K'
# otherwise mykey =  NULL
		if len(keylist) == 0:
			mykey = ''
# refdict stores {output file, corresponding reference file} pair 
		refdict = {}
		inputParlist.append(configfile[testcounter]['inputPar'])
		inputPhylist.append(configfile[testcounter]['inputPhy'])
		inputParlist = list(set(inputParlist)) #delete duplicates
		inputPhylist = list(set(inputPhylist))
		outReflist = list(set(outReflist))
# remove the key words, (like %K, %T...), suppose there is a key existed
		if len(keylist) > 0:
			for P_value in numpelist:
				for Par_value in inputParlist:
					for Phy_value in inputPhylist:
						if (Phy_value.find('$T', 0, len(Phy_value)) >= 0):
							mykey = '$T'
						command = '../charmrun +p$P ++local ../OpenAtom regression/'+Par_value+' '+'regression/'+Phy_value+' > '+'temp.'+configfile[testcounter]['name']+'.'+mykey+'.$P.result'
# THE output files name
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
# store {output file, corresponding reference file} pair into the refdict dictionary
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

# remove the key words, (like %K, %T...), suppose there is no key existed
		if len(keylist) == 0:
			for P_value in numpelist:
				for Par_value in inputParlist:
					for Phy_value in inputPhylist:
						if (Phy_value.find('$T', 0, len(Phy_value)) >= 0):
							mykey = '$T'
						command = command = '../charmrun +p$P ++local ../OpenAtom regression/'+Par_value+' '+'regression/'+Phy_value+' > '+'temp.'+configfile[testcounter]['name']+'.'+mykey+'.$P.result'
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
#delete duplicates
		commandlist = list(set(finalcommandlist))
# sort the results
		outReflist.sort()
		outputlist = list(set(finaloutputlist))
		outputlist.sort()
		for value in commandlist:
			print value
# prepare to execute the command, change the working directory to the openatom working directory and prepare output directory
			os.chdir('..')
			os.chdir('build-O3')
			os.popen('mkdir test-output')
			os.chdir('test-output')
			os.popen('ln -s ../../' + configfile[testcounter]['inputDir'] + '* .')
			os.popen('ln -s ../data/DATABASE ..')
			os.popen('sh ../../utils/setup')
			os.popen('sh ./tidy '+set_name)
			os.popen(value)
# change the working directory back to python directory, prepare next iteration
			os.chdir('..')
			os.chdir('..')
			os.chdir('python')
# finished all the openatom commands, starts to compare the results
		lengh = len(outputlist)
		counter = 0
		while counter < lengh:
# output files' path
			testoutput = '../build-O3/test-output/' + outputlist[counter]
# reference files' path
			reffile = '../build-O3/test-output/regression/' + refdict[outputlist[counter]]
			print 'comparing file :' + outputlist[counter]
# check the results
			test_result = checkresult.checkresult(testoutput, reffile, configfile[0]['iternum'], configfile[0]['numSigDigits'])
# prepare the summary
			if test_result == True:
				summary_True.append(outputlist[counter])
			if test_result == False:
				summary_False.append(outputlist[counter])
			counter = counter + 1
		testcounter=testcounter + 1
# clear the {output file, corresponding reference file} pair dictionary, prepare for the next test iteration
		refdict.clear()
		os.popen('rm -rf ../build-O3/test-output/regression')
		os.popen('rm ../build-O3/DATABASE')
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
import sys
configfilename = sys.argv[1]
set_name = sys.argv[2]
config_reader(configfilename, set_name)
				


					
