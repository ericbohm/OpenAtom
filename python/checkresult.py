def checkresult(testoutput, reffile, iternum, magnum):
	import snippet
	import compare_iteration
	import os
#prepare the temp file's name
	testtemp = testoutput+'temp'
#snips the specified iteration into the temp test file
	snippet.snippet(iternum, testoutput, testtemp)
#compare the the iteration output with the reference file
	returnval = compare_iteration.compare_iteration(testtemp, reffile, magnum)
	print returnval
	os.remove(testtemp)
	return returnval

