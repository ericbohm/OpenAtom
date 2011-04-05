def checkresult(testoutput, reffile, iternum, magnum):
	import snippet
	import compare_2
	import os
	testtemp = testoutput+'temp'
	snippet.snippet(iternum, testoutput, testtemp)
	print compare_2.compare_2(testtemp, reffile, magnum)
	os.remove(testtemp)

