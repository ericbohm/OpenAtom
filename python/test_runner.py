def print_summary(results):
	import string

	header = string.ljust('TEST NAME',20)
	header = header + string.ljust('VARIANT',20)
	header = header + string.ljust('PES',4) + 'RESULT'

	print string.ljust('', 80, '=')
	print header
	print string.ljust('', 80, '-')

	passed = 0
	failed = 0

	for i in range(len(results)):
		if i != 0 and results[i][0] == results[i-1][0]:
			summ_str = make_summary_string(results[i],False)
		else:
			summ_str = make_summary_string(results[i],True)
			if i != 0:
				print ''
		print summ_str

		if results[i][3]:
			passed = passed + 1
		else:
			failed = failed + 1

	print string.ljust('', 80, '=')
	print `passed` + ' PASSED\t' + `failed` + ' FAILED'

def make_summary_string(result, print_name):
	import string
	summ_str = ''
	if print_name:
		summ_str = string.ljust(result[0],20)
	else:
		summ_str = string.ljust('',20)
	summ_str = summ_str + string.ljust(result[1],20) + string.ljust(result[2],4)
	if result[3]:
		summ_str = summ_str + 'PASSED'
	else:
		summ_str = summ_str + 'FAILED: ' + result[4]
	return summ_str

import config_parser
import subprocess
import sys
import os
runs = []
f = open(sys.argv[1], 'r')
for line in f:
	runs.append(line.split())
f.close()

base_directory = os.getcwd()

results = []

for run in runs:
	os.chdir(base_directory)
	results = results + config_parser.config_parser(run[0], run[1], run[2])

print_summary(results)
