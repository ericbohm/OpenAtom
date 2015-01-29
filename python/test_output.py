def print_table_header():
	import string

	header = string.ljust('NUMBER', 10) + string.ljust('TEST NAME',20)
	header = header + string.ljust('VARIANT',20)
	header = header + string.ljust('PES',4) + 'RESULT'

	print string.ljust('', 80, '=')
	print header
	print string.ljust('', 80, '-')

def print_summary_string(test_num, total_num, test_dict):
	import string
	import sys
	summ_string = '[' + string.rjust(str(test_num),len(str(total_num)))
	summ_string = summ_string + '/' + str(total_num) + ']'
	summ_string = string.ljust(summ_string, 10)

	summ_string = summ_string + string.ljust(test_dict['test_name'],20)
	summ_string = summ_string + string.ljust(test_dict['variant'],20)
	summ_string = summ_string + string.ljust(test_dict['numpe'],4)

	sys.stdout.write(summ_string)
	sys.stdout.flush()

def print_result_string(test_dict):
	if test_dict['result'][0]:
		print 'PASSED'
	else:
		print 'FAILED - ' + test_dict['result'][1]

def print_table_footer(passed, failed):
	import string

	print string.ljust('', 80, '=')
	print 'PASSED ' + `passed` + '\tFAILED ' + `failed` + '\n'	
