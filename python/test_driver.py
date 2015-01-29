import os
import subprocess
import sys

import test_parser
import test_runner

tests = []

# Read the tests from the passed in config file
if len(sys.argv) == 2:
	print 'Reading test configuration from ' + sys.argv[1] + '\n'
	f = open(sys.argv[1], 'r')
	for line in f:
		args = line.split()
		tests = tests + test_parser.parse_tests(args[0], args[1], args[2])
	f.close()
elif len(sys.argv) == 4:
	print 'Reading test configuration for (' + sys.argv[1] + ',' + sys.argv[2] + ',' + sys.argv[3] + ')\n'
	tests = test_parser.parse_tests(sys.argv[1],sys.argv[2],sys.argv[3])
else:
	print 'ERROR: Must specify either a config file, or a single test'
	sys.exit(-1)

num_failed = test_runner.execute_tests(tests)
sys.exit(num_failed)
