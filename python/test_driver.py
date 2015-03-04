import os
import subprocess
import sys

import test_parser
import test_runner

tests = []

# Read the tests from the passed in config file
if len(sys.argv) == 3:
	print 'Reading test configuration from ' + sys.argv[2] + '\n'
	f = open(sys.argv[2], 'r')
	for line in f:
		args = line.split()
		tests = tests + test_parser.parse_tests(sys.argv[1], args)
	f.close()
else:
	print 'ERROR: Must specify a build path and a config file!'
	print 'EXAMPLE: python test_driver.py ../build-O3 make_test.config'
	sys.exit(1)

num_failed = test_runner.execute_tests(tests)
sys.exit(num_failed)
