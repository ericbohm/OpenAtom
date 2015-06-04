def parse_tests(exe_dir, test_args):
	import os
	import yaml

	# config_dict will contain entries common to every test in this config file
	config_dict = {}

	exe_path = os.path.abspath(exe_dir)
	data_path = ''
	test_path = ''
	test_name = ''

	if (len(test_args) > 0):
		data_path = os.path.abspath(test_args[0])
	if (len(test_args) > 1):
		test_path = os.path.join(data_path,'tests',test_args[1])
	if (len(test_args) > 2):
		test_name = test_args[2]
	if (len(test_args) <= 0 or len(test_args) > 3):
		print "ERROR: Bad number of test arguments: " + str(test_args)
		return 1

	if (len(test_args) != 2):
		print "ERROR: Currently we only support running tests with two arguments"
		print "EXAMPLE: ../data/water_32M_10Ry regression"
		return 1

	config_dict['data_path'] = data_path
	config_dict['exe_path'] = exe_path
	config_dict['charmrun'] = os.path.join(exe_path,'charmrun')
	config_dict['exe'] = os.path.join(exe_path,'OpenAtom')
	config_dict['test_path'] = test_path
	config_dict['output_path'] = os.path.join(test_path,'output')
	config_dict['config_file'] = os.path.join(test_path,'testConfig.yml')

	filename = config_dict['config_file']
	output_path = config_dict['output_path']

	# Open the config file and stream it in with yaml
	stream = file(filename)
	config_file = yaml.load(stream)
	stream.close()

  # Append the general info from the config file to our overall config
	config_dict = dict(config_dict.items() + config_file[0].items())

	tests = []
	test_count = len(config_file)
	counter = 1

	# Make a dictionary entry in tests for every test we read from the config
	while counter < test_count:
		variants = []
		numpes = []
		iteration = 1
		sig_figs = 1

		cfg = config_file[counter]

		test_name = cfg['name']
		par_cfg = os.path.join(test_path,cfg['par_cfg'])
		phy_cfg = os.path.join(test_path,cfg['phy_cfg'])
		out_ref = os.path.join(test_path,cfg['out_ref'])

		if cfg.has_key('variants'):
			variants = cfg['variants']

		if cfg.has_key('numpes'):
			numpes = cfg['numpes']

		if cfg.has_key('iteration'):
			iteration = cfg['iteration']

		if cfg.has_key('sig_figs'):
			sig_figs = cfg['sig_figs']

		if len(variants) == 0:
			variants.append('')
		if len(numpes) == 0:
			numpes.append(1)

		for var in variants:
			for numpe in numpes:
				# Each test starts with all the general info. May be overwritten.
				test_dict = {}
				test_dict = dict(test_dict.items() + config_dict.items())

				if cfg.has_key('iteration'):
					test_dict['iteration'] = cfg['iteration']

				if cfg.has_key('sig_figs'):
					test_dict['sig_figs'] = cfg['sig_figs']

				output_file = test_name
				if var != '':
					output_file = output_file + '_' + var
				if len(numpes) > 1:
					output_file = output_file + '_' + str(numpe)
				output_file = output_file + '.result'

				test_dict['test_name'] = test_name
				test_dict['variant'] = var
				test_dict['numpe'] = str(numpe)
				test_dict['par_cfg'] = par_cfg.replace('$V',var).replace('$P',str(numpe))
				test_dict['phy_cfg'] = phy_cfg.replace('$V',var).replace('$P',str(numpe))
				test_dict['out_ref'] = out_ref.replace('$V',var).replace('$P',str(numpe))
				test_dict['output_file'] = os.path.join(output_path, output_file)

				tests.append(test_dict)
		counter = counter + 1

	print 'Read ' + str(len(tests)) + ' tests from ' + config_dict['config_file']
	return tests
