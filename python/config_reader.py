def config_reader(config_dict):
	import yaml
	import os

	filename = config_dict['config_file']
	output_path = config_dict['output_path']
	test_path = config_dict['test_path']

	stream = file(filename)
	config_file = yaml.load(stream)
	stream.close()

	tests = []

	test_count = len(config_file)
	counter = 1

	while counter < test_count:
		variants = []
		numpes = []

		cfg = config_file[counter]

		test_name = cfg['name']
		par_cfg = os.path.join(test_path,cfg['par_cfg'])
		phy_cfg = os.path.join(test_path,cfg['phy_cfg'])
		out_ref = os.path.join(test_path,cfg['out_ref'])

		if cfg.has_key('variants'):
			variants = cfg['variants']

		if cfg.has_key('numpes'):
			numpes = cfg['numpes']

		if len(variants) == 0:
			variants.append('')
		if len(numpes) == 0:
			numpes.append(1)

		for var in variants:
			for numpe in numpes:
				test_dict = {}

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
				print 'Read test: ' + `test_dict`
		counter = counter + 1
	return tests
