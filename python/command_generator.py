def command_generator(config_dict):
	import os
	tests = config_dict['tests']

	data_path = config_dict['data_path']
	charmrun = os.path.relpath(config_dict['charmrun'],data_path)
	exe = os.path.relpath(config_dict['exe'],data_path)

	for test_dict in tests:
		numpe = test_dict['numpe']
		par_cfg = os.path.relpath(test_dict['par_cfg'],data_path)
		phy_cfg = os.path.relpath(test_dict['phy_cfg'],data_path)

		command = [charmrun, '++local', '+p'+numpe, exe, par_cfg, phy_cfg]

		test_dict['command'] = command
	return config_dict



