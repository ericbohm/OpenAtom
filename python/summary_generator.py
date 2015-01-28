def summary_generator(config_dict):
	import yaml
	import checkresult
	import string

	tests = config_dict['tests']
	iteration = config_dict['iteration']
	sig_figs = config_dict['sig_figs']

	summaries = []

	for test_dict in tests:
		summary = [test_dict['test_name'],test_dict['variant'],test_dict['numpe']]

		output_file = test_dict['output_file']
		out_ref = test_dict['out_ref']

		test_result = checkresult.checkresult(output_file, out_ref, iteration, sig_figs)
		summary = summary + test_result
		summaries.append(summary)

	return summaries
