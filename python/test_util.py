def clean_dataset():
	import os
	import subprocess
	devnull = open(os.devnull, 'w')
	subprocess.call('../../utils/tidy', stdout=devnull, stderr=devnull)

def setup_dataset():
	import os
	import subprocess
	devnull = open(os.devnull, 'w')
	subprocess.call('../../utils/setup', stdout=devnull, stderr=devnull)
