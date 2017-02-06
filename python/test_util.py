def clean_dataset():
	import os
	import subprocess
	devnull = open(os.devnull, 'w')
        tidy = os.getcwd() + '/tidy'
        subprocess.call([tidy, 'water'], stdout=devnull, stderr=devnull)


def setup_dataset():
	import os
	import subprocess
	devnull = open(os.devnull, 'w')
	subprocess.call('../../utils/setup', stdout=devnull, stderr=devnull)
