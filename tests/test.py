import sys
import os
sys.path.append(os.path.abspath('../scripts'))
import numpy as np
import spec
import helpers
from matplotlib import pyplot as pp


if __name__ == '__main__':

	test_data = np.loadtxt('test_ft.txt')
	# pp.plot(test_data[:,0],test_data[:,1])
	# pp.show()

	# cut_list = spec.peakpick(test_data, threshold_min=0.001)[:,0]
	# width = 0.7

	# cut_spec = spec.cutspec(test_data, cut_list,width=width)

	# from matplotlib import pyplot as pp 

	# pp.plot(cut_spec[:,0],cut_spec[:,1])
	# pp.show()

	# peakpick = helpers._pp_for_cutting(test_data, 0.1)

	noise_guess_vals = np.linspace(0.0002, 0.1,50)
	for val in noise_guess_vals:
		noise_level = helpers._noise_calc(test_data, val)

		print "GUESS: %.3f   ACTUAL: %.3f" %(val,noise_level)

	
