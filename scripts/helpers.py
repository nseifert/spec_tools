import numpy as np
import spec as sp



def _pp_for_cutting(data, threshold_min):
	max_check = lambda x: data[x,1] >= data[x-1,1] and data[x,1] >= data[x+1,1]
	thres_min_check = lambda x: data[x,1] > threshold_min

	return np.array([[data[i,0],data[i,1], i] for i in range(1,len(data)-1) if 
		max_check(i) and thres_min_check(i)])


def _linewidth(spec, transition, noise_level):
	idx = int(transition[2])
	c_low = 0
	c_up = 0

	chk = lambda x: spec[x,1] >= noise_level
	x_chk = lambda x,y: spec[x,1] > spec[y,1]
	bound_chk = lambda x: x > spec.shape[0]

	while chk(idx-c_low):
		if x_chk(idx-c_low-1,idx-c_low) and x_chk(idx-c_low+1, idx-c_low):
			break
		c_low += 1

	while chk(idx+c_up):
		if bound_chk(idx+c_up):
			c_up = (spec.shape[0]-1) - idx
			break
			
		if x_chk(idx+c_up-1, idx+c_up) and x_chk(idx+c_up+1, idx+c_up):
			break

	
	return np.array([idx-c_low,idx+c_up],dtype=int)



def _noise_calc(spec, noise_level, spline_res=0.002):
	spectrum = sp.spline(spec, spline_res)
	peakpick = _pp_for_cutting(spectrum,noise_level)

	statlist = []

	for i in range(len(peakpick)-1):
		c1 = _linewidth(spectrum,peakpick[i],noise_level)[1]
		c2 = _linewidth(spectrum,peakpick[i+1],noise_level)[0]
		#print "C1 %d C2 %d" % (c1,c2)

		if spec[c2,0] > spec[c1,0]:
			statlist.extend([spec[j,1] for j in range(c1,c2)])

	avg = np.average(np.array(statlist))
	std = np.std(np.array(statlist))

	return np.sqrt(2.0) * np.sqrt(avg**2 + std**2)


def _cut_fixed_width(spec,cut_list,width):

	mask = np.ones(spec.shape[0])

	res = spec[1,0]- spec[0,0] 

	for val in cut_list:

		low_bound = val - width/2
		high_bound = val + width/2

		idx_low = np.floor((low_bound-spec[0,0])/res) + 1
		idx_high = np.floor((high_bound-spec[0,0])/res) - 1

		mask[idx_low:idx_high] = 0

	spec[:,1] = spec[:,1] * mask
	return spec

def _cut_to_noise(spec,cut_list,noise_level):

	# First calculate noise level
	pass
	#real_noise_lvl = __noise_calc(spec,noise_level)



