import numpy as np 
import scipy.fftpack as sfft
import scipy.interpolate as intp
import os



def fft(data, start_freq, end_freq, srate, freq_mult = 1.0E6, usecol=-1):
	"""
	fft(*args, **kwargs): 
	Fourier transforms an input time-domain file or time-domain data array.
	A Kaiser-Bessel window with beta = 9.5 is also applied on the time-domain file.

	---------
	Required arguments are:

	time_dom: Either a valid filename/path string, or 1 or 2D numpy array.

	start_freq = Starting frequency for Fourier transform cut (in MHz default)

	end_freq = Ending frequency for FT cut (in MHz default)

	srate = Sampling rate of time-domain data in samples/second

	---------
	Optional arguments:
	freq_mult (default 1.0E6): Sets your frequency multiplier for start/end_freq. 
		Use 1.0E6 for MHz, 1.0E9 for GHz, 1.0 for Hz, etc. 

	use_col (default -1):  Selects column for intensity data when importing 2D array. 

	---
	Returns:

	2-dimensional numpy array with 2 columns; 1st column is frequencies, 
	2nd column is intensities.

	"""


	if isinstance(data,basestring):
		path = os.path.abspath(data)
		data = np.loadtxt(path,usecols=(usecol,))


	else:	
		#print data

		try: # _FUTURE_ Potential problem: If usecol > max column in data. Needs thrown exception
			data = data[:,usecol]
		except SystemError:
			pass # 1D array

	# _FUTURE_ Code in some way to let use choose which window to use?
	data = np.append(np.kaiser(len(data),9.5)*data, np.zeros(len(data)))
	data = np.column_stack((sfft.fftfreq(len(data),1.0/srate)/(freq_mult),abs(sfft.fft(data))/100.0))


	data = data[(start_freq <= data[:,0]) & (data[:,0] <= end_freq)]

	return data


def peakpick(data, threshold_min, threshold_max = 0):
	"""
	peakpick(data, threshold_min, threshold_max = 0):

	Takes input Fourier transform, and returns a peakpick of all peaks within 
	the given threshold.

	----
	Required arguments are:

	data: 2D Fourier transform array, with frequency in 1st column and intensity in 2nd.
	threshold_min: Minimum intensity cutoff for peakpick
	threshold_max (default 0): If not zero, maximum intensity cutoff for peakpick. 
	If left at 0, peakpick will have no maximum cutoff.


	----
	Returns:

	2-dimensional numpy array with 2 columns; 1st column is frequencies, 
	2nd column is intensities.
	"""

	max_check = lambda x,y,z: x >= y and x >= z
	thres_min_check = lambda x: x > threshold_min
	thres_max_check = lambda x: x < threshold_max

	if not threshold_max:
		return np.array([[data[i,0],data[i,1]] for i in range(1,len(data)-1) if 
			max_check(data[i,1], data[i-1,1], data[i+1,1]) and 
			thres_min_check(data[i,1])])
	else:
	 	return np.array([[data[i,0],data[i,1]] for i in range(1,len(data)-1) if 
	 		max_check(data[i,1], data[i-1,1], data[i+1,1]) and 
	 		thres_min_check(data[i,1]) and thres_max_check(data[i,1])])


def spline(data, resolution):
	""" 
	spline(data, resolution):

	Splines an input Fourier transform to a given resolution 

	---
	Required arguments are:

	data: Input 2-column Fourier transform, with frequency in column 0, 
		and intensity in column 1. 

	resolution: Desired splined resolution, in units of the frequency axis. 

	----
	Returns:

	2-dimensional numpy array with 2 columns; 1st column is frequencies, 
	2nd column is intensities.
	"""

	old_res = (data[-1,0]-data[0,0])/len(data)
	scale = old_res/resolution
	new_len = int(np.ceil(scale*len(data)))

	output = np.zeros((new_len,2))
	output[:,0] = np.arange(data[0,0],data[-1,0],resolution)
	output[:,1] = intp.splev(output[:,0], intp.splrep(data[:,0],data[:,1],s=0),der=0)

	return output


#if __name__ == "__main__":

	#from matplotlib import pyplot as pp

