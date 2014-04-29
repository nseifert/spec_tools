import numpy as np 
import scipy.fftpack as sfft
import scipy.interpolate as intp
import os
from helpers import *


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

	max_check = lambda x: data[x,1] >= data[x-1,1] and data[x,1] >= data[x+1,1]
	thres_min_check = lambda x: data[x,1] > threshold_min
	thres_max_check = lambda x: data[x,1] < threshold_max

	if not threshold_max:
		return np.array([[data[i,0],data[i,1]] for i in range(1,len(data)-1) if 
			max_check(i) and thres_min_check(i)])
	else:
	 	return np.array([[data[i,0],data[i,1]] for i in range(1,len(data)-1) if 
	 		max_check(i) and thres_min_check(i) and thres_max_check(i)])


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




def cutspec(spec,cut_list,**kwargs):
	"""
	cutspec(spec,cut_list,width):

	Cuts a set of frequencies from an input Fourier transformed spectrum. 
	Assumes spectrum is equally spaced across entire frequency range.

	---
	Required arguments are:

	spec: Input 2-column Fourier transform, either as a filename (string)
	or 2D ndarray

	cut_list: List or 1D ndarray of frequencies to cut. 

	cut_to_noise (default False): Determines if frequencies will be cut with
	a fixed frequency window, or "smart cut" to the noise level.

		if cut_to_noise is False, you must supply cutspec() with the variable
		"width":

		width (float): Specified frequency window size for transition cuts, in units
		of the frequency axis of input FT. 

		if cut_to_noise is True, you must supply cutspec() with the variable 
		"noise_lvl":

		noise_lvl: Estimated noise level in units of input FT's intensity axis.


	-- 
	Notes: If the user decides to use cut_to_noise, noise_lvl must be a reasonable
	estimate of the noise level. However, cutspec() will call helper functions
	to statistically ascertain the noise level. Using cut_to_noise is advantageous,
	as it will dynamically allocate the frequency window to cut on a line-to-line basis.

	"""
	if 'cut_to_noise' not in kwargs:
		cut_to_noise = False


	if isinstance(spec,basestring):
		try:
			path = os.path.abspath(spec)
		except:
			print 'cutspec() CRITICAL ERROR:'
			print 'Cannot find spectrum with specified path.'
		try:
			spec = np.loadtxt(path)
		except: 
			print 'Spectrum at specified path is malformed.'
	else:
		if not len(spec): 
			raise Exception("Spectrum is empty.")

	if not len(cut_list):
		raise Exception("Cut List is empty.")

	if cut_to_noise:

		try: 
			if not isinstance(noise_level, float):
				raise TypeError("noise_level must be a float.")
		except:
			print 'cutspec() CRITICAL ERROR:'
			print 'You did not specify the required \"noise_level\" variable to a float'  
			return None

		return __cut_to_noise(spec,cut_list,noise_level)

	if not cut_to_noise:

		try:
			width = kwargs['width']		

		except:
			print 'cutspec() CRITICAL ERROR:'
			print 'You did not specify the required "width" variable [cut_to_noise is False]'
			return None

		return __cut_fixed_width(spec,cut_list,width)




