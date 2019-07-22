#########################################################
### Compare simulated and observed X-ray image arrays ###
#########################################################

import numpy as np
import glob, gc 
from scipy.interpolate import interp2d
from scipy.ndimage.interpolation import rotate, shift
from multiprocessing import Pool
import matplotlib.pylab as plt
from read import open_fits, calibrate, init
from align_BCGs import align_bcgs

data, errorsq, x, y, xp = init(obsfile, errfile)   

def shift_image(image, offset):
	xstep, ystep, startshift = offset
	imagecut = np.roll(np.roll(image, xstep-startshift, axis=0), ystep-startshift, axis=1)
	imagecut[np.isnan(imagecut)] = 0 #otherwise rotation gets fucked
	return imagecut

def shift_compare(image, offset=(0,0,0), plot=False, time=None):
	shiftimage = shift_image(image, offset)
	masked = np.ma.masked_equal(shiftimage, 0)
	diff = (shiftimage - data)**2 / errorsq
	chisq = np.nansum(diff[~masked.mask])
	print("offset (%d, %d) complete" % (xstep, ystep))
	if plot:
		plot(data, shiftimage, diff, time)
	del(shiftimage, shiftimage)
	gc.collect()
	return chisq

def match(simfiles, potfiles, filenum, xraypeak, peak_threshold = 0.7, pixels_to_shift = 10, suffix=''):
	chisqs = np.empty(len(simfiles))
	
	file = simfiles[filenum]
	potfile = potfiles[filenum]

	rot_image = align_bcgs(potfile, file)
	binned_data = rot_image
	binned_data[binned_data == 0] = np.nan

	#optimise rescaling within circle 
	l = np.nansum(data*binned_data/errorsq)/np.nansum(binned_data**2 / errorsq)
	rot_image *= l
	
	if filenum > filemin: 
		chisqs_around_peak = np.zeros([pixels_to_shift, pixels_to_shift])
		for xstep in range(pixels_to_shift):
			for ystep in range(pixels_to_shift):
				chisqs_around_peak[xstep,ystep] = shift_compare(rot_image,offset=(xstep, ystep, pixels_to_shift/2))
				print("xstep = %d, ystep = %d done!" % (xstep, ystep))

		minx, miny = np.argmin(chisqs_around_peak)[0]
		chisqs[filenum] = shift_compare(rot_image, (minx, miny, pixels_to_shift//2))
		print("done!")

	else:
		print("halos not merged yet. done!")
		chisqs[filenum] = shift_compare(rot_image)

	np.save('chisq%s' % suffix, chisqs)

def compare():
	pool = Pool(10)
	print(len(simfiles))
	for filenum in range(len(simfiles)):
		pool.apply_async(match, args=(simfiles, potfiles, filenum, xraypeak))
