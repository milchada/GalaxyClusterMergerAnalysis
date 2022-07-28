#########################################
####### may merge with align_BCGs #######
## currently: finding BCGs in obs, sim ##
#########################################

import numpy as np
from astropy.io import fits
from skimage.feature import peak_local_max
from scipy.ndimage.interpolation import rotate

def find_peak(file, ret_peaks=False, xmin=800, xmax=1200, ymin = 800, ymax = 1200, angle=0):	
	#read in FITS slice of projection
	snapnum = file.split('slice_')[1].split('.fits')[0]
	#potential is always negative, while this algorithm finds the maximum, so make everything positive
	img = -1*fits.getdata(file)
	
	#if you look far from the cluster center, you may catch other local peaks, so focus on the center
	#the default values above work for a FITS image of size 2048x2048. Check the size of the image, select central ~1Mpc
	imcut = img[xmin:xmax, ymin:ymax]
	
	#just in case you want to work with rotated images, the coordinates of the minima will of course move
	imcut = rotate(imcut, angle)
	pts = peak_local_max(imcut, num_peaks=2)

	pts += (xmin, ymin) #coz that's xmin, ymin of imcut
	if ret_peaks:
		return pts
	else:
		#this just gives the separation between peaks in kpc
		#use this to match the BCG separation in observation
		header = fits.getheader(file)
		return np.linalg.norm(pts[1] - pts[0])*header['CDELT1']