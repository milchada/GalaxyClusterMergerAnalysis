#########################################
####### may merge with align_BCGs #######
## currently: finding BCGs in obs, sim ##
#########################################

import numpy as np
from astropy.io import fits
from skimage.feature import peak_local_max
from scipy.ndimage.interpolation import rotate

def find_peak(file, ret_peaks=False, xmin=800, xmax=1200, ymin = 800, ymax = 1200, angle=0):	
	snapnum = file.split('slice_')[1].split('.fits')[0]
	img = -1*fits.getdata(file)
	imcut = img[xmin:xmax, ymin:ymax]
	imcut = rotate(imcut, angle)
	pts = peak_local_max(imcut, num_peaks=2)

	pts += (xmin, ymin) #coz that's xmin, ymin of imcut
	if ret_peaks:
		return pts
	else:
		header = fits.getheader(file)
		return np.linalg.norm(pts[1] - pts[0])*header['CDELT1']