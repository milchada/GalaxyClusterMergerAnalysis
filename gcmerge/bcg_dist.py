#########################################
####### may merge with align_BCGs #######
## currently: finding BCGs in obs, sim ##
#########################################

import numpy as np
import glob, os
from astropy.io import fits
from sklearn.cluster import KMeans
from skimage.feature import peak_local_max
from matplotlib import colors
from scipy.ndimage.interpolation import rotate

def find_peak(file, ret_peaks=False, xmin=800, xmax=1200, ymin = 800, ymax = 1200, axis=None, angle=0):	
	snapnum = file.split('slice_')[1].split('.fits')[0]
	img = -1*fits.getdata(file)
	imcut = img[xmin:xmax, ymin:ymax]
	imcut = rotate(imcut, angle)
	if axis:
		axis.cla()
		axis.imshow(imcut, norm=colors.LogNorm(imcut[imcut>0].min(),imcut.max()),origin='bottom left')
		axis.text(imcut.shape[0]/2,410,snapnum, color='w')
		# axis.set_ylim(0,ymax-ymin)
	pts = peak_local_max(imcut, num_peaks=2)

	if axis:
		axis.scatter(pts[:,1],pts[:,0],s=0.03,c='r')

	pts += (xmin, ymin) #coz that's xmin, ymin of imcut
	if ret_peaks:
		return pts
	else:
		header = fits.getheader(file)
		return np.linalg.norm(pts[1] - pts[0])*header['CDELT1']