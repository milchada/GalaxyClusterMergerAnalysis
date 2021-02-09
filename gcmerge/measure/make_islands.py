##################################################
# Takes images, finds edges, sorts into features #
##################################################

from find_features import filter_edge
from islands import Pixel, mkLines, mkIslands
from skimage.feature import peak_local_max
from astropy.io import fits
import numpy as np

def make_islands(tempfile, xrayfile, mfp_kpc= 23, xmin=800, xmax=1200, ymin=800, ymax=1200, edgecontrast=8): 
	
	"""
	tempfile: Temperature Mazzotta weighted projected FITS file
	mfp_kpc:  Mean free path in kpc for the system
	"""
	resolution = fits.getheader(tempfile)['CDELT1']

	temp = fits.getdata(tempfile)[xmin:xmax, ymin:ymax] * constants.k_B.to('keV K**-1').value
	xray = fits.getdata(xrayfile)[xmin:xmax, ymin:ymax]
	
	img_edges, temp_img, pts, peak = filter_edge(temp, edgecontrast=edgecontrast, sigma=mfp_kpc/resolution)
	
	ind = np.lexsort((pts[:,0],pts[:,1])) #i.e. sorted first by y, then by x
	pixlist = [Pixel(pt) for pt in pts[ind]]
	lines = mkLines (pixlist)
	islandlist = mkIslands (lines, 3)

	xraypeak = peak_local_max(xray, min_distance = 15)[0]
	return islandlist, temp, xraypeak