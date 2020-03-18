##################################################
# Takes images, finds edges, sorts into features #
##################################################

from align_bcgs import align_bcgs
from find_features import filter_edge
from islands import Pixel, mkLines, mkIslands
from skimage.feature import peak_local_max
from astropy.io import fits
import numpy as np

def make_islands(obsfile, errfile, potfile, xrayfile, tempfile, bcg1_pix=None, bcg2_pix=None, mfp_kpc= 23): 
	
	"""
	potfile:  Potential slice FITS file
	xrayfile: Xray SB emission weighted projected FITS file
	tempfile: Temperature Mazzotta weighted projected FITS file
	bcg1_pix: Pixel position of BCG1 - pixel positions of potential minima
	bcg2_pix: Pixel position of BCG2
	mfp_kpc:  Mean free path in kpc for the system
	"""

	resolution = fits.getheader(tempfile)['CDELT1']

	temp_aligned = align_bcgs(obsfile, errfile, potfile, tempfile, bcg1_pix, bcg2_pix)
	temp_aligned *= constants.k_B.to('keV K**-1').value
	xray_aligned = align_bcgs(obsfile, errfile, potfile, xrayfile, bcg1_pix, bcg2_pix)
	img_edges, temp_img,   pts,   peak = filter_edge(temp_aligned, edgecontrast=4, sigma=mfp_kpc/resolution)
	sb_edges,  sb_img, sb_pts, sb_peak = filter_edge(xray_aligned, edgecontrast=4, sigma=mfp_kpc/resolution)

	ind = np.lexsort((pts[:,0],pts[:,1])) #i.e. sorted first by y, then by x
	pixlist = [Pixel(pt) for pt in pts[ind]]
	lines = mkLines (pixlist)
	islandlist = mkIslands (lines, 3)

	xraypeak = peak_local_max(sb_img, min_distance = 15)[0]
	return islandlist, temp_img, sb_img, xraypeak