############################################
# Plot radial profiles across each feature #
############################################

import glob
import numpy as np 
import matplotlib.pylab as plt
from matplotlib.path import Path
from matplotlib import patches
from astropy.io import fits
from skimage.feature import peak_local_max
from find_features import filter_edge
from fit_arcs import * 
from islands import *

def make_islands(filenum, bcg1_pix, bcg2_pix):
	potfile = potfiles[filenum]
	xrayfile = xrayfiles[filenum]
	tempfile = tempfiles[filenum]
	
	mfp_a2146 = 23 #mean free path in kpc
	resolution = fits.getheader(tempfile)['CDELT1']

	temp_aligned = align_bcgs(potfile, tempfile, bcg1_pix, bcg2_pix)
	temp_aligned *= constants.k_B.to('keV K**-1').value
	xray_aligned = align_bcgs(potfile, xrayfile, bcg1_pix, bcg2_pix)
	img_edges, temp_img,   pts,   peak = filter_edge(temp_aligned, 
		edgecontrast=4, sigma=mfp_a2146/resolution)
	sb_edges,  sb_img, sb_pts, sb_peak = filter_edge(xray_aligned, 
		edgecontrast=4, sigma=mfp_a2146/resolution)

	ind = np.lexsort((pts[:,0],pts[:,1])) #i.e. sorted first by y, then by x
	pixlist = [Pixel(pt) for pt in pts[ind]]
	lines = mkLines (pixlist)
	islandlist = mkIslands (lines, 3)

	xraypeak = peak_local_max(sb_img, min_distance = 15)[0]
	return islandlist, temp_img, sb_img, xraypeak

def make_profile(temp_img, other_img, islandlist, island,label,centre,ax1, ax2):
	if centre == 'xray':
		centre = xraypeak
	elif centre == 'arc':
		arcfit = fit_arc(temp_img, islandlist, island)
		#mincontrast, len_feature, xc, yc, rad, resid
		chisq_red = arcfit[:,-1]/arcfit[:,1]
		best_fit = np.argmin(np.gradient(chisq_red)[np.gradient(chisq_red) > 0])
		threshold, length, xc, yc, rad, chisq = arcfit[best_fit]
		centre = (xc, yc)

	feature = find_points_above_contrast(temp_img, islandlist, island, threshold)[:,0]

	theta = np.rad2deg(np.arctan2(feature[:,0] - centre[0], feature[:,1] - centre[1]))
	theta[theta < 0] += 360
	rad = np.mean(np.linalg.norm(feature - centre, axis=1))
	w1 = patches.Wedge((centre[1], centre[0]), 1.2*rad, theta.min(), theta.max(), color='w', alpha=0.4)

	X,Y=np.mgrid[0:temp_img.shape[1],0:temp_img.shape[0]]
	points = np.vstack((X.ravel(),Y.ravel())).T

	path = Path(w1.properties()['path'].vertices)
	grid = path.contains_points(points).reshape(temp_img.shape)
	wedge = temp_img*grid.T
	wedge[wedge == 0] = np.nan
	profile = radial_profile(wedge, centre)
	ax1.plot(np.arange(len(profile)), profile,label=label)
	
	wedge = other_img*grid.T
	wedge[wedge == 0] = np.nan
	profile = radial_profile(wedge, centre)
	ax2.plot(np.arange(len(profile)), profile,label=label)
	return ax1, ax2
	#ok now actually run it this way for a feature

