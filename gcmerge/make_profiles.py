############################################
# Plot radial profiles across each feature #
############################################

import numpy as np 
from matplotlib.path import Path
from matplotlib import patches
from astropy.io import fits
from skimage.feature import peak_local_max
from find_features import filter_edge
from fit_arcs import fit_arc
from islands import *
from align_BCGs import align_bcgs

def make_islands(potfile, xrayfile, tempfile, bcg1_pix, bcg2_pix):
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

def select_feature(temp_img, islandlist, island, gwidth, centre=None, chisqmax=10):
	arcfit = fit_arc(temp_img, islandlist, island)
	chisq_red = arcfit[:,-1]/arcfit[:,1]
	# best_fit = np.argmin(np.gradient(chisq_red)[np.gradient(chisq_red) > 0])
	best_fit = np.argmax(np.where(chisq_red < chisqmax))
	threshold, length, xc, yc, rad, chisq = arcfit[best_fit]
	if not centre:
		centre = (yc, xc)
	feature = find_points_above(sb_img, islandlist, island, threshold, type='sb')[:,0]
	return feature, centre, rad

def make_profile(temp_img, other_img, islandlist, island, label, centre, ax, ax1, ax2, gwidth=5, chisqmax=10):
	if centre == 'xray':
		feature, centre, rad = select_feature(temp_img, islandlist, island, gwidth, centre=xraypeak, chisqmax=chisqmax)
	elif centre == 'arc':
		feature, centre, rad = select_feature(temp_img, islandlist, island, gwidth, chisqmax=chisqmax)

	color = cm.gist_rainbow(island/float(len(islandlist)))
	ax.scatter(feature[:,1], feature[:,0], color=color,s=2, lw=0, label='%d' % (rad*dx))

	theta = np.rad2deg(np.arctan2(feature[:,0] - centre[0], feature[:,1] - centre[1]))
	theta[theta < 0] += 360
	# rad = np.mean(np.linalg.norm(feature - centre, axis=1)) #radius of curvature

	#THIS IS WRONG. I think it's because I stupidly redefined rad instead of using the one from the fit. 
	#but also the angular dependence seems wrong. see if that goes away with correct r
	w1 = patches.Wedge((centre[1], centre[0]),rad+45, theta.min(), theta.max(), width= 60,color='w', alpha=0.4)
	ax.add_patch(w1)

	X,Y=np.mgrid[0:temp_img.shape[1],0:temp_img.shape[0]]
	points = np.vstack((X.ravel(),Y.ravel())).T

	path = Path(w1.properties()['path'].vertices)
	grid = path.contains_points(points).reshape(temp_img.shape)
	wedge = temp_img*grid.T
	wedge[wedge == 0] = np.nan
	profile = radial_profile(wedge, centre)
	try:
		startind = np.argwhere(np.isnan(profile)).max()
	except ValueError:
		startind = 0
	profile = profile[startind:]
	ax1.plot(np.arange(len(profile)), profile,label=label, color=color, lw=2)
	ymin, ymax = ax1.get_ylim()
	# ax1.vlines(rad - startind, ymin, ymax, color=color, linestyle='dotted')
	
	wedge = other_img*grid.T
	wedge[wedge == 0] = np.nan
	profile = radial_profile(wedge, centre)
	profile = profile[startind:]
	ax2.plot(np.arange(len(profile)), profile,label=label, color=color, lw=2)
	ymin, ymax = ax2.get_ylim()
	# ax2.vlines(rad - startind, ymin, ymax, color=color, linestyle='dotted')	

	#ok now actually run it this way for a feature

def standoff_distance(islandlist, isle1, isle2, dx=None, gwidth=5, chisqmax=5):
	feature1, c1, r1 = select_feature(temp_img, islandlist, isle1, gwidth=gwidth, chisqmax=chisqmax)
	feature2, c2, r2 = select_feature(temp_img, islandlist, isle2, gwidth=gwidth, chisqmax=chisqmax)
	dist = np.empty([len(feature1), len(feature2)])
	for i in range(len(feature1)):
		dist[i] = np.linalg.norm(feature1[i] - feature2, axis = 1)
	if dx == None:
		mp=np.argwhere(dist==dist.min())[0]
		return feature1[mp[0]], feature2[mp[1]]
	else:
		return dist.min()*dx
