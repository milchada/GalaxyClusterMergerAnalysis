############################################
# Plot radial profiles across each feature #
############################################

import numpy as np 
from matplotlib.path import Path
from matplotlib import patches
from astropy.io import fits
from measure_feature import select_feature
from matplotlib import colors, cm

def radial_profile(data, center):
	"""These have the same weighting as the FITS projection"""
	y, x = np.indices((data.shape))
	r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
	r = r.astype(int)
	keep = ~np.isnan(data.ravel())
	tbin = np.bincount(r.ravel()[keep], data.ravel()[keep])
	#counts number of points at a given radius, weighted by the temperature at that point
	nr = np.bincount(r.ravel()[keep])
	#counts number of points at each radius, not weighted
	radialprofile = tbin / nr
	#ratio of these two gives the profile. yes makes sense.
	return radialprofile 

def make_profile(temp_img, island, label, ax, centre = 'arc', ax1=None, gwidth=5, chisqmax=10, ret=True):
	if centre == 'xray':
		feature, centre, rad = select_feature(temp_img, island, gwidth, centre=xraypeak, chisqmax=chisqmax)
	elif centre == 'arc':
		feature, centre, rad = select_feature(temp_img, island, gwidth, chisqmax=chisqmax)

	ax.scatter(feature[:,1], feature[:,0],s=2, lw=0, label='%d' % (rad*dx))

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
	if ax1:
		ax1.plot(np.arange(len(profile)), profile,label=label, color=color, lw=2)
		ymin, ymax = ax1.get_ylim()
	if ret:
		return profile 
