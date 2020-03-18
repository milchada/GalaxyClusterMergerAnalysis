############################################
# Takes features, fits arc, returns radius #
#### Measures distances between features ###
############################################

import numpy as np
from points_above_gradient import find_points_above
from fit_arcs import fit_arc
from align_bcgs import angle

def select_feature(temp_img, islandlist, island, gwidth, centre=None, chisqmax=10, threshold = None):
	arcfit = fit_arc(temp_img, islandlist, island)
	chisq_red = arcfit[:,-1]/arcfit[:,1]
	best_fit = np.argmax(np.where(chisq_red < chisqmax))
	if not threshold:
		threshold, length, xc, yc, rad, chisq = arcfit[best_fit]
	if not centre:
		centre = (yc, xc)
	feature = find_points_above(sb_img, islandlist, island, threshold, type='sb')[:,0]
	return feature, centre, rad

"this directly gives the radius of curvature of the cold front"

def standoff_distance(temp_img, islandlist, isle1, isle2, dx=None, gwidth=5, chisqmax=5):
	#caveat here is I need to know which isle # is the a given feature
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

def upstream_shock_bcg(temp_img, islandlist, isle, bcg_xy, dx = 6.8, gwidth = 5, chisqmax = 5):
	feature, c, r = select_feature(temp_img, islandlist, isle, gwidth=gwidth, chisqmax=chisqmax, threshold = 0.9)
	#manual threshold to select only highest contrast part of shock
	dist = np.linalg.norm(feature  - bcg_xy, axis = 1)
	return dist.mean()*dx

def angles(temp_img, islandlist, isle_bs, isle_us, bcg_xy, gwidth = 5, chisqmax = 5):
	bs, cb, rb = select_feature(temp_img, islandlist, isle_bs, gwidth=gwidth, chisqmax=chisqmax, threshold = 0.9)
	us, cu, ru = select_feature(temp_img, islandlist, isle_us, gwidth=gwidth, chisqmax=chisqmax, threshold = 0.9)
	#how to find the angle between three points/two lines? 
	bs_mean = np.mean(bs, axis = 0)
	us_mean = np.mean(us, axis = 0)
	
