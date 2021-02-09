############################################
# Takes features, fits arc, returns radius #
#### Measures distances between features ###
############################################

import numpy as np
from points_above_gradient import find_points_above
from fit_arcs import fit_arc
from align_bcgs import angle

def select_feature(temp_img, sb_img, island, gwidth=5, centre=None, chisqmax=10, threshold = None):
	arcfit = fit_arc(temp_img, island)
	chisq_red = arcfit[:,-1]/arcfit[:,1]
	best_fit = np.argmax(np.where(chisq_red < chisqmax))
	if not threshold:
		threshold, length, xc, yc, rad, chisq = arcfit[best_fit]
	if not centre:
		centre = (yc, xc)
	feature = find_points_above(sb_img, island, threshold, type='sb')[:,0]
	return feature, centre, rad

"this directly gives the radius of curvature of the cold front"

def standoff_distance(feature1, feature2, dx=None):
	#caveat here is I need to know which isle # is the a given feature
	dist = np.empty([len(feature1), len(feature2)])
	for i in range(len(feature1)):
		dist[i] = np.linalg.norm(feature1[i] - feature2, axis = 1)
	if dx == None:
		mp=np.argwhere(dist==dist.min())[0]
		return feature1[mp[0]], feature2[mp[1]]
	else:
		return dist.min()*dx

def upstream_shock_bcg(feature, bcg_xy, dx = 6.8):
	#manual threshold to select only highest contrast part of shock
	dist = np.linalg.norm(feature  - bcg_xy, axis = 1)
	return dist.mean()*dx

def angles(bs, us, bcg_xy):
	bs_mean = np.mean(bs, axis = 0)
	us_mean = np.mean(us, axis = 0)
	bs_to_bcg = angle(bcg_xy, bs_mean)
	us_to_bcg = angle(bcg_cy, us_mean)
	return np.rad2deg(bs_to_bcg - us_to_bcg)
