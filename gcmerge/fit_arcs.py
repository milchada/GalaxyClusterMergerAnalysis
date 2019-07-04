##############################
# Fit arcs to sharp features #
##############################

import glob
import numpy as np
from matplotlib import pylab as plt, colors, cm
from scipy import optimize
import gc

def calc_R(x,y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f(c, x, y):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()

def leastsq_circle(x,y):
    # coordinates of the barycenter
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
    xc, yc = center
    Ri       = calc_R(x, y, *center)
    R        = Ri.mean()
    residu   = np.sum((Ri - R)**2)
    return xc, yc, R, residu

def fit_arc(island, time,ax1=None, ax2=None):
	cmap = cm.seismic
	i = -1
	arcfit = np.empty((18,6))
	for mincontrast in np.arange(.9,0,-.05):
		i += 1
		arcfit[i,0] = mincontrast
		try:
			#select n points on either side of the central point
			feature = find_points_above_contrast(island, mincontrast)[:,0]
			xdata = feature[:,1]
			ydata = feature[:,0]
			arcfit[i, 1] = len(feature)
			fit = leastsq_circle(xdata, ydata)[:3]
			arcfit[i, 2:5] = fit[:3]
			#sum of distances of points from fit
			arcfit[i, 5] = fit[-1]
			del(fit, xdata, ydata, feature)
			gc.collect()
			print(mincontrast, " done")
		
		except IndexError:
			print(mincontrast, " not enough pts")
			continue
		except (RuntimeError):
			print(mincontrast, "no good fit")
			continue
	return arcfit

def radial_profile(data, center):
	"""Make these emission weighted"""
	y, x = np.indices((data.shape))
	r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
	r = r.astype(np.int)
	keep = ~np.isnan(data.ravel())
	tbin = np.bincount(r.ravel()[keep], data.ravel()[keep])
	nr = np.bincount(r.ravel()[keep])
	radialprofile = tbin / nr
	return radialprofile 

