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
			# if ax1:
			# 	ax1.scatter(mincontrast, sum(leastdist)/len(leastdist), c = 'k', marker = 'x')
			# 	ax2.scatter(cy, cx, c = cmap(mincontrast), lw=0)
			# 	ax1.set_title('t = %0.1f Gyr' % time)
			# 	ax2.set_title('t = %0.1f Gyr' % time)
			print(mincontrast, " done")
		except IndexError:
			print(mincontrast, " not enough pts")
			continue
		except (RuntimeError):
			print(mincontrast, "no good fit")
			continue
	return arcfit

def test_arcfits(island, showpts=True):
	arcfit = fit_arc(island, 1.5)
	fig, ax = plt.subplots()
	sm = cm.ScalarMappable(cmap=cm.seismic, norm=colors.Normalize(0,1))
	if showpts:
		feature = find_points_above_contrast(island, 0)[:,0]
		xdata = feature[:,1]
		ydata = feature[:,0]
		plt.scatter(ydata, xdata, marker='x', c='k')
	for row in range(len(arcfit)):
		cx, cy, r = arcfit[row, 2:5]
		feature = find_points_above_contrast(island, arcfit[row,0])[:,0]
		xdata = feature[:,1]
		ydata = feature[:,0]
		thetas = np.arctan2((ydata - cy),(xdata - cx))
		xfit = r*np.cos(thetas) + cx
		yfit = r*np.sin(thetas) + cy
		plt.scatter(yfit, xfit, c = sm.to_rgba(arcfit[row,0]), lw=0)
	sm.set_array([])
	fig.colorbar(sm)
	plt.xlim(300,600)
	plt.ylim(300,600)

names = name = {1:'Cold Front', 3: 'Bow', 6: 'Swirl', 9: 'Upstream'}
def arcfits_stability(islandnums=[1,3,6,9],names=names, time=1.5):
	fig1, ax1 = plt.subplots()
	for island in islandnums:
		fig, ax = plt.subplots(ncols= 2)
		arcfit = fit_arc(island, 1.5)
		ax1.plot(arcfit[:,0], arcfit[:,-1], label=names[island])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('Contrast threshold')
	ax1.set_ylabel(r'$\Sigma d_{min}$')
	return fig1, ax1