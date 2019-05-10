import glob
import numpy as np
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy import constants
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_gradient_magnitude

mfp_a2146 = 23 #mean free path in kpc
resolution = fits.getheader(files[0])['CDELT1']

halfwidth = 256
center = (halfwidth, halfwidth)
max_feature_length_kpc = 500

def circle_ypos(x, xm, ym, r):
	t = np.arccos((x - xm)/r) #t in [0, pi]
	t[(x - xm)/r > 1] = np.arccos(1)
	t[(x - xm)/r < -1] = np.arccos(-1)
	return r*np.sin(t) + ym

def circle_yneg(x, xm, ym, r):
	t = np.arccos((x - xm)/r) #t in [0, pi]
	t[(x - xm)/r > 1] = np.arccos(1)
	t[(x - xm)/r < -1] = np.arccos(-1)
	return - r*np.sin(t) + ym

def fit_half(fn, xdata, ydata, p0):
	fit = curve_fit(fn, xdata, ydata, p0 = p0)
	xfit = np.arange(xdata.min(), xdata.max(), 1)
	xyfit = np.empty([len(xfit), 2])
	if np.inf not in fit[1]: #i.e. fit impossible coz not an arc
		cx, cy, r = fit[0]
		xyfit[:,0] = xfit
		xyfit[:,1] = fn(xfit, cx, cy, r)
	else:
		xyfit[:,1] = np.nan
		cx, cy, r = np.nan
		print "arc is a bad fit"
	return xyfit, cx, cy, r 

def fit_arc(ax1, ax2, island, resolution, time):
	peak = find_points_above_contrast(island, 1)
	cmap = cm.seismic
	i = -1
	arcfit = np.empty((18,6))
	for mincontrast in np.arange(.9,0,-.05):
		i += 1
		arcfit[i,0] = mincontrast
		try:
			#select n points on either side of the central point
			feature = find_points_above_contrast(island, mincontrast)[:,0]
			guess = np.mean(np.linalg.norm(feature - center, axis = 1))
			xdata = feature[:,1]
			ydata = feature[:,0]
			if len(ydata[ydata > halfwidth]):
			#fit an arc to these points:
				posfit, cx, cy, r = fit_half(circle_ypos, xdata[ydata > halfwidth], ydata[ydata > halfwidth], p0 = [halfwidth, halfwidth, guess])
				xyfit = posfit
			if len(ydata[ydata < halfwidth]):
				negfit, cxn, cyn, rn = fit_half(circle_yneg, xdata[ydata < halfwidth], ydata[ydata < halfwidth], p0 = [halfwidth, halfwidth, guess])
			if len(ydata[ydata > halfwidth]) & len(ydata[ydata < halfwidth]):
				xyfit = np.vstack((posfit, negfit))
				cx = (cx + cxn)/2.
				cy = (cy + cyn)/2.
				r = (r + rn)/2.
			elif len(ydata[ydata < halfwidth]):
				xyfit = negfit
				cx, cy, r = cxn, cyn, rn
			arcfit[i, 1:5] = len(feature), cx, cy, r 
			#sum of distances of points from fit

			leastdist = np.array([np.min(np.linalg.norm(xyfit - pt,axis=1)) for pt in feature])
			arcfit[i, 5] = sum(leastdist)/len(leastdist)
			ax1.scatter(mincontrast, sum(leastdist)/len(leastdist), c = 'k', marker = 'x')
			ax2.scatter(cy, cx, c = cmap(mincontrast), lw=0)
			ax1.set_title('t = %0.1f Gyr' % time)
			ax2.set_title('t = %0.1f Gyr' % time)
			print(mincontrast, " done")
		except IndexError:
			print(mincontrast, " not enough pts")
			continue
		except (RuntimeError, TypeError):
			print(mincontrast, "no good fit")
			continue
	return arcfit

def test_arcfits(island):
	fig, ax = plt.subplots(ncols=2)
	arcfit = fit_arc(ax[0], ax[1], island, resolution, 1.5)
	plt.close()
	fig, ax = plt.subplots()
	sm = cm.ScalarMappable(cmap=cm.seismic, norm=colors.Normalize(0,1))
	for row in range(len(arcfit)):
		cx, cy, r = arcfit[row, 2:5]
		feature = find_points_above_contrast(island, arcfit[row,0])[:,0]
		xdata = feature[:,1]
		ydata = feature[:,0]
		xpos = xdata[ydata > halfwidth]
		xneg = xdata[ydata < halfwidth]
		if len(xpos):
			xfit = np.arange(xpos.min(),xpos.max(),.1)
			yfit_pos = circle_ypos(xfit, cx, cy, r)
			plt.scatter(yfit_pos, xfit, c = sm.to_rgba(arcfit[row,0]), lw=0)
		if len(xneg):
			xfit = np.arange(xneg.min(),xneg.max(),.1)
			yfit_neg = circle_yneg(xfit, cx, cy, r)
			plt.scatter(yfit_neg, xfit, c = sm.to_rgba(arcfit[row,0]), lw=0)
	plt.scatter(ydata, xdata, marker='x', c='k')
	sm.set_array([])
	fig.colorbar(sm)
	plt.xlim(100,500)
	plt.ylim(100,500)

names = name = {1:'Cold Front', 3: 'Bow', 6: 'Swirl', 9: 'Upstream'}
def arcfits_stability(islandnums=[1,3,6,9],names=names, time=1.5):
	fig1, ax1 = plt.subplots()
	for island in islandnums:
		fig, ax = plt.subplots(ncols= 2)
		arcfit = fit_arc(ax[0],ax[1],island, resolution, 1.5)
		ax1.plot(arcfit[:,0], arcfit[:,-1], label=names[island])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('Contrast threshold')
	ax1.set_ylabel(r'$\Sigma d_{min}$')
	return fig1, ax1