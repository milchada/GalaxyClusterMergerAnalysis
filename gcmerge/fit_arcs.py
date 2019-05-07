import glob
import numpy as np
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy import constants
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_gradient_magnitude
from scipy.interpolate import CubicSpline
# from find_features import halfwidth, filter_edge

files = glob.glob('fitsfiles/temp/*fits')
files.sort()

mfp_a2146 = 23 #mean free path in kpc
resolution = fits.getheader(files[0])['CDELT1']

def find_features(filenum,isfile=True,type='temp'):
	if isfile:
		file = files[filenum]
	else:
		file = filenum
	img_edges, img = filter_edge(file, isfile=isfile)
	pts = np.argwhere(img_edges)
	
	if type == 'temp':
		img *= constants.k_B.to('keV K**-1').value
	ggm = gaussian_gradient_magnitude(img[img_edges], sigma=mfp_a2146/resolution)

	return img_edges, img, pts, ggm, resolution 

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

def fit_arc(ax1, ax2, island, resolution, time):
	
	peak = find_points_above_contrast(island, 1)

	numpoints = int(max_feature_length_kpc/resolution)	
	norm = colors.Normalize(vmin = 1, vmax = numpoints)
	cmap = cm.seismic

	arcfit = np.empty((18, 6))
	i = -1
	for mincontrast in np.arange(.9,0,-.05):
		i += 1
		arcfit[i,0] = mincontrast
		try:
			#select n points on either side of the central point
			feature = find_points_above_contrast(island, mincontrast)[:,0]
			xdata = feature[:,1]
			ydata = feature[:,0]
			ypos = ydata[ydata > halfwidth]
			guess = np.mean(np.linalg.norm(feature - center, axis = 1))
			if len(ypos):
			#fit an arc to these points:
				fit = curve_fit(circle_ypos, xdata[ydata > halfwidth], ydata[ydata > halfwidth], p0 = [halfwidth, halfwidth, guess])
			else:
				fit = curve_fit(circle_yneg, xdata[ydata < halfwidth], ydata[ydata < halfwidth], p0 = [halfwidth, halfwidth, guess])
			if np.inf not in fit[1]: #i.e. fit impossible coz not an arc
				cx, cy, r = fit[0]
				arcfit[i, 1:5] = len(xdata), cx, cy, r
				if (cx < img.shape[0]) & (cy < img.shape[1]):
					xfit = np.arange(xdata.min(), xdata.max(), 0.1)
					xyfit = np.empty([len(xfit), 2])
					xyfit[:,0] = xfit
					if len(ypos):
						xyfit[:,1] = circle_ypos(xfit, cx, cy, r)
					else:
						xyfit[:,1] = circle_ypos(xfit, cx, cy, r)
					del(xfit)
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

def main():
	fig1, ax1 = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
	fig2, ax2 = plt.subplots(nrows = 2, ncols = 3, sharex = True, sharey = True)
	for filenum in range(3,9):
		pts, peak, resolution = find_features(filenum)
		fit_arc(ax1.flatten()[filenum-3], ax2.flatten()[filenum-3], pts, peak, resolution, filenum/10. + 1)

	ax1[1][1].set_xlabel('# points on either side of feature centre')
	ax1[0][1].set_ylabel(r'$\chi^2_{red}$')
	sm = cm.ScalarMappable(cmap=cm.seismic, norm=colors.Normalize(1, max_feature_length_kpc/resolution))
	sm.set_array([])
	fig2.colorbar(sm)
	ax2[1][1].set_xlabel('X')
	ax2[0][1].set_ylabel('Y')