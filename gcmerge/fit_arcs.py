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

	peak = np.argmax(ggm) #can do for more points than just this one 
	return img_edges, img, pts, peak, resolution 

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

def fit_arc(ax1, ax2, pts, peak, resolution, time):
	
#this does the job but looks so ugly
#is there really no simpler way?

	numpoints = int(max_feature_length_kpc/resolution)	
	norm = colors.Normalize(vmin = 1, vmax = numpoints)
	cmap = cm.seismic

	for n in range(1, numpoints):
		try:
			#select n points on either side of the central point
			feature = pts[peak - n:peak + n + 1]
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

					leastdist = np.linalg.norm(xyfit - feature[i], axis = 1)
					ax1.scatter(n, sum(leastdist)/len(leastdist), c = 'k', marker = 'x')
					ax2.scatter(cy, cx, c = cmap(norm(n)))
					ax1.set_title('t = %0.1f Gyr' % time)
					ax2.set_title('t = %0.1f Gyr' % time)
					print(n, " done")
		except (RuntimeError, TypeError):
			print(n, "no good fit")
			continue

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