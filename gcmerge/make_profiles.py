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
from astropy import constants
from scipy.ndimage import gaussian_gradient_magnitude
from find_features import filter_edge

mfp_a2146 = 23 #mean free path in kpc
resolution = fits.getheader(files[0])['CDELT1']

def find_features(filenum, edgecontrast=4, isfile=True, type='temp', halfwidth = 256, peak_threshold = 0.9):
	if isfile:
		file = files[filenum]
	else:
		file = filenum
	img_edges, img = filter_edge(file, edgecontrast=edgecontrast,isfile=isfile)
	pts = np.argwhere(img_edges)
	
	if type == 'temp':
		img *= constants.k_B.to('keV K**-1').value
	ggm = gaussian_gradient_magnitude(img[img_edges], sigma=mfp_a2146/resolution)

	peak = np.argmax(ggm) #can do for more points than just this one 
	return img_edges, img, pts, peak, resolution 

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

tempfiles = glob.glob('tempproj/*fits')#('fitsfiles/temp/*fits')
tempfiles.sort()
img_edges, img, pts, peak, resolution = find_features(18)

files = glob.glob('xray_sb/*fits')#fitsfiles
files.sort()
sb_edges, sb_img, sb_pts, sb_peak, resolution = find_features(18)

ind = np.lexsort((pts[:,0],pts[:,1])) #i.e. sorted first by y, then by x
pixlist = [Pixel(pt) for pt in pts[ind]]
lines = mkLines (pixlist)
islandlist = mkIslands (lines, 3)

def find_points_above_contrast(island, mincontrast=1):
	feature = islandlist[island]
	points = np.array([feature.lines[0].pixlist[0].rawx, feature.lines[0].pixlist[0].rawy])
	for line in feature.lines:
		for pix in line.pixlist:
			points = np.insert(points, -1, (pix.rawx,pix.rawy))
	points = points[1:-1]
	points = np.reshape(points, (len(points)/2, 2))
	ggm = gaussian_gradient_magnitude(img, 1)
	ggms = []
	for point in points:
		ggms.append(ggm[point[0]][point[1]])
	ggms = np.array(ggms)
	return points[np.argwhere( ggms >= mincontrast*ggms.max())]

centre = peak_local_max(sb_img, min_distance = 15)
centre = centre[0]

def main(feature,label,centre=centre):
	theta = np.rad2deg(np.arctan2(feature[:,0] - centre[0], feature[:,1] - centre[1]))
	theta[theta < 0] += 360
	rad = np.mean(np.linalg.norm(feature - centre, axis=1))
	w1 = patches.Wedge((centre[1], centre[0]), 1.2*rad, theta.min(), theta.max(), color='w', alpha=0.4)
	
	X,Y=np.mgrid[0:img.shape[1],0:img.shape[0]]
	points = np.vstack((X.ravel(),Y.ravel())).T

	path = Path(w1.properties()['path'].vertices)
	grid = path.contains_points(points).reshape(img.shape)
	wedge = img*grid.T
	wedge[wedge == 0] = np.nan
	profile = radial_profile(wedge, centre)
	plt.plot(np.arange(len(profile)), profile,label=label)
