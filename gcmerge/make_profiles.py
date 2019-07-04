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
from find_features import filter_edge
from fit_arcs import * 
from islands import *

mfp_a2146 = 23 #mean free path in kpc
resolution = fits.getheader(files[0])['CDELT1']

sbfiles = glob.glob('xray_sb/*fits')#fitsfiles
sbfiles.sort()
tempfiles = glob.glob('tempproj/*fits')#('fitsfiles/temp/*fits')
tempfiles.sort()

time_of_interest = 1.8
snapnum = int(10*time_of_interest)
temptimes = [file.split('_')[-1].split('.')[0] for file in tempfiles]
sbtimes   = [file.split('_')[-1].split('.')[0] for file in sbfiles]

img_edges, img,    pts,    peak,    resolution = filter_edge(tempfiles[temptimes.index(snapnum)], edgecontrast=4,isfile=True, type='temp', sigma=mfp_a2146/resolution)
sb_edges,  sb_img, sb_pts, sb_peak, resolution = filter_edge(sbfiles[sbtimes.index(snapnum)], edgecontrast=4,isfile=True, type='sb', sigma=mfp_a2146/resolution)

ind = np.lexsort((pts[:,0],pts[:,1])) #i.e. sorted first by y, then by x
pixlist = [Pixel(pt) for pt in pts[ind]]
lines = mkLines (pixlist)
islandlist = mkIslands (lines, 3)

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
