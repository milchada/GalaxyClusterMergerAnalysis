#########################################
####### may merge with align_BCGs #######
## currently: finding BCGs in obs, sim ##
#########################################

import numpy as np
import glob, os
from astropy.io import fits
from sklearn.cluster import KMeans
from skimage.feature import peak_local_max
import matplotlib.pylab as plt
from matplotlib import colors

def find_peak(file, axis=None, ret=True, xmin=400, xmax=600, ymin = 400, ymax = 600, min_distance=3, threshold_rel=0.01):	
	snapnum = file.split('slice_')[1].split('.fits')[0]
	img = -1*fits.getdata(file)
	imcut = img[xmin:xmax, ymin:ymax]
	if axis:
		axis.cla()
		axis.imshow(imcut, norm=colors.LogNorm(imcut[imcut>0].min(),imcut.max()),origin='bottom left')
		axis.text(imcut.shape[0]/2,410,snapnum, color='w')
		# axis.set_ylim(0,ymax-ymin)
	pts = peak_local_max(imcut, min_distance=min_distance, threshold_rel=threshold_rel)
	if axis:
		axis.scatter(pts[:,1],pts[:,0],s=0.03,c='r')

	pts += (xmin, ymin) #coz that's xmin, ymin of imcut
	
	if ret:
		return pts 

def bcg_separation(file, distmax_kpc=470, ret_peaks=False, xmin=300, xmax=1200, ymin = 300, ymax = 1200, axis=None):	
	data = fits.getdata(file)
	header = fits.getheader(file)
	md = 3
	maxsep = distmax_kpc/header['CDELT1']
	oldlen = len(data)**2 
	while md < maxsep:
		peaks = find_peak(file, axis=axis, min_distance = md, xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax) #x1, y1, x2, y2 in pix  
		md += 1 
		if len(peaks) == 2: 
			break 
		elif len(peaks) == oldlen: 
			break 
		oldlen = len(peaks) 
	
	if len(peaks) == 1:
		print( "cluster relaxed")
	else:
		if len(peaks) > 2: 
			kmeans = KMeans(n_clusters=2)
			kmeans.fit(peaks)
			cid = kmeans.predict(peaks)
			peaks = (kmeans.cluster_centers_).astype(int)
		
		if ret_peaks:
			return peaks
		else:
			return np.linalg.norm(peaks[1] - peaks[0])*header['CDELT1']