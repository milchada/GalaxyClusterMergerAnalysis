############################################
# Find cold + shock fronts in X-ray images #
############################################

import numpy as np 
import matplotlib.pylab as plt
from matplotlib import cm, colors
from astropy.io import fits
from skimage.feature import canny, peak_local_max
from scipy.ndimage import gaussian_gradient_magnitude
from sklearn.cluster import KMeans
from astropy import constants

halfwidth = 256
peak_threshold = 0.9 #this really finds the right x-ray peak

def filter_edge(img, plot=False, edgecontrast = 4, edge_threshold=0, cut = False, sigma=1):
	centre = peak_local_max(img,threshold_rel=peak_threshold)
	if cut:
		if len(centre) > 1:
			#is there one peak or are there two?
			kmeans = KMeans(n_clusters=1)
			kmeans.fit(centre)
			ccentres = kmeans.cluster_centers_
			wss1 = np.mean(np.linalg.norm(centre - ccentres, axis=1))

			kmeans = KMeans(n_clusters=2)
			kmeans.fit(centre)
			cid = kmeans.predict(centre)
			ccentres = kmeans.cluster_centers_
			wss2 = np.mean(np.linalg.norm(centre[cid==0] - ccentres[0], axis=1))
			+np.mean(np.linalg.norm(centre[cid==1] - ccentres[1], axis=1))

			if wss2 < wss1: #i.e. 2 clusters better fit
				print("2 bright peaks")
				bigger_cluster = np.argmax(np.bincount(cid))
				centre = centre[cid == bigger_cluster] #subselect points around main cluster
			centre = np.mean(centre,axis = 0)
			centre = (int(centre[0]), int(centre[1]))

		if type(centre[0]) != int:
			centre = centre[0]
		if (centre[0]-halfwidth) < 0:
			left = 0
			right = 2*halfwidth
		elif (centre[0]+halfwidth) > img.shape[0]:
			right = img.shape[0]
			left = right - 2*halfwidth
		else:
			left, right = centre[0] - halfwidth, centre[0] + halfwidth
		if (centre[1]-halfwidth) < 0:
			bottom = 0
			top = 2*halfwidth
		elif (centre[1]+halfwidth) > img.shape[0]:
			top = img.shape[0]
			bottom = right - 2*halfwidth
		else:
			bottom, top = centre[1] - halfwidth, centre[1] + halfwidth
		
		imcut = img[left:right, bottom:top]
	else:
		imcut = img
		
	edge = canny(np.log10(imcut), sigma=edgecontrast, high_threshold=edge_threshold)

	if plot:
		plt.imshow(imcut, cmap = cm.viridis, norm = colors.LogNorm(imcut.max()/1e4, imcut.max()))
		plt.imshow(np.ma.masked_less(edge*imcut, imcut.max()*peak_threshold), cmap = cm.gray_r, norm = colors.LogNorm(imcut.max()/1e4, imcut.max()))
		
		time = float(file.split('proj_')[1].split('.')[0])/10.
		plt.title('%0.2f Gyr' % time)
		plt.savefig('featuremaps/sq_weighted_%0.2f.png' % time)
	print("Edges found!")
	
	pts = np.argwhere(edge)

	ggm = gaussian_gradient_magnitude(imcut[edge], sigma=sigma)

	peak = np.argmax(ggm) #can do for more points than just this one 
	return edge,  imcut, pts, peak

	