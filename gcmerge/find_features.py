import glob
import numpy as np 
import matplotlib.pylab as plt
from matplotlib import cm, colors, patches
from astropy.io import fits
from skimage.feature import canny, peak_local_max
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from sklearn.cluster import KMeans
# from scipy.ndimage import gaussian_filter

files = glob.glob('xray_sb/*')
files.sort()

pressure = glob.glob('pseudopressure/*')
pressure.sort()

halfwidth = 256
peak_threshold = 0.9 #this really finds the right x-ray peak
# smoothed_img = gaussian_filter(imcut, 1)

def filter_edge(file,plot=False, isfile=True,edgecontrast = 4, edge_threshold=0, cut = False):
	if isfile:
		img = fits.open(file)[0].data
	else:
		img = file
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
		dim = img.shape[0]
		imcut = img[int(dim/4):int(3*dim/4), int(dim/4):int(3*dim/4)]
		
	edge1 = canny(np.log10(imcut), sigma=edgecontrast, high_threshold=edge_threshold)

	if plot:
		plt.imshow(imcut, cmap = cm.viridis, norm = colors.LogNorm(imcut.max()/1e4, imcut.max()))
		plt.imshow(np.ma.masked_less(edge1*imcut, imcut.max()*peak_threshold), cmap = cm.gray_r, norm = colors.LogNorm(imcut.max()/1e4, imcut.max()))
		
		time = float(file.split('proj_')[1].split('.')[0])/10.
		plt.title('%0.2f Gyr' % time)
		plt.savefig('featuremaps/sq_weighted_%0.2f.png' % time)
	print("Edges found!")
	return edge1, imcut

from sklearn.cluster import SpectralClustering
from sklearn.preprocessing import StandardScaler

def cluster(features, nclusters=5):
	X = StandardScaler().fit_transform(features)
	model = SpectralClustering(n_clusters=nclusters, affinity='nearest_neighbors',
		assign_labels='kmeans')
	return model.fit_predict(X)

#This is by far the best clustering method I have found! but sensitive to:
#1) edge brightness cut (0.1 works well)
#2) edgecontrast (5 works well)
#3) nclusters - currently guessing as len(points)/40 

def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    keep = ~np.isnan(data.ravel())
    tbin = np.bincount(r.ravel()[keep], data.ravel()[keep])
    nr = np.bincount(r.ravel()[keep])
    radialprofile = tbin / nr
    return radialprofile 

from circle import circle

def select_wedge(image, feature, axis, color):
	center = (halfwidth, halfwidth)
	guess = np.mean(np.linalg.norm(feature - center, axis = 1))
	fit = curve_fit(circle, feature[:,0], (feature[:,1]-halfwidth)**2, p0 = [guess])
	rad = fit[0]
	yfit = np.sqrt(circle(feature[:,0], rad)) + halfwidth
	print(rad)
	#note - if you do imshow(image), you do scatter(y, x) for things to line up
	#how does this affect wedges?
	try:
		thetas = np.rad2deg(np.arcsin((feature[:,1]-center[0])/rad)) - 90 #because of limits of arcsin output
		w1 = patches.Wedge(center, 2*rad, np.nanmin(thetas), np.nanmax(thetas))
		if sum([w1.contains_point(point) for point in feature]) == 0: #because adding 180deg maintains sin and flips sign of cos
			w1 = patches.Wedge(center, 2*rad, np.nanmin(thetas)+180, np.nanmax(thetas)+180)
		X,Y=np.mgrid[0:image.shape[1],0:image.shape[0]]
		in_wedge=np.array([w1.contains_point(point) for point in np.vstack((X.ravel(),Y.ravel())).T]).reshape(image.shape)
		wedge = image*in_wedge
		wedge[wedge == 0] = np.nan
		profile = radial_profile(wedge, center)
		# print(profile)
		axis.plot(np.arange(len(profile)), profile, label = "%0.1f kpc" % rad, color = color)
		axis.vlines(rad, 0,np.nanmax(profile), color = color, linestyle = 'dashed')
		axis.set_ylim(np.nanmin(profile),np.nanmax(profile))
		axis.legend()
		print("Wedge made")
		return wedge
	except ValueError:
		print("This feature isn't an arc")

from scipy.signal import find_peaks

def profile_towards_point(img_edges, imcut):
	peak = np.empty(len(img_edges))
	peakvalue = np.empty(len(img_edges))
	for ind in range(len(img_edges)):
		point = img_edges[ind]
		try:
			if point[0] < halfwidth:
				x = np.arange(point[0] - 10,halfwidth)
			else:
				x = np.arange(halfwidth, point[0] + 10)
			slope = (point[1]-halfwidth)/float(point[0]-halfwidth)
			if not np.isinf(slope):
				y = (slope*(x-halfwidth) + halfwidth).astype(int)
				x = x[abs(y) < 512]
				y = y[abs(y) < 512]
				profile = imcut[x, y]
				r = np.sqrt((x-halfwidth)**2+(y-halfwidth)**2)
				
				peaks = find_peaks(profile)[0]
				if peaks.any():
					mostpeak=np.argmax(profile[peaks])
					plt.plot(r, profile, c='k', alpha=0.5, linewidth=.5)
					plt.scatter(r[peaks[mostpeak]],profile[peaks[mostpeak]], c='r', s=1)
					plt.scatter(np.linalg.norm(point-halfwidth), 
						profile[np.argmin(abs(r-np.linalg.norm(point-halfwidth)))], c='b', s=1)
					peak[ind] = r[peaks[mostpeak]]
					peakvalue[ind] = profile[peaks[mostpeak]]
		except IndexError:
			print("Error %d" % ind)
			continue
	return peak, peakvalue

#now i want to cluster by the location of the peak 
	#this made it worse

def fit_wedge(file, threshold=0.1, distcut = 250, mincurve=0): #really minfeaturelength should be set by angular size at z=0.2323
	img_edge, imcut = filter_edge(file, edge_threshold=threshold)
	scale = 3.3556099403329 #fits.getheader(file)['CDELT1'] #pix to kpc
	pts = np.argwhere(img_edge)
	edgelum=img_edge*imcut
	edgelum = edgelum[edgelum>0]
	distance = np.linalg.norm(pts - (halfwidth,halfwidth), axis=1)*scale
	closepts = pts[distance < distcut]
	clusters = cluster(closepts, nclusters = len(closepts)/int(150/scale)) 
								#because shock fronts known to be ~400kpc long, i'm oversplitting rn
	splinecurve = np.empty(len(closepts))
	for cl in np.unique(clusters):
		x = np.unique(closepts[:,0][clusters == cl])
		y = closepts[:,1][[np.argwhere(closepts[:,0] == xi)[0][0] for xi in x]]
		cs = CubicSpline(x,y)
		splinecurve[clusters==cl] = cs(closepts[:,0][clusters == cl], 2)
	print("Curvature from spline fit")

	plt.clf()
	plt.imshow(imcut, norm = cm.colors.LogNorm(imcut.max()/1e4, imcut.max()))
	plt.ylim(0,512) #DO THIS BEFORE PATCHES
	halfrange = np.std(splinecurve[abs(splinecurve) < 10])
	norm = colors.Normalize(-halfrange, halfrange)
	m = cm.ScalarMappable(cmap = cm.RdBu, norm=norm)
	cols = m.to_rgba(splinecurve)
	plt.scatter(closepts[:,1], closepts[:, 0], color = cols, s=0.05)
	plt.savefig('splinecurvetest.png')
	
	
for file in files[14:50]:
	try:
		fit_wedge(file)
	except (ValueError, IndexError,IOError):
		continue
