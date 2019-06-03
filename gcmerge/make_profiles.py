import glob
import numpy as np 
import matplotlib.pylab as plt
from matplotlib.path import Path
from matplotlib import cm, colors, patches
from astropy.io import fits
from skimage.feature import canny, peak_local_max
from sklearn.cluster import KMeans
from astropy import constants
from scipy.ndimage import gaussian_gradient_magnitude

halfwidth = 256
peak_threshold = 0.9 #this really finds the right x-ray peak
# smoothed_img = gaussian_filter(imcut, 1)

def filter_edge(file,plot=False, isfile=True,edgecontrast = 2, edge_threshold=0, cut = False):
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

files = glob.glob('tempproj/*fits')#('fitsfiles/temp/*fits')
files.sort()

mfp_a2146 = 23 #mean free path in kpc
resolution = fits.getheader(files[0])['CDELT1']

def find_features(filenum, edgecontrast=4, isfile=True, type='temp'):
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

files = glob.glob('tempproj/*fits')#('fitsfiles/temp/*fits')
files.sort()
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
