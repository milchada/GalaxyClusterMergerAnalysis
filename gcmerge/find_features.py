import numpy as np 
import matplotlib.pylab as plt
from matplotlib import cm, colors
from astropy.io import fits
from scipy import ndimage as ndi 
from skimage.feature import canny, peak_local_max
from skimage.filters import sobel
import glob

files = glob.glob('xray_sb/*')
files.sort()
edgecontrasts = (5,7)
halfwidth = 256
peak_threshold = 0.1

def main(file):
	img = fits.open(file)[0].data
	centre = peak_local_max(img,threshold_rel=peak_threshold)
	if len(centre) > 1:
		centre = np.mean(centre,axis = 0)
		centre = (int(centre[0]), int(centre[1]))

	try:
		imcut = img[centre[0]-halfwidth:centre[0]+halfwidth, centre[1]-halfwidth:centre[1]+halfwidth]
	except IndexError:
		imcut = img[centre[0][0]-halfwidth:centre[0][0]+halfwidth, centre[0][1]-halfwidth:centre[0][1]+halfwidth]
	ggm = ndi.filters.gaussian_gradient_magnitude(imcut, 1)
	edge1 = canny(np.log10(imcut), sigma=edgecontrasts[0])
	edge2 = canny(np.log10(ggm), sigma=edgecontrasts[1])
	plt.imshow(imcut, cmap = cm.viridis, norm = colors.LogNorm(imcut.max()/1e4, imcut.max()))
	plt.imshow(np.ma.masked_less(edge1*imcut, imcut.max()*peak_threshold), cmap = cm.gray_r, norm = colors.LogNorm(imcut.max()/1e4, imcut.max()))
	plt.imshow(np.ma.masked_less(edge2*imcut, imcut.max()*peak_threshold), cmap = cm.magma, norm = colors.LogNorm(ggm.max()/1e4, ggm.max()))
	time = float(file.split('proj_')[1].split('.')[0])/10.
	plt.title('%0.2f Gyr' % time)
	plt.savefig('featuremaps/sq_weighted_%0.2f.png' % time)
	print "done!"

for file in files:
	print file
	try:
		main(file)
	except (ValueError, IOError):
		print "file corrupt"
	plt.close()