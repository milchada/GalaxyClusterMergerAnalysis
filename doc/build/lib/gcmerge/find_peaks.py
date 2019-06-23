##identify minima in gravitational potential to compare to BCG positions

import numpy as np 
import matplotlib.pylab as plt
from matplotlib import colors
from astropy.io import fits
from skimage.feature import peak_local_max
import glob

def find_peak(file, axis=None, ret=True):	
	snapnum = file.split('slice_')[1].split('.fits')[0]
	img = -1*fits.getdata(file)
	scale = fits.getheader(file)['CDELT1'] #kpc/pix
	imcut = img[400:600,400:600]
	if axis:
		axis.cla()
		axis.imshow(imcut, norm=colors.LogNorm(imcut[imcut>0].min(),imcut.max()),origin='bottom left')
		axis.text(imcut.shape[0]/2,410,snapnum, color='w')
		axis.set_ylim(0,200)
	pts = peak_local_max(imcut, min_distance=3, threshold_rel=.01)
	pts += 400 #coz that's xmin, ymin of imcut
	
	if len(pts) == 2:
		peaks = np.concatenate((pts[0],pts[1]))
	else:
		peaks = np.concatenate((pts, pts))
	if axis:
		axis.scatter(pts[:,1],pts[:,0],s=0.03,c='r')
	if ret:
		return peaks

if __name__=="__main__":
	files = glob.glob('*fits')
	files.sort()
	fig, ax = plt.subplots(nrows=3,ncols=3,sharex=True, sharey=True,start=9,end=18)
	peaks = np.empty([end-start+1, 4])
	for axis, file in zip(ax.flatten(), files[start:end]):
		try:
			find_peak(file, axis)
		except(IOError,IndexError,ValueError):
			continue
		fig.savefig('potential_z.png')