import yt, os, glob
import numpy as np

def rotate_proj(ds, snapnum, theta=0, phi=0, width=(1,'Mpc')):
	t = np.deg2rad(theta)
	p = np.deg2rad(phi)
	n = [np.sin(t)*np.cos(p), -np.sin(t)*np.sin(p), np.cos(t)]
	s = yt.FITSOffAxisProjection(ds, normal=n, fields=('gas', 'temperature'), center='c', weight_field='mazzotta_weighting', width=width)
	s.writeto('snapnum_%d_theta_%d_phi%d_temp.fits' % (snapnum, theta, phi))

def compute(dir, snapnum, thetas=[0, 15, 30, 45, 60, 90], phis=[0, 15, 30, 45, 60, 90]):
	ds = yt.load(basedir+dir+'/Data_0000'+str(snapnum))
	for theta in thetas:
		for phi in phis:
			rotate_proj(ds, snapnum, theta, phi)


from astropy.constants import k_B
from astropy.io import fits
import matplotlib, gc
matplotlib.use('agg')
import matplotlib.pylab as plt 
from matplotlib import cm, colors
from scipy.ndimage import rotate
kB = k_B.to('keV/K').value

phi15 = glob.glob('*145*phi14*fits')
phi30 = glob.glob('*145*phi29*fits')
phi45 = glob.glob('*145*phi45*fits')
phi15.sort(); phi30.sort(); phi45.sort()

def plot(phi0, phi45, phi90, ncols=3, rot1=0, rot2=0, rot3=0):
	fig, ax = plt.subplots(nrows=3, ncols=ncols, sharex=True, sharey=True)
	for i in range(ncols):
		rot0 = rotate(fits.getdata(phi0[i])*kB, rot1, reshape=False)
		rot45 = rotate(fits.getdata(phi45[i])*kB, rot2, reshape=False)
		rot90 = rotate(fits.getdata(phi90[i])*kB, rot3, reshape=False)
		ax[0][i].imshow(rot0, origin='lower', cmap=cm.afmhot, norm = colors.Normalize(0.2,20))
		ax[1][i].imshow(rot45, origin='lower', cmap=cm.afmhot, norm = colors.Normalize(0.2,20))
		ax[2][i].imshow(rot90, origin='lower', cmap=cm.afmhot, norm = colors.Normalize(0.2,20))
		del(rot0, rot45, rot90)
		gc.collect()

	return fig, ax