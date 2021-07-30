#BCG velocities

import numpy as np
import yt, gc, glob
from astropy.io import fits
from skimage.feature import peak_local_max
from sklearn.cluster import KMeans  
import matplotlib.pylab as plt
from matplotlib import colors, cm

def v_sphere(ds, center, radius):
	s = ds.sphere(center=center, radius=(radius,'kpc'))
	v = s['particle_velocity'].to('km/s')#
	print(len(v), ' particles collected')
	return v.mean(axis=0)
	
def v_bcg(filename, r1 = 50, min_distance_kpc=50,width=5):
	ds = yt.load(filename)
	pot = yt.FITSSlice(ds, 'z', ('gravitational_potential'), width=(width,'Mpc')).hdulist[0]
	print("potential read in")
	dx = pot.header['CDELT1']
	p = peak_local_max(-1*pot.data, min_distance=int(min_distance_kpc/dx))
	if len(p) > 2:
		kmeans = KMeans(n_clusters=2)
		fit = kmeans.fit(p)  
		p1, p2 = fit.cluster_centers_ 
	else:
		p1, p2 = p
	p1 = p1*dx
	p2 = p2*dx
	print("peaks found")
	
	cx1 = pot.header['CRVAL1'] - dx*pot.header['CRPIX1']
	p1[:2] += cx1
	p2[:2] += cx1
	c1 = np.insert(np.array([p1[1],p1[0]]), 2, 7000.)/1000. #kpc to Mpc
	c2 = np.insert(np.array([p2[1],p2[0]]), 2, 7000.)/1000.

	#now use c1, c2 to select a certain # of particles around each pot min
	#they're currently ~ 1Mpc apart so should be safe to draw a circle around each
	v1 = v_sphere(ds, c1, r1)
	v2 = v_sphere(ds, c2, r1)
	print("Done!")
	return (v1 - v2)

def n(theta, phi): 
	t = np.deg2rad(theta); p = np.deg2rad(phi)
	return np.array([np.sin(t)*np.cos(p), np.sin(t)*np.sin(p), np.cos(t)]) 

def compare(vlos_obs, filename, r1=50, min_distance_kpc=50, width=2):
	vrel = v_bcg(filename, r1, min_distance_kpc, width).value
	a = np.arange(0,95,5)
	theta, phi = np.meshgrid(a, a)
	vlos = np.dot(n(theta, phi).T, vrel)
	delta = abs(vlos - vlos_obs)/vlos_obs
	plt.imshow(delta, origin='lower', cmap=cm.viridis)
	plt.xticks(np.arange(len(a))[::3],a[::3])
	plt.yticks(np.arange(len(a))[::3],a[::3])
	plt.ylabel(r'$\theta(^\circ)$')  
	plt.xlabel(r'$\phi(^\circ)$')  
	plt.xlim(0,18); plt.ylim(0,18)
	plt.colorbar()
	plt.title(r'|$v_{los,pred} - v_{los,obs}|/v_{los, obs}$')
	fig = plt.gcf()
	ax = plt.gca()
	return fig, ax

