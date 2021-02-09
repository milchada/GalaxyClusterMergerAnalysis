#BCG velocities

import numpy as np
import yt, gc, glob
from astropy.io import fits
from skimage.feature import peak_local_max
from sklearn.cluster import KMeans  

def v_sphere(ds, center, radius):
	s = ds.sphere(center=center, radius=(radius,'kpc'))
	v = s['particle_velocity'].to('km/s')#
	print(len(v), ' particles collected')
	return v.mean(axis=0)
	
def v_bcg(filename, r1 = 300, r2 = 1000, min_distance_kpc=50,width=5):
	ds = yt.load(filename)
	pot = yt.FITSSlice(ds, 'z', ('gravitational_potential'), width=(width,'Mpc')).hdulist[0]
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
	
	cx1 = pot.header['CRVAL1'] - dx*pot.header['CRPIX1']
	p1[:2] += cx1
	p2[:2] += cx1
	c1 = np.insert(np.array([p1[1],p1[0]]), 2, 7000.)/1000. #kpc to Mpc
	c2 = np.insert(np.array([p2[1],p2[0]]), 2, 7000.)/1000.

	#now use c1, c2 to select a certain # of particles around each pot min
	#they're currently ~ 1Mpc apart so should be safe to draw a circle around each
	v1 = v_sphere(ds, c1, r1)
	v2 = v_sphere(ds, c2, r1)
	v1l = v_sphere(ds, c1, r2)
	v2l = v_sphere(ds, c2, r2)

	vrel12 = v1 - v2l
	vrel21 = v2 - v1l
	return vrel12, vrel21, c1, c2

from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

files = glob.glob('Data*')
vrels = Parallel(n_jobs=num_cores)(delayed(v_bcg)(file) for file in files)
