#BCG velocities

import yt
from astropy.io import fits
from bcg_dist import find_peak

snapnum = 10
filename = 'Data_0000%d' % snapnum
potfile = 'fitsfiles/potential/potential_slice_%d.fits' % snapnum
ds = yt.load(filename)
p1, p2 = find_peak(potfile, ret_peaks = True)
p1 = np.insert(p1, 2, 0)
p2 = np.insert(p2, 2, 0)
dx = fits.getheader(potfile)['CDELT1']
c1 = p1*dx/1000.
c2 = p2*dx/1000.

#now use c1, c2 to select a certain # of particles around each pot min
#they're currently ~ 1Mpc apart so should be safe to draw a circle around each
s1 = ds.sphere(center=c1, radius=(100,'kpc'))
s2 = ds.sphere(center=c2, radius=(100,'kpc'))
#now get the particle IDs within these

v1 = s1.ds.all_data('particle_velocity').to('km/s')
v2 = s2.ds.all_data('particle_velocity').to('km/s')
v1 = np.mean(v1, axis = 1)
v2 = np.mean(v2, axis = 1)
