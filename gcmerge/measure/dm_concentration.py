import numpy as np 
from astropy.io import fits
from matplotlib import patches
from matplotlib.path import Path
from astropy import cosmology
from sklearn.cluster import KMeans
from skimage.feature import peak_local_max
from scipy.optimize import curve_fit
from astropy import units

def find_peak(file, distmax_kpc=470, ret_peaks=False, xmin=300, xmax=1200, ymin = 300, ymax = 1200):	
	snapnum = file.split('slice_')[1].split('.fits')[0]
	img = -1*fits.getdata(file)
	imcut = img[xmin:xmax, ymin:ymax]
	pts = peak_local_max(imcut, num_peaks=2)

	if axis:
		axis.scatter(pts[:,1],pts[:,0],s=0.03,c='r')

	pts += (xmin, ymin) #coz that's xmin, ymin of imcut
	if ret_peaks:
		return pts
	else:
		header = fits.getheader(file)
		return np.linalg.norm(pts[1] - pts[0])*header['CDELT1']

def radial_profile(data, center):
	"""These have the same weighting as the FITS projection"""
	y, x = np.indices((data.shape))
	r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
	r = r.astype(int)
	keep = ~np.isnan(data.ravel())
	tbin = np.bincount(r.ravel()[keep], data.ravel()[keep])
	#counts number of points at a given radius, weighted by the temperature at that point
	nr = np.bincount(r.ravel()[keep])
	#counts number of points at each radius, not weighted
	radialprofile = tbin / nr
	#ratio of these two gives the profile. yes makes sense.
	return radialprofile 

cosmo = cosmology.default_cosmology.get()
rho_crit = cosmo.critical_density(.2323).to('Msun kpc**-3').value

potfile = 'fitsfiles/potential/potential_slice_13.fits'
rhofile = 'fitsfiles/rhoproj/rhoproj_13.fits'

c1, c2  = find_peak(potfile, ret_peaks = True)
dx = fits.getheader(potfile)['CDELT1']
rad = np.linalg.norm(c1 - c2)*dx
c1 = yt.YTArray((c1[1]*dx,c1[0]*dx,7000.), 'kpc')	
c2 = yt.YTArray((c2[1]*dx,c2[0]*dx,7000.), 'kpc')

#lokas 2002
def sigma(R, c, rv):
	"""
	R, rv in kpc
	"""
	Rp = R/rv 
	Rp[Rp == 0] = dx/(2*rv)
	gc = 1/(np.log(1+c) - c/(1+c))
	
	Mv = 4*np.pi*(rv**3)*200*rho_crit/3 #Msun

	f1 = (c*Rp)**2 - 1

	rhonorm = c**2 * gc * Mv /(2*np.pi * rv**2)
	rhoproj = np.zeros(R.shape)

	rs = rv/c 

	# print(len(Rp), len(f1))
	rhoproj[R > rs] = rhonorm*np.abs(1 - np.arccos(1/(c*Rp[R>rs]))/np.sqrt(np.abs(f1[R>rs])))/f1[R>rs]**2
	rhoproj[R < rs] = rhonorm*np.abs(1 - np.arccosh(1/(c*Rp[R<rs]))/np.sqrt(np.abs(f1[R<rs])))/f1[R<rs]**2
	singularity = np.argmin(abs(R - rs))
	# # print singularity
	rhoproj[singularity] = rhonorm/3.
	
	return rhoproj * units.Msun.to('g')/(units.kpc.to('cm')**2)#msun/kpc**2


def jointfit_rhoproj(rhofile, potfile, xmin=800, xmax=1200):
	p1, p2  = find_peak(potfile, ret_peaks = True)
	dx = fits.getheader(potfile)['CDELT1']
	rho = fits.getdata(rhofile)[xmin:xmax, xmin:xmax]
	p1 -= xmin
	p2 -= xmin

	x = np.arange(rho.shape[0])
	y = np.arange(rho.shape[1])
	X, Y = np.meshgrid(x, y)
	xy = np.vstack((X.ravel()*dx, Y.ravel()*dx))

	def sigma_sum(xy, c1, c2, rv1, rv2):
		x, y = xy
		r1 = np.linalg.norm((x - p1[1]*dx, y - p1[0]*dx), axis = 0) #in kpc
		r2 = np.linalg.norm((x - p2[1]*dx, y - p2[0]*dx), axis = 0)
		return sigma(r1, c1, rv1) + sigma(r2, c2, rv2)

	pfit, pcov = curve_fit(sigma_sum, xy, rho.ravel(), p0 = [3, 6, 1200, 200])
	return pfit, pcov #c, rvir

def singlefit_rhoproj(rhofile, potfile, xmin=800, xmax=1200):
	p1, p2  = find_peak(potfile, ret_peaks = True)
	dx = fits.getheader(potfile)['CDELT1']
	rho = fits.getdata(rhofile)[xmin:xmax, xmin:xmax]
	p1 -= xmin
	p2 -= xmin
	rad = 0.8*np.linalg.norm(p1 - p2).astype(int) #80% of distance between BCGs so profiles don't overlap too much
	def fit(p, rad, rho):
		c = patches.Circle((p[1],p[0]), rad)
		X,Y=np.mgrid[0:rho.shape[1],0:rho.shape[0]]
		points = np.vstack((X.ravel(),Y.ravel())).T
		path = Path(rad*(c.properties()['path'].vertices) + (p[1], p[0]))
		grid = path.contains_points(points).reshape(rho.shape)
		cutout = rho*grid.T
		cutout[cutout==0] = np.nan
		profile = radial_profile(cutout, (p2[1],p2[0]))
		r = np.arange(len(profile))*dx

		return curve_fit(sigma, r, profile, p0 = [3, 1200])[0]

	 fit1 = fit(p1, rad, rho)
	 fit2 = fit(p2, rad, rho)
	 return fit1, fit2

def nfw(r, rhoc, rs):
    return rhoc/(r/rs * (1+ r/rs)**2)

def Mvir(rho0, c, rs):
	rho0 *= units.g.to('Msun')/(units.cm.to('kpc')**3)
	return 4*np.pi*(rs**3)*rho0*(np.log(1+c) - c/(1+c))

file = 'Data_000013'

def density_profile(file, center, rad, filename=None):
	ds = yt.load(file)
	sph = ds.sphere(center, rad)
	prof = yt.create_profile(sph, ["radius"], n_bins=100,fields=["all_density"], weight_field='cell_volume')
	if filename:
		prof.save_as_dataset(filename)
	return prof.field_data.values(), prof.x_bins

def fit_density(filename, file=None, center=None, potfile=None, rmax = 500):
	if file:
		dens, r = density_profile(file, center, rad, filename)
	else:
		rho = h5py.File(filename+'.h5')
		r = rho['data']['radius'].value
		dens = rho['data']['all_density'].value
	
	r *=units.cm.in_units('kpc')
	keep = (r < rmax)*(dens > 0)
	rho_0, rs = curve_fit(nfw, r[keep], dens[keep])[0]
	rvir = r[np.argmin(abs(dens - 200*rho_crit))]
	c = rvir/rs
	return rho_0, c, rs

