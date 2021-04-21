import numpy as np 
# import matplotlib
# matplotlib.use('agg')
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy.constants import k_B, m_p
import glob, gc
from scipy.ndimage import gaussian_gradient_magnitude as ggm
from astropy import units as u

kB = k_B.to('keV/K').value
m_p = m_p.to('g').value

cmaps = {'temperature':cm.afmhot, 'photon_emissivity':cm.magma, 'rm':cm.jet, 'beta':cm.inferno_r, 'ggm':cm.magma, 'magnetic_field':cm.viridis}

dirs = ['beta=inf/', 'beta=200/', 'beta=100/', 'beta=50/']

def get_tfile(dir, snapnum, field='temperature'):
	if field =='ggm':
		tfiles = glob.glob(dir+'fitsfiles/photon_emissivity/*fits')
	else:
		tfiles = glob.glob(dir+'fitsfiles/%s/*fits' % field)
	# print(tfiles)
	tfile = [file for file in tfiles if str(snapnum) in file][0]
	return tfile

def plot(ax, tfile=None,dir='', snapnum='', xmin=0, xmax=1,cmap=cm.afmhot, norm=colors.Normalize(0.2, 20), field='temperature', label=False, ggm_sigma=10):
	if not len(tfile):
		tfile = get_file(dir, snapnum, field)
	t = fits.getdata(tfile)
	if field == 'temperature':
		t *= kB
	if field in ['alfven', 'velocity']:
		t /= 1e5 #cm/s --> km/s
	if field == 'ggm':
		g = ggm(t, ggm_sigma)
		t = g/t 
		del(g); gc.collect()
	if field =='ggm':
		norm = norm(t)
	xmin = int(xmin*t.shape[0])
	xmax = int(xmax*t.shape[0])
	t = t[xmin:xmax, xmin:xmax]
	a = ax.imshow(t, origin='lower', cmap=cmap, norm=norm)#, interpolation='nearest')
	if label:
		plt.colorbar(a, ax=ax)
	del(t, a)
	gc.collect()


def axisticks(ax, tfile, xmin, xmax, x=True, y=True, tickmax=400, ntix=5):
	xtix = np.linspace(-tickmax, tickmax, ntix)
	dx = fits.getheader(tfile)['CDELT1']
	xpos = xtix/dx + (xmin+xmax)/2.
	if x:
		ax.set_xticks(xpos)
		ax.set_xticklabels(['%d' % x for x in xtix])
	if y:
		ax.set_yticks(xpos)
		ax.set_yticklabels(['%d' % x for x in xtix])

def plotall(dir, snapnum, fields=['temperature','photon_emissivity', 'rm', 'beta'], xmin=0.25, xmax=0.75, label=True):

	for field in fields:
		plt.clf()
		dfile = plot(dir, snapnum, plt.gca(), xmin, xmax, cmap = cmaps[field], field=field, norm = norm[field], label=label)
		axisticks(plt.gca(), dfile, 0, plt.gca().get_xlim()[1])
		plt.savefig(dir+field+str(snapnum)+'.png')

# def movie(dir, snapnums, title, cmap, norm, label):
# 	for snapnum in snapnums:
# 		fig, ax = plt.subplots()
# 		plot(dir, snapnum, ax, cmap=cmap, norm=norm, field=title, label=label)
# 		fig.savefig('%s_%d.png' % (title, snapnum))
# 		del(fig, ax)
# 		gc.collect()

#change these norms for different beta
norms = {'temperature':colors.Normalize(0,22),
		'photon_emissivity': colors.LogNorm(1e-21, 3e-17),
		'ggm':lambda t: colors.LogNorm(t.max()/1e5, t.max()),
		'beta':colors.LogNorm(20,1000),
		'rm':colors.Normalize(-1000,1000),
		'magnetic_field':colors.LogNorm(1e-7,2e-5)}

def evolution():
	snaps = [142,146,150,154]
	dir = 'beta=100/'
	fig, ax = plt.subplots(ncols = 4, nrows = 4, figsize=(8,8))
	fields = ['photon_emissivity', 'temperature', 'magnetic_field', 'beta']
	for i in range(len(fields)):
		for j in range(len(snaps)):
			plot(get_tfile(dir, snaps[j], fields[i]), ax = ax[i][j], field=fields[i], cmap=cmaps[fields[i]], norm = norms[fields[i]], xmin=0.25, xmax=0.75)

def allbeta():
	snapnum=156
	dirs = ['beta=200/', 'beta=100/', 'beta=50/']
	fig, ax = plt.subplots(ncols = 4, nrows = 3, figsize=(8,6))
	fields = ['photon_emissivity', 'temperature', 'magnetic_field', 'beta']
	for i in range(len(dirs)):
		for j in range(len(fields)):
			plot(get_tfile(dirs[i], snapnum, fields[j]), ax = ax[i][j], field=fields[j], cmap=cmaps[fields[j]], norm = norms[fields[j]], xmin=0.25, xmax=0.75)

def mhd_vs_turb():
	dirs = ['beta=inf/', 'beta=50/', 'beta=50/turnoff_at_0_7Gyr/']
	snapnum = 156
	fig, ax = plt.subplots(ncols = 3, nrows = 2, figsize = (8,4))
	fields = ['photon_emissivity', 'temperature']
	for i in range(len(fields)):
		for j in range(len(dirs)):
			plot(get_tfile(dirs[i], snapnum, fields[j]), ax = ax[i][j], field=fields[i], cmap=cmaps[fields[i]], norm = norms[fields[i]], xmin=0.25, xmax=0.75)

def rms():
	snaps = [150,154,160]
	dirs = ['beta=200/', 'beta=100/', 'beta=50/']
	fig, ax = plt.subplots(ncols = 3, nrows = 3, figsize=(8,8))
	for i in range(len(dirs)):
		for j in range(len(snaps)):
				plot(get_tfile(dirs[i], snaps[j], 'rm'), ax = ax[i][j], field='rm', cmap=cmaps['rm'], norm = norms['rm'], xmin=0.25, xmax=0.75)

def rmrot():
	from astropy.units import cm 
	snapnum = 150
	dir = 'beta=100/fitsfiles/rm/'
	files = glob.glob(dir+'*150*theta*'); files.sort()
	fig, ax = plt.subplots(ncols = 3, nrows = 2, figsize=(6,4))
	for i in range(len(files)):
		ax.flatten()[i+1].imshow(fits.getdata(files[i])*cm.to('kpc'), origin='lower', cmap=cmaps['rm'], norm=norms['rm'])

def dip_profiles():
	basedir = '/Users/Mila/Documents/Research/2015-2021 Thesis/Paper3-Abell2146-MHD/FITS/'
	dx = 0.48828125
	betas = ['200','100','50']
	fig, ax = plt.subplots(nrows = 3,ncols=2, sharex=True)
	fields = ['beta', 'kT', 'pressure', 'sb', 'rho']
	colors = ['tab:blue', 'tab:green', 'tab:orange']
	for i in range(len(fields)):
		for j in range(len(betas)):
			file = basedir+'b'+betas[j]+'_'+fields[i]+'.dat'
			a = np.genfromtxt(file)
			x = a[:,0]
			y = a[:,1]*(dx**2)
			ax.flatten()[i].plot(x,y,color=colors[j], label=r'$\beta$ = %s' % betas[j])
	for j in range(len(betas)):
		a = np.genfromtxt(basedir+'b'+betas[j]+'_kT.dat')
		x = a[:,0]
		t = a[:,1]*(dx**2)
		b = np.genfromtxt(basedir+'b'+betas[j]+'_rho.dat')
		n = b[:,1]*(dx**2)/m_p 
		k = t*pow(n, -2./3)
		ax.flatten()[-1].plot(x,k,color=colors[j], label=r'$\beta$ = %s' % betas[j])