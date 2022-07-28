import numpy as np 
# import matplotlib
# matplotlib.use('agg')
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy.constants import k_B, m_p
import glob, gc
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
	tfile = [file for file in tfiles if '_'+str(snapnum) in file][0]
	return tfile

def plot(ax, tfile=None,dir='', snapnum='', xmin=0, xmax=1, field='temperature', label=False, ggm_sigma=10):
	cmap = cmaps[field]
	norm = norms[field]
	if not tfile:
		tfile = get_tfile(dir, snapnum, field)
	t = fits.getdata(tfile)
	if field == 'temperature':
		t *= kB
	if field in ['alfven', 'velocity']:
		t /= 1e5 #cm/s --> km/s
	if field == 'ggm':
		g = ggm(t, ggm_sigma)
		t = g/t 
		norm = colors.LogNorm(t.min(), t.max())
		del(g); gc.collect()
	xmin = int(xmin*t.shape[0])
	xmax = int(xmax*t.shape[0])
	t = t[xmin:xmax, xmin:xmax]
	a = ax.imshow(t, origin='lower', cmap=cmap, norm=norm)#, interpolation='nearest')
	if label:
		plt.colorbar(a, ax=ax)
	del(t, a)
	gc.collect()


def axisticks(ax, tfile, xmin, xmax, x=True, y=True, tickmax=400, ntix=5, fontsize=14):
	xtix = np.linspace(-tickmax, tickmax, ntix)
	dx = fits.getheader(tfile)['CDELT1']
	xpos = xtix/dx + (xmin+xmax)/2.
	if x:
		ax.set_xticks(xpos)
		ax.set_xticklabels(['%d' % x for x in xtix], fontsize=fontsize)
	if y:
		ax.set_yticks(xpos)
		ax.set_yticklabels(['%d' % x for x in xtix], fontsize=fontsize)

def movie(dir, snapnums, title, cmap, norm, label):
	for snapnum in snapnums:
		fig, ax = plt.subplots()
		plot(dir, snapnum, ax, cmap=cmap, norm=norm, field=title, label=label)
		fig.savefig('%s_%d.png' % (title, snapnum))
		del(fig, ax)
		gc.collect()

#change these norms for different beta
norms = {'temperature':colors.Normalize(2,18),
		'photon_emissivity': colors.LogNorm(2e3, 5e6),
		'ggm': colors.LogNorm(3e-23, 3e-18),
		'beta':colors.LogNorm(20,1000),
		'rm':colors.Normalize(-10000,10000),
		'magnetic_field':colors.LogNorm(1e-7,2e-5)}


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
				plot(tfile=get_tfile(dirs[i], snaps[j], 'rm'), ax = ax[i][j], field='rm', cmap=cmaps['rm'], norm = norms['rm'], xmin=0.25, xmax=0.75)

def rmrot():
	from astropy.units import cm 
	snapnum = 150
	dir = 'beta=100/fitsfiles/rm/'
	files = glob.glob(dir+'*150*theta*'); files.sort()
	fig, ax = plt.subplots(ncols = 3, nrows = 2, figsize=(6,4))
	for i in range(len(files)):
		ax.flatten()[i+1].imshow(fits.getdata(files[i])*cm.to('kpc'), origin='lower', cmap=cmaps['rm'], norm=norms['rm'])

def dip_profiles():
	basedir = '/Users/Mila/Documents/Research/PhD/Paper3/dips/sector_'
	dx = 0.48828125
	betas = ['200','100','50', 'inf']
	
	fields = ['tproj', 'sb']#, 'beta', pressure', 'entropy',  'rho']
	labels = [r'$\beta$', 'kT (keV)', r'$\epsilon (10^{5} cm^{-2}s^{-1})$',r'$\rho_g (10^{-26} g/cm^{-3}$)', 
				r'P$_{pseudo} (\sqrt{\rm{SB}}\times$ kT)', r'K$_{pseudo} ($ kT \times {\rm SB}^{-1/3}$)']
	colors = ['tab:blue', 'tab:green', 'tab:orange', 'tab:red']
	fig, ax = plt.subplots(nrows = 3,ncols=2, sharex=True)

	
	xmin = np.zeros(len(betas))
	for j in range(len(betas)):
		# sb = np.genfromtxt(basedir+'beta'+betas[j]+'_sb.dat')[:,1] *(dx**2)
		t = np.genfromtxt(basedir+'beta'+betas[j]+'_temp.dat')[:,1]*kB *(dx**2)
		# k = t* sb**(-1./3) #sb ~ rho^2
		xmin[j] = np.argmin(t[:40]) #K flattening

	xmin -= xmin[-1]

	for i in range(len(fields)):
		for j in range(len(betas)):
			file = basedir+'beta'+betas[j]+'_'+fields[i]+'.dat'
			a = np.genfromtxt(file)
			x = a[:,0] - xmin[j]
			y = a[:,1]*(dx**2)
			if fields[i] == 'temp':
				y *= kB
			ax.flatten()[i].plot(x,y,color=colors[j], label=r'$\beta$ = %s' % betas[j])

	for j in range(len(betas)):
		sb = np.genfromtxt(basedir+'beta'+betas[j]+'_sb.dat')[:,1] *(dx**2)
		t = np.genfromtxt(basedir+'beta'+betas[j]+'_temp.dat')[:,1]*kB *(dx**2)
		p = np.sqrt(sb) * t
		k = t* sb**(-1./3) #sb ~ rho^2
		x = np.genfromtxt(basedir+'beta'+betas[j]+'_temp.dat')[:,0] - xmin[j]
		ax[1][0].plot(x,p,color=colors[j])
		ax[1][1].plot(x,k,color=colors[j])
	i=0
	for a in ax.flatten(): 
		a.vlines(32, a.get_ylim()[0], a.get_ylim()[1], color='k', linestyle='dashed')
		a.vlines(75, a.get_ylim()[0], a.get_ylim()[1], color='k', linestyle='dotted')
		a.vlines(43, a.get_ylim()[0], a.get_ylim()[1], color='k', linestyle='dotted')

		a.set_ylabel(labels[i])
		i += 1
	ax[0][0].set_yscale('log')
	plt.xlim(2,35)

	plt.tight_layout()

from scipy.ndimage import gaussian_filter
def unsharp_mask(img, sigma=1, contrast = 2):
	smooth = gaussian_filter(img, sigma)
	return contrast*img - (contrast-1)*smooth


def beta_proj(ds, filename, width=(2, 'Mpc'), ret=False):
	def rhosq(field, data): 
		return data['gas','density']* data['gas','density'] 
	 ds.add_field(("gas","rhosq"), rhosq, units='g**2/cm**6', force_override=True,sampling_type='cell') 
	pth = yt.FITSProjection(ds, 'z', ('gas','pressure'), width=width, weight_field='rhosq',center=('max','emission_measure')) 
	pmag = yt.FITSProjection(ds, 'z', ('gas','magnetic_pressure'), width=width, weight_field='rhosq',center=('max','emission_measure')) 
	bproj = pth['pressure'].data/pmag['magnetic_pressure'].data 
	hdu = fits.PrimaryHDU(bproj) 
	hdu.writeto(filename, overwrite=True)
	if ret:
		return bproj

