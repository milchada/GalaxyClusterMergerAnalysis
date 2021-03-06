import numpy as np 
# import matplotlib
# matplotlib.use('agg')
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy.constants import k_B
import glob, gc
from scipy.ndimage import gaussian_gradient_magnitude as ggm

kB = k_B.to('keV/K').value

cmaps = {'temperature':cm.afmhot, 'photon_emissivity':cm.magma, 'rm':cm.jet, 'beta':cm.inferno_r, 'ggm':cm.magma, 'magnetic_field':cm.viridis}

dirs = ['beta=inf/', 'beta=200/', 'beta=100/', 'beta=50/']

def plot(dir, snapnum, ax, xmin=0, xmax=1,cmap=cm.afmhot, norm=colors.Normalize(0.2, 20), field='temperature', label=False, ggm_sigma=10):
	if field =='ggm':
		tfiles = glob.glob(dir+'fitsfiles/photon_emissivity/*fits')
	else:
		tfiles = glob.glob(dir+'fitsfiles/%s/*fits' % field)
	# print(tfiles)
	tfile = [file for file in tfiles if str(snapnum) in file][0]
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
	del(t, tfiles, a)
	gc.collect()
	return tfile

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

def fig1():
	snaps = [142,146,150,154]
	dir = 'beta=100/'
	fig, ax = plt.subplots(ncols = 4, nrows = 4, figsize=(8,8))
	fields = ['photon_emissivity', 'temperature', 'magnetic_field', 'beta']
	for i in range(len(fields)):
		for j in range(len(snaps)):
			plot(dir, snaps[j], ax = ax[i][j], field=fields[i], cmap=cmaps[fields[i]], norm = norms[fields[i]], xmin=0.25, xmax=0.75)

def fig2():
	snapnum=156
	dirs = ['beta=200/', 'beta=100/', 'beta=50/']
	fig, ax = plt.subplots(ncols = 4, nrows = 3, figsize=(8,6))
	fields = ['photon_emissivity', 'temperature', 'magnetic_field', 'beta']
	for i in range(len(dirs)):
		for j in range(len(fields)):
			plot(dirs[i], snapnum, ax = ax[i][j], field=fields[j], cmap=cmaps[fields[j]], norm = norms[fields[j]], xmin=0.25, xmax=0.75)

def fig3():
	dirs = ['beta=inf/', 'beta=50/', 'beta=50/turnoff_at_0_7Gyr/']
	snapnum = 156
	fig, ax = plt.subplots(ncols = 3, nrows = 2, figsize = (8,4))
	fields = ['photon_emissivity', 'temperature']
	for i in range(len(fields)):
		for j in range(len(dirs)):
			plot(dirs[j], snapnum, ax = ax[i][j], field=fields[i], cmap=cmaps[fields[i]], norm = norms[fields[i]], xmin=0.25, xmax=0.75)


def fig4():
	snaps = [142,146,150,154]
	dirs = ['beta=200/', 'beta=100/', 'beta=50/']
	fig, ax = plt.subplots(ncols = 4, nrows = 3, figsize=(8,6))
	for i in range(len(dirs)):
		for j in range(len(snaps)):
			try:
				plot(dirs[i], snaps[j], ax = ax[i][j], field='rm', cmap=cmaps['rm'], norm = norms['rm'], xmin=0.25, xmax=0.75)
			except IndexError:
				continue
