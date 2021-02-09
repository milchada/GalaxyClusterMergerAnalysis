import numpy as np 
# import matplotlib
# matplotlib.use('agg')
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy.constants import k_B
import glob, gc

kB = k_B.to('keV/K').value

def plot(dir, snapnum, ax, xmin=0, xmax=None,cmap=cm.afmhot, norm=colors.Normalize(0.2, 20), field='temperature', label=False):
	tfiles = glob.glob(dir+'fitsfiles/%s/*fits' % field)
	tfile = [file for file in tfiles if str(snapnum) in file][0]
	t = fits.getdata(tfile)
	if field == 'temperature':
		t *= kB
	if field in ['alfven', 'velocity']:
		t /= 1e5 #cm/s --> km/s
	if not xmax:
		xmax = t.shape[0]
	t = t[:, xmin:xmax]
	a = ax.imshow(t, origin='lower', cmap=cmap, norm=norm)
	if label:
		plt.colorbar(a, ax=ax)
	del(t, tfiles, tfile, a)
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

def movie(dir, snapnums, title, cmap, norm, label):
	for snapnum in snapnums:
		fig, ax = plt.subplots()
		plot(dir, snapnum, ax, cmap=cmap, norm=norm, field=title, label=label)
		fig.savefig('%s_%d.png' % (title, snapnum))
		del(fig, ax)
		gc.collect()

# if __name__=="__main__":
# 	movie('', [132, 134, 138, 140,142, 144], 'temperature', cmap=cm.afmhot, norm=colors.Normalize(0.2, 20))
# 	movie('', [132, 134, 138, 140,142, 144], 'magnetic_field', cmap=cm.magma, norm=colors.LogNorm(1e-7,1e-4))
	
def mass():
	pre = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/'
	suff = 'c1=4.1/a1=0/c2=5.2/a2=1.0/vrel=1452/b=100/'
	files = [pre+'m1=5e14/1.6e14/'+suff, 
			 pre+'m1=6e14/2.1e14/'+suff, 
			 pre+'m1=7e14/2.4e14/'+suff]
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 164, ax1)
	plot(files[1], 142, ax2)
	plot(files[2], 152, ax3)
	for a in ax.flatten():
		axisticks(a, files[0]+'fitsfiles/temperature/temperature_proj_164.fits', xmin=150, xmax= 350, tickmax=300, y=False)
	axisticks(ax1, files[0]+'fitsfiles/temperature/temperature_proj_164.fits', xmin=150, xmax= 350, tickmax=300, x=False)
	return fig, ax

def massratio():
	pre = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=6e14/'
	suff = '/c1=4.1/a1=0/c2=5.2/a2=2.0/vrel=1452/b=100/'
	files = [pre+'1e14'+suff,
			 pre+'1.6e14'+suff+'rs=0.5a/',
			 pre+'2.1e14'+suff]
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 162, ax1)
	plot(files[1], 16, ax2)
	plot(files[2], 142, ax3)
	for a in ax.flatten():
		axisticks(a, files[0]+'fitsfiles/temperature/temperature_proj_162.fits', xmin=175, xmax= 375, tickmax=300, y=False)
	axisticks(ax1, files[0]+'fitsfiles/temperature/temperature_proj_162.fits', xmin=175, xmax= 375, tickmax=300, x=False)
	
	return fig, ax

def impactparam():
	basedir='/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=6e14/1.6e14/c1=4.1/a1=0/c2=5.2/a2=2.0/vrel=1452/b='
	files = [basedir+'0/',
			 basedir+'100/rs=0.5a/',
			 basedir+'250/']
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 14, ax1)
	plot(files[1], 15, ax2)
	plot(files[2], 15, ax3)
	for a in ax.flatten():
		axisticks(a, files[0]+'fitsfiles/temperature/temperature_proj_14.fits', xmin=150, xmax= 350, tickmax=300, y=False)
	axisticks(ax1, files[0]+'fitsfiles/temperature/temperature_proj_14.fits', xmin=150, xmax= 350, tickmax=300, x=False)
	return fig, ax

def vrel():
	basedir = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=5e14/1.6e14/c1=4.1/a1=0/c2=5.2/a2=1.0/vrel='
	files = [basedir+'1252/b=100/',
			 basedir+'1752/b=100/',
			 basedir+'2200/b=100/']
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 174, ax1)
	plot(files[1], 144, ax2)
	plot(files[2], 122, ax3)
	for a in ax.flatten():
		axisticks(a, files[0]+'fitsfiles/temperature/temperature_proj_154.fits', xmin=174, xmax= 350, tickmax=300, y=False)
	axisticks(ax1, files[0]+'fitsfiles/temperature/temperature_proj_154.fits', xmin=174, xmax= 350, tickmax=300, x=False)
	return fig, ax


def c1():
	pre = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=5e14/1.6e14/'
	suff = '/a1=0/c2=5.2/a2=2.0/vrel=1452/b=100/'
	files = [pre + 'c1=3.1' + suff,
			 pre + 'c1=4.1' + suff,
			 pre + 'c1=5.1' + suff]
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 162, ax1)
	plot(files[1], 162, ax2)
	plot(files[2], 162, ax3)
	return fig, ax

def rs1():
	basedir = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=6e14/1.6e14/c1=4.1/a1=0/c2=5.2/a2=2.0/vrel=1452/b=100/rs='
	files = [basedir + '0.5a',
			 basedir + 'a',
			 basedir + '2a']
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 146, ax1)
	plot(files[1], 158, ax2)
	plot(files[2], 160, ax3)
	return fig, ax

def c2():
	pre = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=5e14/1.6e14/c1=4.1/a1=0/'
	suff = '/a2=2.0/vrel=1452/b=100/'
	files = [pre + 'c2=5.2' + suff +'rs2=0.75a/',
			 pre + 'c2=6.3' + suff,
			 pre + 'c2=7.2' + suff]
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 162, ax1)
	plot(files[1], 162, ax2)
	plot(files[2], 162, ax3)

	fs = glob.glob(files[2]+'fitsfiles/photon_emissivity/*162.fits')[0] 
	f = fits.getdata(fs)
	plot(files[0], 162, ax1, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	plot(files[1], 162, ax2, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	plot(files[2], 162, ax3, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	return fig, ax

def a1():
	pre = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=5e14/1.6e14/c1=4.1/'
	suff = '/c2=5.2/a2=2.0/vrel=1452/b=100/'
	files = [pre + 'a1=0' + suff,
			 pre + 'a1=1' + suff,
			 pre + 'a1=2' + suff]
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 162, ax1)
	plot(files[1], 162, ax2)
	plot(files[2], 162, ax3)

	fs = glob.glob(files[2]+'fitsfiles/photon_emissivity/*162.fits')[0] 
	f = fits.getdata(fs)
	plot(files[0], 162, ax1, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	plot(files[1], 162, ax2, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	plot(files[2], 162, ax3, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	return fig, ax

def a2():
	pre = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=5e14/1.6e14/c1=4.1/a1=0/c2=5.2/'
	suff = '/vrel=1452/b=100/'
	files = [pre + 'a2=0.0' + suff,
			 pre + 'a2=2.0' + suff +'rs2=0.75a/',
			 pre + 'a2=4.0' + suff]
	fig, ax = plt.subplots(ncols=3)
	ax1, ax2, ax3 = ax.flatten()
	plot(files[0], 160, ax1)
	plot(files[1], 160, ax2)
	plot(files[2], 160, ax3)


	fs = glob.glob(files[2]+'fitsfiles/photon_emissivity/*160.fits')[0] 
	f = fits.getdata(fs)
	plot(files[0], 160, ax1, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	plot(files[1], 160, ax2, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	plot(files[2], 160, ax3, field='photon_emissivity', cmap=cm.magma, norm=colors.LogNorm(f.max()/1e4, f.max()))
	return fig, ax

def theta_phi(snapnum):
	dir = '/gpfs/gibbs/pi/nagai/uc24/a2146_gamer/Hydro/m1=6e14/2.1e14/c1=4.1/a1=0/c2=5.2/a2=1.0/vrel=1452/b=100'
	theta = [0, 14, 30, 45, 60, 90]
	phi = [0, 29, 45, 60, 90]
	fig, ax = plt.subplots(ncols=3, nrows=2)
	i = 0

	for zip(t, ls) in (theta, lss):
		for p in phi:
			im = fits.getdata('snapnum_%d_theta_%d_phi%d_temp.fits' % (snapnum, t, p))*kB
			color = m.to_rgba(p)# ax.flatten()[i].imshow(im, origin='lower', cmap=cm.afmhot, norm=colors.Normalize(0.2, 20))
			ax.plot()
			i += 1
			del(im)
			gc.collect()
	return fig, ax 
