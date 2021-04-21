from paramtests_plot import *

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

def plotic(file, ax, clr, linestyle = 'solid'):
	ax1, ax2, ax3 = ax.flatten()
	rs = file.split('rs')[1].split('rvir')[0]
	rc = file.split('rc')[1].split('rvir')[0]
	label = r'$r_s$ = %sa, $r_c$ = %s a' % (rs, rc)
	f = h5py.File(file)['fields']
	r = f['radius'].value*u.cm.to('kpc')
	t = f['temperature'].value*kB
	n = f['density'].value/m_p
	k = t*pow(n, -2./3)
	ax1.loglog(r, n, label = label, c=clr, linestyle=linestyle)
	ax2.plot(r, t, c=clr, linestyle=linestyle)
	ax3.loglog(r, k, c=clr, linestyle=linestyle)
	
def ic_compare(cluster=1, alpha='1.5'):
	if cluster == 1:
		files = glob.glob('sNFW_profile_m3.3e*a%s*' % alpha)
	elif cluster == 2:
		files = glob.glob('sNFW_profile_m9.0e*a%s*' % alpha)
	files.sort() 
	fig, ax = plt.subplots(ncols=3, sharex=True)
	
	clrs = cm.tab10
	i = 0
	for file in files:
		if ('rc0.1rvir' in file):
			plotic(file,ax, clrs(i))
			i += 1
	for file in files:
		if ('rs0.6' in file) and ('rc0.01' in file):
			plotic(file, ax, clrs(2), linestyle= 'dotted')
		if ('rs0.6' in file) and ('rc0.05' in file):
			plotic(file, ax, clrs(2), linestyle='dashed')

	h, l = ax[0].get_legend_handles_labels()
	ax[0].legend(h,l)
	for a in ax.flatten(): a.set_xlim(1,1e3)
	ax[0].set_ylabel(r'$n_g (cm^{-3})$')
	ax[1].set_ylabel(r'$T (keV)$')
	ax[2].set_ylabel(r'$K (keV cm^{2})$')
	return fig, ax 