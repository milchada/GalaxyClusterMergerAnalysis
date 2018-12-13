import numpy as np
import matplotlib.pylab as plt
from matplotlib import colors
from read_fits import open_fits, calibrate
from scipy.interpolate import interp2d

def plot(data, sim, time, angle, suffix=''):
	fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey = True)	
	plot1 = ax1.imshow(data, norm = matplotlib.colors.LogNorm(data[data>0].min(),data.max()))
	rotimage, diff = rotate_image(data, sim, angle)
	plot2 = ax2.imshow(rotimage, norm = matplotlib.colors.LogNorm(data[data>0].min(),data.max()))
	plot3 = ax3.imshow(diff, norm = matplotlib.colors.LogNorm(np.nanmax(diff)/1e3, np.nanmax(diff)))
	ax1.set_title('Obs')
	ax2.set_title(r'$\lambda\times$ Sim')
	ax3.set_title(r'$\left(\frac{\lambda\times Sim - Obs}{\sigma_{obs}}\right)^2$')
	fig.tight_layout()
	print( "Plots complete")
	fig.colorbar(plot1, ax=ax1, shrink = 0.33)
	fig.colorbar(plot1, ax=ax2, shrink = 0.33)
	fig.colorbar(plot3, ax=ax3, shrink = 0.33)
	fig.savefig('comp/%0.2f_Gyr_%d_deg%s.png' % (time, angle, suffix))
	plt.show(block=False)
	del(fig, ax1, ax2, ax3)
	gc.collect()

def plotall(simfiles, startsnap=0, endsnap = None, suffix=''):
	if endsnap == None:
		endsnap=len(simfiles)
	bs = np.load('best_offset%s.npy' % suffix)
	ba = np.load('best_angle%s.npy' % suffix)
	for filenum in xrange(startsnap, endsnap):
		file = simfiles[filenum]
		sdata, sheader = open_fits(file,0)
		sx = calibrate(sdata.shape[1],sheader,axis=1)
		sy = calibrate(sdata.shape[1],sheader,axis=2)
		sx -= sx.mean()
		sy -= sy.mean()
		print( "Snap %d read in" % filenum)

		#interpolate
		f = interp2d(sx, sy, sdata)
		binned_data = f(x,y)
		print( "Data binned")
		rolled_data = shift(binned_data, shift = bs[filenum])
		l = np.nansum(data*rolled_data/errorsq)/np.nansum(rolled_data**2 / errorsq)
		rolled_data *= l 
		if select_halo:
			rolled_data *= halomask
			rolled_data[rolled_data == 0] = np.nan
		plot(rolled_data, times[filenum], ba[filenum], suffix)
		print( "plot complete")

def plot_c2(label):
	cs = np.load('chisq.npy')
	cs[cs==0]=np.nan
	times = np.load('times.npy')
	simfiles = glob.glob('xrayprojz/*fits')
	simfiles.sort()
	snapnums = [int(file.split('proj_')[1].split('.fits')[0]) for file in simfiles]
	times = times[snapnums]
	plt.plot(times, cs, label=label)
	plt.xlabel('t (Gyr)')
	plt.ylabel(r'$\chi^2$')
	plt.tight_layout()
# GGM = Gaussian Gradient Magnitude - look up literature