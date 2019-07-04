#########################################################
### Compare simulated and observed X-ray image arrays ###
#########################################################

import numpy as np
import glob, gc 
from scipy.interpolate import interp2d
from scipy.ndimage.interpolation import rotate, shift
from multiprocessing import Pool
import matplotlib.pylab as plt
from read import open_fits, calibrate, init
from align_bcgs import offset_tilt

basedir = '/home/fas/nagai/uc24/scratch60/'
obsfile = basedir+'kT_out.fits'
errfile = basedir+'err_kT_out.fits'
fitsdir = 'fitsfiles/temp'
xraysb_obs = basedir+'ff.img.e300_7000_bin3_all_obsids_box_excl_point_sources.fits'
xraysb_err = basedir+'ff.img.e300_7000_bin3_all_obsids_box_0s_to_1s_no_bkg_subtr2_thresh.fits'

simfiles, times, data, errorsq, x, y, xp = init(obsfile, errfile, fitsdir)

potfiles = glob.glob('fitsfiles/potential_z/*fits')
potfiles.sort()

xraypeak = init(xraysb_obs, xraysb_err, 'fitsfiles/xray_sb', peak_only=True)

def rotate_image(image, angle):
	image[np.isnan(image)] = 0 #otherwise rotation gets fucked
	rotimage = rotate(image, angle = angle)
	datasize = image.shape[0]
	rotsize = rotimage.shape[0]

	if rotsize != datasize:
		xstart = (rotsize - datasize)//2
		rotimage = rotimage[xstart:xstart+datasize, xstart:xstart+datasize]
		rotimage[rotimage == 0] = np.nan
	return rotimage

def shift_image(image, offset):
	xstep, ystep, startshift = offset
	imagecut = np.roll(np.roll(image, xstep-startshift, axis=0), ystep-startshift, axis=1)
	imagecut[np.isnan(imagecut)] = 0 #otherwise rotation gets fucked
	return imagecut

def plot(data, rotimage, diff, time, suffix=''):
	fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey = True)	
	
	plot1 = ax1.imshow(data, norm = matplotlib.colors.LogNorm(data[data>0].min(),data.max()))
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

def shift_rotate_compare(data, image, offset, angle, plot=False, time=None):
	shiftimage = shift_image(image, offset)
	rotimage = rotate_image(shiftimage, angle)
	nvalidpix = len(np.ma.masked_invalid(rotimage).compressed())
	diff = (rotimage - data)**2 / errorsq
	chisq = np.nansum(diff)/nvalidpix
	print("offset (%d, %d) complete" % (xstep, ystep))
	del(shiftimage, rotimage, nvalidpix)
	gc.collect()
	if plot:
		plot(data, rotimage, diff, time)
	return chisq

def match(simfiles, potfiles, filenum, xraypeak, peak_threshold = 0.7, pixels_to_shift = 10, suffix=''):
	chisqs = np.empty(len(simfiles))
	
	file = simfiles[filenum]
	potfile = potfiles[filenum]

	sdata, sheader = open_fits(file,0)
	sx = calibrate(sdata.shape[1],sheader,axis=1)
	sy = calibrate(sdata.shape[1],sheader,axis=2)
	sx -= sx.mean()
	sy -= sy.mean()
	print("Snap %d read in" % filenum)

	#interpolate
	f = interp2d(sx, sy, sdata)
	binned_data = f(x,y)
	print("Data binned")

	#align x-ray peaks
	offset, rot_angle = offset_tilt(potfile)
	coords -= bcg_peak 
	coords *= -1 #how much to shift by
	print("Xray peaks aligned")
	
	#optimise rescaling within circle 
	l = np.nansum(data*binned_data/errorsq)/np.nansum(binned_data**2 / errorsq)
	binned_data *= l
	
	if filenum > filemin: 
		chisqs_around_peak = np.zeros([pixels_to_shift, pixels_to_shift])
		for xstep in range(pixels_to_shift):
			for ystep in range(pixels_to_shift):
				chisqs_around_peak[xstep,ystep] = shift_rotate_compare(data, binned_data,offset=(xstep, ystep, pixels_to_shift/2))
				print("xstep = %d, ystep = %d done!" % (xstep, ystep))

		minx, miny = np.argmin(chisqs_around_peak)[0]
		chisqs[filenum] = best_rotation(data, rerolled_data, (minx, miny, pixels_to_shift//2))
		print("done!")

	else:
		print("halos not merged yet. done!")
		chisqs[filenum] = shift_rotate_compare(data, binned_data, (0,0,0))

	np.save('chisq%s' % suffix, chisqs)

def compare():
	pool = Pool(10)
	print(len(simfiles))
	for filenum in range(len(simfiles)):
		pool.apply_async(match, args=(simfiles, potfiles, filenum, xraypeak))
