import numpy as np
import gc, glob
from image_matching import *
from scipy.ndimage.measurements import center_of_mass as com 


obsfile = '../kT_out.fits'
errfile = '../err_kT_out.fits'
simdir = 'tempproj'
simfiles, times, data, errorsq, x, y, xp = init(obsfile, errfile, simdir)
#xray peak should come from surface brightness, coz higher res
s1, t1, d1, er1, x1, y1, xraypeak = init(xraysb_obs, xraysb_err, 'rhoproj')
del(s1, t1, d1, er1, x1, y1, xp)
gc.collect()

def rotate_image(data, image, angle):
	image[np.isnan(image)] = 0 #otherwise rotation gets fucked
	rotimage = rotate(image, angle = angle)
	datasize = data.shape[0]
		#no we have to start from the center
	rotsize = rotimage.shape[0]
	if rotsize != datasize:
		xstart = (rotsize - datasize)/2
		ystart = (rotsize-datasize)/2
		rotimage = rotimage[xstart:xstart+datasize, ystart:ystart+datasize]
		rotimage[rotimage == 0] = np.nan
	diff = ((rotimage - data))**2/errorsq
	diff[diff==np.inf]=np.nan
	return rotimage, diff

angles = np.linspace(0, 350, 36)
def best_rotation(data, image, offset, makeplot=False):
	xstep, ystep, startshift = offset
	chisq = []
	for angle in angles:
		imagecut = np.roll(np.roll(image, xstep-startshift, axis=0), ystep-startshift, axis=1)
		imagecut[np.isnan(imagecut)] = 0 #otherwise rotation gets fucked
		rotimagecut, diff = rotate_image(data, image, angle)
		nvalidpix = len(np.ma.masked_invalid(diff).compressed())
		chisq.append(np.nansum(diff)/nvalidpix)
		# print "angle %d complete" % angle
	print "offset (%d, %d) complete" % (xstep, ystep)
	best_angle = angles[chisq.index(min(chisq))]
	if makeplot:
		plot(image, times[filenum], best_angle) #so the angle info is stored in the image name
	del(rotimagecut, imagecut)
	gc.collect()
	return min(chisq), best_angle

def match(simfiles, filenum, xraypeak, peak_threshold = 0.1, pixels_to_shift = 10, makeplot=False, suffix=''):
	if glob.glob('chisq%s.npy' % suffix):
		chisqs = np.load('chisq%s.npy' % suffix)
		best_angles = np.load('best_angle%s.npy' % suffix)
		best_shifts = np.load('best_offset%s.npy' % suffix)
	else:
		best_shifts = np.empty([len(simfiles), 2])
		best_angles = np.empty(len(simfiles))
		chisqs = np.empty(len(simfiles))
	try:
		file = simfiles[filenum]
		sdata, sheader = open_fits(file,0)
		sx = calibrate(sdata.shape[1],sheader,axis=1)
		sy = calibrate(sdata.shape[1],sheader,axis=2)
		sx -= sx.mean()
		sy -= sy.mean()
		print "Snap %d read in" % filenum

		#interpolate
		f = interp2d(sx, sy, sdata)
		binned_data = f(x,y)
		print "Data binned"

		#align x-ray peaks
		coords = peak_local_max(binned_data,threshold_rel=peak_threshold)
		print len(coords), " peaks"
		coords -= xraypeak 
		coords *= -1 #how much to shift by
		rolled_data = shift(binned_data, shift = coords[0])
		print "Xray peaks aligned"
		if select_halo:
			rolled_data *= halomask
			min_brightness = np.nanmax(rolled_data)/range
			rolled_data[rolled_data<min_brightness] = min_brightness
			rolled_data[rolled_data == 0] = np.nan
		#optimise rescaling within circle 
		l = np.nansum(data*rolled_data/errorsq)/np.nansum(rolled_data**2 / errorsq)
		rolled_data *= l
		
		if len(coords) == 1:
			print "System relaxed at snap %d. done!" % filenum
			best_shifts[filenum] = coords
			chisqs[filenum], best_angles[filenum] = best_rotation(data, rolled_data, (0,0,0), makeplot=makeplot)

		elif len(coords) == 2: 
			if filenum > 12: #temp fix, by eye
				chisqs_around_peak1 = np.zeros([pixels_to_shift,pixels_to_shift])
				chisqs_around_peak2 = np.zeros([pixels_to_shift,pixels_to_shift])
				best_rot = 0 #really just a placeholder
				rerolled_data = shift(binned_data, shift = coords[1])
				if select_halo:
					rerolled_data *= halomask
					min_brightness = np.nanmax(rolled_data)/range
					rerolled_data[rerolled_data<min_brightness] = min_brightness
					rerolled_data[rerolled_data == 0] = np.nan
				l = np.nansum(data*rerolled_data/errorsq)/np.nansum(rerolled_data**2 / errorsq)
				rerolled_data *= l

				"""cythonize/vectorise this nested loop"""

				for xstep in xrange(pixels_to_shift):
					for ystep in xrange(pixels_to_shift):
						chisqs_around_peak1[xstep,ystep], best_rot =best_rotation(data, rolled_data,offset=(xstep, ystep, pixels_to_shift/2))
						chisqs_around_peak2[xstep,ystep], best_rot =best_rotation(data, rerolled_data,offset=(xstep, ystep, pixels_to_shift/2))

				if chisqs_around_peak2.min() < chisqs_around_peak1.min():
					minx, miny = np.argwhere(chisqs_around_peak2 == chisqs_around_peak2.min())[0]
					best_shifts[filenum] = coords[1] - np.array([pixels_to_shift/2,pixels_to_shift/2]) + np.array([minx, miny])
					chisqs[filenum], best_angles[filenum] = best_rotation(data, rerolled_data, (minx, miny, pixels_to_shift/2), makeplot = makeplot)
				else:
					minx, miny = np.argwhere(chisqs_around_peak1 == chisqs_around_peak1.min())[0]
					best_shifts[filenum] = coords[0] - np.array([pixels_to_shift/2,pixels_to_shift/2]) + np.array([minx, miny])
					chisqs[filenum], best_angles[filenum] = best_rotation(data, rolled_data, (minx, miny, pixels_to_shift/2), makeplot = makeplot)
				print "done!"
			else:
				print "halos not merged yet. done!"
				chisqs[filenum], best_angles[filenum] = best_rotation(data, rolled_data, (0,0,0), makeplot=makeplot)
				best_shifts[filenum] = coords[0]

		else: #> 2 local peaks
			centroid = np.mean(coords,axis = 0) #remember, the xray peak has already been subtracted from this
			rolled_data = shift(binned_data, [int(round(centroid[0])), int(round(centroid[1]))])
			if select_halo:
				rolled_data *= halomask
				min_brightness = np.nanmax(rolled_data)/range
				rolled_data[rolled_data<min_brightness] = min_brightness
				rolled_data[rolled_data == 0] = np.nan

			l = np.nansum(data*rolled_data/errorsq)/np.nansum(rolled_data**2 / errorsq)
			rolled_data *= l
			chisqs_around_centroid = np.empty([pixels_to_shift*2,pixels_to_shift*2])
			
			"""cythonize/vectorise this nested loop"""

			for xstep in xrange(pixels_to_shift*2):
				for ystep in xrange(pixels_to_shift*2):
					chisqs_around_centroid[xstep,ystep], best_rot =best_rotation(data, rolled_data,offset=(xstep, ystep, pixels_to_shift))

			minx, miny = np.argwhere(chisqs_around_centroid == chisqs_around_centroid.min())[0]
			best_shifts[filenum] = centroid - np.array([pixels_to_shift,pixels_to_shift]) + np.array([minx,miny])
			chisqs[filenum], best_angles[filenum] = best_rotation(data, rolled_data, (minx, miny, pixels_to_shift), makeplot = makeplot)

		np.save('chisq%s' % suffix, chisqs)
		np.save('best_offset%s' % suffix, best_shifts)	
		np.save('best_angle%s' % suffix, best_angles)
		
	except (IOError, ValueError, IndexError):
		print "Error in snap # ", filenum

cs = np.load('chisq_temp.npy')
filestogo = np.argwhere(cs<1)
filenums = [filestogo[i][0] for i in xrange(len(filestogo))]
pool = Pool(10)
for filenum in filenums:#np.argwhere(cs==0):
	pool.apply_async(match, args=(simfiles, filenum, xraypeak), kwds=dict(suffix='_temp'))


def comparelens():
	obsfile = 'A2146_WL_HSTcatalog/a2146.moderate_2Dbin_px170.dat'
	errfile = '?'
	