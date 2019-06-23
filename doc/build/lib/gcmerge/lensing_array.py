from simulated_WLfield import convertToCatalog as cc
import numpy as np
import pandas as pd
import matplotlib.pylab as plt 
from astropy import units, cosmology
from scipy.interpolate import interp2d
import os, glob

lensing_cat = 'a2146.moderate_2Dbin_px291.dat'
lensing_data = pd.DataFrame(data=np.genfromtxt(lensing_cat), columns = ['x','y','g1','g2','wi','sigma'])
arcsecPerPix = 202./4096 #HST-WFC param
zobs = 0.2323
mpcPerArcsec = cosmology.Planck15.angular_diameter_distance(zobs)*units.arcsec.in_units('radian')

peak = np.argmax(lensing_data['g1']**2 + lensing_data['g2']**2)
bcga_pos = (2712.764, 6332.1972)  #in pix

lensing_data['x'] -= lensing_data['x'][peak]
lensing_data['y'] -= lensing_data['y'][peak]
lensing_data['x'] *= arcsecPerPix*mpcPerArcsec
lensing_data['y'] *= arcsecPerPix*mpcPerArcsec

simfiles = glob.glob('rhoproj/*fits')
simfiles.sort()

def cat_to_array(lensing_cat):
	xs = lensing_cat['x'].unique()
	ys = lensing_cat['y'].unique()
	g1_array = np.empty([len(xs), len(ys)])
	g2_array = np.empty([len(xs), len(ys)])
	for row in range(len(lensing_data)):
		xind = np.argwhere(xs == lensing_cat['x'][row])
		yind = np.argwhere(ys == lensing_cat['y'][row])
		print(xind, yind)
		g1_array[xind, yind] = lensing_cat['g1'][row]
		g2_array[xind, yind] = lensing_cat['g2'][row]
	return g1_array, g2_array

g1_array, g2_array = cat_to_array(lensing_data)

best_angles = np.load('best_angle.npy')
best_shifts = np.load('best_offset.npy')
chisqs = np.empty([len(simfiles), 3])

def rotate_image(data, image, angle):
	image[np.isnan(image)] = 0 #otherwise rotation gets fucked
	rotimage = rotate(image, angle = angle)
	datasize = data.shape[0]
		#no we have to start from the center
	rotsize = rotimage.shape[0]
	if rotsize != datasize:
		xstart = (rotsize - datasize)//2
		ystart = (rotsize-datasize)//2
		rotimage = rotimage[xstart:xstart+datasize, ystart:ystart+datasize]
		rotimage[rotimage == 0] = np.nan
	diff = ((rotimage - data))**2/errorsq
	diff[diff==np.inf]=np.nan
	return rotimage, diff

def shift_and_rotate(data, image, offset, angle):
	xstep, ystep, startshift = offset
	imagecut = np.roll(np.roll(image, xstep-startshift, axis=0), ystep-startshift, axis=1)
	imagecut[np.isnan(imagecut)] = 0 #otherwise rotation gets fucked
	rotimagecut, diff = rotate_image(data, image, angle)
	nvalidpix = len(np.ma.masked_invalid(diff).compressed())
	return np.nansum(diff)/nvalidpix

def match(simfiles, filenum, peak_threshold = 0.8, pixels_to_shift = 10, suffix='_lens'):
	file = simfiles[filenum]
	# if not os.path.exists(file.split('.fits')[0]+'_cat.dat'):
	# 	cc.processInput(file) #i.e. convert to catalog 
	
	simcat = pd.DataFrame(data=np.genfromtxt(file.split('.fits')[0]+'_cat.dat'), 
		columns = ['x','y','k','g1','g2'])
	simpeak = np.argmax(simcat['g1']**2 + simcat['g2']**2)
	simcat['x'] -= simcat['x'][simpeak]
	simcat['y'] -= simcat['y'][simpeak] 

	xmask = (simcat['x'] > lensing_data['x'].min()) * (simcat['y'] < lensing_data['x'].max())
	ymask = (simcat['x'] > lensing_data['y'].min()) * (simcat['y'] < lensing_data['y'].max())
	inds = np.arange(len(simcat))*xmask*ymask
	simcat = simcat.iloc[inds[inds>0]]
	simcat.index = np.arange(len(simcat))
	
	simcat['xbins'] = np.digitize(simcat['x'].unique(), lensing_data['x'].unique())
	simcat['ybins'] = np.digitize(simcat['y'].unique(), lensing_data['y'].unique())
	binned_cat = pd.DataFrame()
	binned_cat['x'] = lensing_data['x']
	binned_cat['y'] = lensing_data['y']
	for row in range(len(binned_cat)):
		xbin = np.argwhere(lensing_data['x'].unique == binned_cat['x'][row])
		ybin = np.argwhere(lensing_data['y'].unique == binned_cat['y'][row])
		binned_cat['g1'] = np.mean(simcat['g1'][simcat['xbins']==xbin][simcat['ybins'] == ybin])
		binned_cat['g2'] = np.mean(simcat['g2'][simcat['xbins']==xbin][simcat['ybins'] == ybin])	
	
	g1_sim, g2_sim = cat_to_array(binned_cat)
	
	#everything is in Mpc
	print("Snap %d read in" % filenum)

	# xbest, ybest = best_shifts[filenum] #no - this is in x-ray image pixels, ≠ lensing image pixels argh
	xbest, ybest = (0,0) #since i've lined up peaks
	chisqs[filenum][0] = shift_and_rotate(g1_array, g1_sim, (xbest, ybest, pixels_to_shift/2), best_angles[filenum])
	chisqs[filenum][1] = shift_and_rotate(g2_array, g2_sim, (xbest, ybest, pixels_to_shift/2), best_angles[filenum])
	chisqs[filenum][1] = shift_and_rotate(g1_array**2+g2_array**2, g1_sim**2,g2_sim**2, (xbest, ybest, pixels_to_shift/2), best_angles[filenum])

	np.save('chisq_lens.npy' % suffix, chisqs)
	print("files saved")

if __name__=="__main__":
	for filenum in range(len(simfiles)):
		match(simfiles, filenum)
