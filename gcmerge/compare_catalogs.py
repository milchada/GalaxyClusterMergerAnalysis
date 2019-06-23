##########################################
# Compare lensing observations with sims #
##########################################

import numpy as np
import pandas as pd
import os, glob
from astropy import units, cosmology, convolution
from lensing_array import dataFrame
from find_peaks import find_peak

lensing_cat = 'A2146_WL_HSTcatalog/a2146.moderate_2Dbin_px291.dat'
bcga_pos = (2712.764, 6332.1972)
bcgb_pos = (-1839, -4146)#relative to BCG-A

lensing_data = dataFrame(lensing_cat)
obs_res = lensing_data['x'][1] - lensing_data['x'][0] #in Mpc

rhoprojfile = 'rhoproj_z_14.fits'
potential_file = 'potential_z_14.fits'

def angle(pt1, pt2):
	relative_pos = pt2 - pt1
	return np.angle(relative_pos[0] + 1j*relative_pos[1], deg=True)

def rotate(simcat, angle):
	theta = np.deg2rad(angle)
	rot = np.array([[np.cos(theta), -np.sin(theta)],
					[np.sin(theta),np.cos(theta)]])
	xy = np.stack(simcat['x'],simcat['y'],axis=1)
	return np.dot(xy,rot)

def shift_rotate_trim(rhoprojfile, potential_file, obscat):
	simcat = pd.DataFrame(data=np.genfromtxt(simfile), columns = ['x','y','k','g1','g2'])
		#kappa is essentially the projected density, so just use its maximum as peak
	peaks = find_peaks(potential_file)
	scale = fits.getheader(potential_file)['CDELT1'] #kpc/pix
	peaks *= scale/1000. #Mpc
	simcat['x'] -= peaks[0][0]
	simcat['y'] -= peaks[0][1]

	theta = angle(peaks[1], peaks[0])
	xy = rotate(simcat, theta)
	simcat['x'] = xy[:,0]
	simcat['y'] = xy[:,1]

	xobs = obscat['x']
	yobs = obscat['y']

	simsub = simcat[(simcat[:,0] > xobs.min()) & (simcat[:,0] < xobs.max()) 
			& (simcat[:,1] > yobs.min()) & (simcat[:,1] < yobs.max())]
	# np.save('%s_contracted' % simfile.split('_cat')[0], simsub)
	return simsub

def convolve(array, sigma= 50/14.):
	kernel = convolution.Gaussian2DKernel(stddev=sigma)
	return convolution.convolve(array, kernel)

def match(simfiles, filenum, peak_threshold = 0.8, pixels_to_shift = 10, suffix='_lens'):
	file = simfiles[filenum]
	simcat = shift_rotate_trim(file, lensing_data)	

	#bin simulated catalog to match obs
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
	
	#everything is in Mpc
	print("Snap %d read in" % filenum)

	chisq_g1 = ((binned_cat['g1'] - lensing_data['g1'])**2)/lensing_data['sigma']**2
	chisq_g2 = ((binned_cat['g2'] - lensing_data['g2'])**2)/lensing_data['sigma']**2
	#I do want to check that the simulated catalog produces reduced shear, which is what the catalog has *I think*
	
	print("files saved")
