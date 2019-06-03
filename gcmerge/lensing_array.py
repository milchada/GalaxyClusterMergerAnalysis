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

# peak = np.argmax(lensing_data['g1']**2 + lensing_data['g2']**2)
bcga_pos = (2712.764, 6332.1972)  #in pix

lensing_data['x'] -= bcg_pos[0]#lensing_data['x'][peak]
lensing_data['y'] -= bcg_pos[1]#lensing_data['y'][peak]
lensing_data['x'] *= arcsecPerPix*mpcPerArcsec
lensing_data['y'] *= arcsecPerPix*mpcPerArcsec

simfiles = glob.glob('rhoproj/*fits')
simfiles.sort()

def trim(simfile, obscat):
	simcat = pd.DataFrame(data=np.genfromtxt(simfile), columns = ['x','y','k','g1','g2'])
		#kappa is essentially the projected density, so just use its maximum as peak
	peak = np.argmax(simcat['k'])
	peak_x = simcat.iloc[peak]['x']
	peak_y = simcat.iloc[peak]['y']
	simcat['x'] -= peak_x
	simcat['y'] -= peak_y

	xobs = obscat['x']
	yobs = obscat['y']

	simsub = simcat[(simcat[:,0] > xobs.min()) & (simcat[:,0] < xobs.max()) 
			& (simcat[:,1] > yobs.min()) & (simcat[:,1] < yobs.max())]
	# np.save('%s_contracted' % simfile.split('_cat')[0], simsub)
	return simsub

best_angles = np.load('best_angle.npy')
best_shifts = np.load('best_offset.npy')
chisqs = np.empty([len(simfiles), 3])

def shift_and_rotate(simcat, offset, angle):
	xstep, ystep, startshift = offset
	simcat['x'] += xstep-startshift
	simcat['y'] += ystep-startshift
	theta = np.deg2rad(angle)
	rot = np.array([[np.cos(theta), -np.sin(theta)],
					[np.sin(theta),np.cos(theta)]])
	xy = np.stack(simcat['x'],simcat['y'],axis=1)
	return np.dot(xy,rot)

def match(simfiles, filenum, peak_threshold = 0.8, pixels_to_shift = 10, suffix='_lens'):
	file = simfiles[filenum]
	# if not os.path.exists(file.split('.fits')[0]+'_cat.dat'):
	# 	cc.processInput(file) #i.e. convert to catalog 
	
	simcat = trim(file, lensing_data)

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

	# xbest, ybest = best_shifts[filenum] #no - this is in x-ray image pixels, ≠ lensing image pixels argh
	xbest, ybest = (0,0) #since i've lined up peaks
	rotxy = shift_and_rotate(simcat, lensing_data, (xbest, ybest, pixels_to_shift/2), best_angles[filenum])
	simcat['x'] = rotxy[:,0]
	simcat['y'] = rotxy[:,1]

	def chisq_for_cat_lists(simcat, lensing_data):
		"Is there a better way than turning them both into 2D arrays and interpolating?"
		return chisq
	
	np.save('chisq_lens.npy' % suffix, chisqs)
	print("files saved")

if __name__=="__main__":
	for filenum in range(len(simfiles)):
		match(simfiles, filenum)
