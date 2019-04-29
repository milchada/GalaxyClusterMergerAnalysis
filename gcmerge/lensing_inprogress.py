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

angles = np.linspace(0, 350, 36)
def best_rotation(data, image, offset):
	xstep, ystep, startshift = offset
	dx = image['x'][1] - image['x'][0]
	coords = np.empty([len(image), 2])
	coords[:,0] = image['x'] - dx*(startshift+xstep)
	coords[:,1] = image['y'] - dx*(startshift+ystep)
	chisq = []
	for angle in angles:
		theta = np.radians(angle)
		c, s = np.cos(theta), np.sin(theta)
		rot = np.array([(c, -s),(s, c)])
		xrot, yrot = coords.dot(rot)[:,0], coords.dot(rot)[:,1]
		print('rotated by %s degrees' % angle)
		xmask = (xrot > lensing_data['x'].min()) * (xrot < lensing_data['x'].max())
		ymask = (yrot > lensing_data['y'].min()) * (yrot < lensing_data['y'].max())
		interp1 = interp2d(xrot[xmask*ymask][::10], yrot[xmask*ymask][::10], image['g1'][xmask*ymask][::10])
		interp2 = interp2d(xrot[xmask*ymask][::10], yrot[xmask*ymask][::10], image['g2'][xmask*ymask][::10])
		print('interpolation complete') #this is the expensive step, hence mask and every 10th
		fit_g1 = interp1(lensing_data['x'], lensing_data['y'])
		fit_g2 = interp2(lensing_data['x'], lensing_data['y'])
		cs = ((lensing_data['g1'].values - fit_g1)**2 + (lensing_data['g2'].values - fit_g2)**2) / (lensing_data['sigma'].values**2)
		print('chisq computed')
		cs[cs == np.inf] = np.nan
		chisq.append(np.nansum(cs)/(len(cs[cs>0])))
		# print("angle %d complete" % angle
	print("offset (%d, %d) complete" % (xstep, ystep))
	return min(chisq), angles[np.argmin(np.array(chisqs))]

def match(simfiles, filenum, peak_threshold = 0.8, pixels_to_shift = 10, suffix='_lens'):
	if glob.glob('chisq%s.npy' % suffix):
		chisqs = np.load('chisq%s.npy' % suffix)
		best_angles = np.load('best_angle%s.npy' % suffix)
		best_shifts = np.load('best_offset%s.npy' % suffix)
	else:
		best_shifts = np.empty([len(simfiles), 2])
		best_angles = np.empty(len(simfiles))
		chisqs = np.empty(len(simfiles))
	
	file = simfiles[filenum]
	if not os.path.exists(file.split('.fits')[0]+'_cat.dat'):
		cc.processInput(file) #i.e. convert to catalog 
	
	simcat = pd.DataFrame(data=np.genfromtxt(file.split('.fits')[0]+'_cat.dat'), 
		columns = ['x','y','k','g1','g2'])
	simpeak = np.argmax(simcat['g1']**2 + simcat['g2']**2)
	simcat['x'] -= simcat['x'][peak]
	simcat['y'] -= simcat['y'][peak] 
	#everything is in Mpc
	print("Snap %d read in" % filenum)

	for xstep in range(pixels_to_shift):
		for ystep in range(pixels_to_shift):
			chisqs_around_peak[xstep,ystep], best_rot = best_rotation(lensing_data, simcat,offset=(xstep, ystep, pixels_to_shift/2))

	xbest, ybest = np.argwhere(chisqs_around_peak == chisqs_around_peak.min())[0]
	best_shifts[filenum] = xbest, ybest
	chisqs[filenum], best_angles[filenum] = best_rotation(lensing_data, simcat, (xbest, ybest, pixels_to_shift/2))
	np.save('chisq%s.npy' % suffix, chisqs)
	np.save('best_angle%s.npy' % suffix, best_angles)
	np.save('best_offset%s.npy' % suffix, best_shifts)
	print("files saved")

if __name__=="__main__":
	for filenum in range(len(simfiles)):
		match(simfiles, filenum)

#for some reason multiprocessing isn't working
#ok take a break, come back and write A2146 outline. how is it already 3pm???