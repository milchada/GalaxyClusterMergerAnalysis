##########################################
## Read in observed and simulated files ##
##########################################

import numpy as np 
import matplotlib.pylab as plt
import matplotlib, glob, gc 
from astropy import units, cosmology
from astropy.io import fits

def open_fits(filename, tablename):
	f = fits.open(filename)
	data = f[tablename].data
	header = f[tablename].header
	return data, header

def calibrate(length, header, axis=1):
	refpix = header['CRPIX'+str(axis)]
	refval = header['CRVAL'+str(axis)]
	step = header['CDELT'+str(axis)] 
	values = np.empty(length)
	for ind in range(length):
		fitsind = ind+1
		values[ind] = refval+(fitsind-refpix)*step
	return values

def init(obsfile, errfile, simdir, peak_only=False):
	#read data
	data, header = open_fits(obsfile, 0)
	errorsq, errhead = open_fits(errfile, 0)
	x = calibrate(data.shape[1],header,axis=1)
	y = calibrate(data.shape[1],header,axis=2)
	lcdm = cosmology.Planck15
	z = 0.2323 #Abell 2146

	#center
	Mpc_rad = lcdm.angular_diameter_distance(z)
	kpc_deg = Mpc_rad.to('kpc')/units.radian.in_units('degree')
	x *= kpc_deg.value
	y *= kpc_deg.value
	x -= x.mean()
	y -= y.mean()
	
	#degrees to kpc
	range = data.max()/data[data>0].min()
	xraypeak = np.argwhere(data == np.nanmax(data))
	data[data<0] = 0 #excessive background subtraction
	print("Observation read in")

	#read_sim
	simfiles = glob.glob(simdir+'/*fits')
	simfiles.sort()

	times = np.arange(len(simfiles))/10.
	if peak_only:
		return xraypeak
	else:
		return simfiles, times, data, errorsq, x, y, xraypeak
