##########################################
## Read in observed and simulated files ##
##########################################

import numpy as np 
import glob
from astropy import units, cosmology
from astropy.io import fits

def calibrate(length, header, axis=1):
	refpix = header['CRPIX'+str(axis)]
	refval = header['CRVAL'+str(axis)]
	step = header['CDELT'+str(axis)] 
	values = np.empty(length)
	for ind in range(length):
		fitsind = ind+1
		values[ind] = refval+(fitsind-refpix)*step
	return values

def init(obsfile, errfile, simdir=None, peak_only=False):
	#read data
	data, header = open_fits(obsfile)
	errorsq, errhead = open_fits(errfile)
	x = calibrate(data.shape[1],header,axis=1)
	y = calibrate(data.shape[1],header,axis=2)
	lcdm = cosmology.Planck15
	z = 0.2323 #Abell 2146
	data[data<0] = 1e-8 #excessive background subtraction
	print("Observation read in")
	
	#center
	xraypeak = np.argwhere(data == np.nanmax(data))
	x -= x[xraypeak[0][0]]
	y -= y[xraypeak[0][1]]

	#degrees to kpc
	Mpc_rad = lcdm.angular_diameter_distance(z)
	kpc_deg = Mpc_rad.to('kpc')/units.radian.in_units('degree')
	x *= kpc_deg.value
	y *= kpc_deg.value
	
	#read_sim
	if simdir:
		simfiles = glob.glob(simdir+'/*fits')
		simfiles.sort()

		times = np.arange(len(simfiles))/10.
	if peak_only:
		return xraypeak
	elif simdir:
		return simfiles, times, data, errorsq, x, y, xraypeak
	else:
		return data, errorsq, x, y, xraypeak
