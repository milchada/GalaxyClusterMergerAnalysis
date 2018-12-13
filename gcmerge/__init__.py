import numpy as np 
import matplotlib.pylab as plt
import matplotlib, glob, gc 
from astropy import units, cosmology
from read_fits import open_fits, calibrate
from scipy.interpolate import interp2d
from scipy.ndimage.interpolation import rotate, shift
from skimage.feature import peak_local_max
from multiprocessing.pool import ThreadPool as Pool

name = 'gcmerge'
__all__ = ['compare', 'plot']

makeplot = False
select_halo = False
xraysb_obs = '../ff.img.e300_7000_bin3_all_obsids_box_excl_point_sources.fits'
xraysb_err = '../ff.img.e300_7000_bin3_all_obsids_box_0s_to_1s_no_bkg_subtr2_thresh.fits'

def init(obsfile, errfile, simdir):
	#read data
	data, header = open_fits(obsfile, 0)
	errorsq, errhead = open_fits(errfile, 0)
	x = calibrate(data.shape[1],header,axis=1)
	y = calibrate(data.shape[1],header,axis=2)
	lcdm = cosmology.default_cosmology.get()
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

	#select halo - currently manual
	if select_halo:
		minbrightness = 1.25e-5
		halomask = np.ma.masked_greater(data, minbrightness).mask
		data *= halomask
	
	#read_sim
	simfiles = glob.glob(simdir+'/*fits')#xrayprojz/ on wiluna
	simfiles.sort()

	times = np.arange(len(simfiles))/10.
	return simfiles, times, data, errorsq, x, y, xraypeak

# if __name__=="__main__":
# 	simfiles = "xrayprojz"
# 	simfiles, times, data, errorsq, x, y, xraypeak = init(xraysb_obs, xraysb_err, simfiles)

