#########################################
####### may merge with align_BCGs #######
## currently: finding BCGs in obs, sim ##
#########################################

import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob, os
from astropy import cosmology
from astropy.io import fits
from sklearn.cluster import KMeans
from scipy.interpolate import interp2d

lcdm=cosmology.Planck15
z = 0.2323
Mpc_rad = lcdm.angular_diameter_distance(z).value
kpc_deg = Mpc_rad*u.Mpc.to('kpc')/u.rad.to('degree')

homedir = '/charra/uchadaya/GAMER/'
#CENTER THESE COORDS IN THE CHANDRA FOV
obsfile = homedir+'ff.img.e300_7000_bin3_all_obsids_box_excl_point_sources.fits'
errfile = homedir+'ff.img.e300_7000_bin3_all_obsids_box_0s_to_1s_no_bkg_subtr2_thresh.fits'
obshead = fits.getheader(obsfile)                           

coords = [("239.058 +66.3482"), "239.0007 +66.37329"]
bcg_coords = SkyCoord(coords, unit=(u.deg), obstime="J2000")

ptg_coords = SkyCoord([str(obshead['RA_PNT'])+" "+str(obshead['DEC_PNT'])], 
	unit=u.deg, obstime='J2000')

def pixel(degree,axis='1'):
	x = np.array([(ind - obshead['CRPIX'+axis])*obshead['CDELT'+axis] + obshead['CRVAL'+axis] 
		for ind in range(obshead['NAXIS'+axis])])
	return np.argmin(abs(x-degree))

# bcg1_pix = [pixel(bcg_coords[0].ra.value), pixel(bcg_coords[0].dec.value, '2')]
# bcg2_pix = [pixel(bcg_coords[1].ra.value), pixel(bcg_coords[1].dec.value, '2')]

from find_peaks import find_peak

def bcg_separation(file, distmax_kpc=470, ret_peaks=False, xmin=300, xmax=1200, ymin = 300, ymax = 1200, axis=None):	
	data = fits.getdata(file)
	header = fits.getheader(file)
	md = 3
	maxsep = distmax_kpc/header['CDELT1']
	oldlen = len(data)**2 
	while md < maxsep:
		peaks = find_peak(file, axis=axis, min_distance = md, xmin=xmin, xmax=xmax, ymin = ymin, ymax = ymax) #x1, y1, x2, y2 in pix  
		md += 1 
		if len(peaks) == 2: 
			break 
		elif len(peaks) == oldlen: 
			break 
		oldlen = len(peaks) 
	
	if len(peaks) == 1:
		print( "cluster relaxed")
	else:
		if len(peaks) > 2: 
			kmeans = KMeans(n_clusters=2)
			kmeans.fit(peaks)
			cid = kmeans.predict(peaks)
			peaks = (kmeans.cluster_centers_).astype(int)
		
		if ret_peaks:
			return peaks
		else:
			return sep