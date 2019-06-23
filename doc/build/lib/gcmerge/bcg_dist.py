#bcg distance
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob, os
from astropy import cosmology
from astropy.io import fits

lcdm=cosmology.Planck15
z = 0.2323
Mpc_rad = lcdm.angular_diameter_distance(z).value
kpc_deg = Mpc_rad*u.Mpc.to('kpc')/u.rad.to('degree')


#CENTER THESE COORDS IN THE CHANDRA FOV
obsfile = '../ff.img.e300_7000_bin3_all_obsids_box_excl_point_sources.fits'
errfile = 'ff.img.e300_7000_bin3_all_obsids_box_0s_to_1s_no_bkg_subtr2_thresh.fits'
obshead = fits.getheader(obsfile)                           

coords = [str(obshead['RA_TARG'])+" "+str(obshead['DEC_TARG']), "239.0007 +66.37329"]
bcg_coords = SkyCoord(coords, unit=(u.deg), obstime="J2000")

ptg_coords = SkyCoord([str(obshead['RA_PNT'])+" "+str(obshead['DEC_PNT'])], 
	unit=u.deg, obstime='J2000')

def pixel(degree,axis='1'):
	x = np.array([(ind - obshead['CRPIX'+axis])*obshead['CDELT'+axis] + obshead['CRVAL'+axis] 
		for ind in range(obshead['NAXIS'+axis])])
	return np.argmin(abs(x-degree))

bcg1_pix = [pixel(bcg_coords[0].ra.value), pixel(bcg_coords[0].dec.value, '2')]
bcg2_pix = [pixel(bcg_coords[1].ra.value), pixel(bcg_coords[1].dec.value, '2')]

#manually correct for offset between bcgs and xray peak
# from read_db import init

# dist_deg = c[1].separation(c[0]).value #degree
# dist_kpc = dist_deg * kpc_deg

def angle(pt1, pt2):
	relative_pos = pt2 - pt1
	return np.angle(relative_pos[0] + 1j*relative_pos[1], deg=True)

#now find potential minima
from find_peaks import find_peak

files = glob.glob(os.getcwd()+'/fitsfiles/potential/zslice/*fits')
files.sort()
file = files[5]
peaks = find_peak(file) #x1, y1, x2, y2 in pix
if minima[0] < minima[1]:
	peak1 = peaks[:2]
	peak2 = peaks[2:]
else:
	peak1 = peaks[2:]
	peak2 = peaks[:2]

shift_pix = peak1 - bcg1_pix
simfiles, times, data, errorsq, x, y, xp = init(obsfile, errfile, fitsdir)
correction = bcg1_pix - xp
bcg1_pix -= correction
bcg2_pix -= correction

sx = calibrate(img.shape[1],simhead,axis=1)
sy = calibrate(img.shape[1],simhead,axis=2)
sx -= sx.mean()
sy -= sy.mean()
f = interp2d(sx, sy, img)
binned_data = f(x,y)
rolled_data = shift(binned_data, shift = coords[0])
rotate_angle = angle(bcg2_pix, bcg1_pix) - angle(peaks2, peaks1)
#how to map pix in sim to pix in obs
distance = np.linalg.norm(peaks[2:] - peaks[:2], axis=0) * fits.getheader(files[0])['CDELT1'] #kpc