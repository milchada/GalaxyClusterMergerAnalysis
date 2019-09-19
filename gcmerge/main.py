################################################################
## Make FITS files, compare them to obs, and save comparisons ##
################################################################

import glob, os, yt 
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import cosmology
import astropy.units as u
from astropy.io import fits
from read_sim import make_fits
from read import open_fits, calibrate
from bcg_dist import bcg_separation
from lensing_array import obs_res, shift_rotate_trim
from compare import shift_rotate_compare
from find_features import 
from fit_arcs import *
from make_profiles import *

lcdm=cosmology.Planck15
z = 0.2323
Mpc_rad = lcdm.angular_diameter_distance(z).value
kpc_deg = Mpc_rad*u.Mpc.to('kpc')/u.rad.to('degree')

homedir = '/gpfs/loomis/project/fas/nagai/uc24/a2146_gamer/' #on Grace
# on wiluna: '/charra/uchadaya/GAMER/'
#CENTER THESE COORDS IN THE CHANDRA FOV
obsfile = homedir+'ff.img.e300_7000_bin3_all_obsids_box_excl_point_sources.fits'
errfile = homedir+'ff.img.e300_7000_bin3_all_obsids_box_0s_to_1s_no_bkg_subtr2_thresh.fits'
obshead = fits.getheader(obsfile)                           

coords = [("239.058 +66.3482"), "239.0007 +66.37329"]
bcg_coords = SkyCoord(coords, unit=(u.deg), obstime="J2000")

ptg_coords = SkyCoord([str(obshead['RA_PNT'])+" "+str(obshead['DEC_PNT'])], 
	unit=u.deg, obstime='J2000')

from astropy.wcs import utils, WCS

data = fits.getdata(obsfile)
x = utils.pixel_to_skycoord(np.arange(data.shape[0]), np.arange(data.shape[1]), WCS(fits.getheader(obsfile)))

bcg1_pix = bcg_coords[0].match_to_catalog_sky(x)[0]  
bcg2_pix = bcg_coords[1].match_to_catalog_sky(x)[0]  

#2 - BCG dist to select snapshots
potfiles = glob.glob('fitsfiles/potential/*.fits')
potfiles.sort()
xrayfiles = glob.glob('fitsfiles/photon_emissivity/*.fits')  
xrayfiles.sort()
tempfiles = glob.glob('fitsfiles/temperature/*.fits')  
tempfiles.sort()

distmax = 470 #kpc
thetamax = np.deg2rad(20)
distmin = distmax*np.cos(thetamax)
shortlist = []
rot_angle = []
shifts = []
dists = []

for file in potfiles:
	dist=bcg_separation(file) #kpc
	dists.append(dist)
	if (dist > distmin) & (dist < distmax):
		shortlist.append(file.split('_')[-1].split('.')[0])
		
plt.clf()
plt.plot(np.arange(len(dists)), dists)
xtix = plt.xticks()[0]
plt.hlines(distmin, xtix.min(), xtix.max())
plt.hlines(distmax, xtix.min(), xtix.max())
plt.xticks(xtix, ['%0.1f' % x for x in xtix/10.])
plt.xlabel(r't$_{sim}$ (Gyr)')
plt.ylabel('d (kpc)')
plt.ylim(0,1500)
plt.savefig('distance.png')
#3 - shift, rotate, compare
#compare.py already does this

#4 - shift, rotate, make profiles
simfiles, times, data, errorsq, x, y, xp = init(obsfile, errfile, fitsdir)

def compare_arrays(filenum, bcg1_pix, bcg2_pix):
	islandlist, temp_img, sb_img, xraypeak = make_islands(filenum, bcg1_pix, bcg2_pix)
	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()
	
	for island in range(len(islandlist)):
		if islandlist[island].count() > 50:
			ax1, ax2 = make_profile(temp_img,sb_img, islandlist, island,label=str(island),centre='arc', ax1=ax1, ax2=ax2)
	handles, labels = ax1.get_legend_handles_labels()
	ax1.set_xlim(0,150)
	loc = ax1.get_xticks()
	xlabels = loc*fits.getheader(tempfiles[filenum])['CDELT1']
	ax1.set_xticks(loc[::2])
	ax1.set_xticklabels(['%d' % x for x in xlabels[::2]])
	ax1.set_xlabel('R (kpc)')
	ax1.set_ylabel('kT (keV)')
	fig1.legend(handles, labels)
	fig1.savefig('temperature_profiles_%d.png' % filenum)

	ax2.set_xlim(0,150)
	loc = ax2.get_xticks()
	xlabels = loc*fits.getheader(tempfiles[filenum])['CDELT1']
	ax2.set_xticks(loc[::2])
	ax2.set_xticklabels(['%d' % x for x in xlabels[::2]])
	ax2.set_xlabel('R (kpc)')
	ax2.set_ylabel(r'Photon emissivity (cts cm$^{-3} s^{-1}$)')
	fig2.legend(handles, labels)
	ax2.set_yscale('log')
	fig2.savefig('xraysb_profiles_%d.png' % filenum)