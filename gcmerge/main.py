################################################################
## Make FITS files, compare them to obs, and save comparisons ##
################################################################

import glob, os, yt 
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import cosmology
import astropy.units as u
from astropy.io import fits
from astropy.wcs import utils, WCS
from bcg_dist import find_peak
from make_profiles import make_islands, select_feature, make_profiles

from astropy.constants import G, c 

# dirs = glob.glob('*e1*/*/*/*/*/*/*/fitsfiles')
dirs = ['1.1e15_3e14/c1=4.1/a1=1/c2=5.2/a2=1.0/vrel=1852/b=250',
'9e14_2.4e14/c1=4.1/a1=2/c2=5.2/a2=1.0/vrel=1652/b=250',
'9e14_2.4e14/c1=4.1/a1=2/c2=5.2/a2=1.0/vrel=2200/b=250',
'6e14_2.1e14/c1=4.1/a1=1/c2=5.2/a2=1.0/vrel=2200/b=250',
'4e14_2.1e14/c1=4.1/a1=1/c2=5.2/a2=1.0/vrel=2200/b=250',
'2e14_2.1e14/c1=4.1/a1=1/c2=5.2/a2=1.0/vrel=2200/b=250']

lcdm=cosmology.Planck15
cosmo = cosmology.default_cosmology.get()
z = 0.2323
Mpc_rad = lcdm.angular_diameter_distance(z).value
kpc_deg = Mpc_rad*u.Mpc.to('kpc')/u.rad.to('degree')

def einstein_radius(cosmo, M, zl, zs):
	Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
	Dl = cosmo.angular_diameter_distance(zl)
	Ds = cosmo.angular_diameter_distance(zs)
	return ((Dls/(Dl*Ds)).to('1/kpc') * (G*M/c**2).to('kpc'))**0.5 * Dl.to('kpc')

#einstein_radius(cosmo, 1e15*u.Msun, 0.2323, 1) = 162 kpc
#einstein_radius(cosmo, 8e14*u.Msun, 0.2323, 1) = 145 kpc
#einstein_radius(cosmo, 4e14*u.Msun, 0.2323, 1) = 101 kpc

homedir = '/gpfs/loomis/project/fas/nagai/uc24/a2146_gamer/' #on Grace
# on wiluna: '/charra/uchadaya/GAMER/'

obsfile = homedir+'ff.img.e300_7000_bin3_all_obsids_box_excl_point_sources.fits'
errfile = homedir+'ff.img.e300_7000_bin3_all_obsids_box_0s_to_1s_no_bkg_subtr2_thresh.fits'
obshead = fits.getheader(obsfile)                           

coords = [ "239.0007 +66.37329", "239.058 +66.3482"]
bcg_coords = SkyCoord(coords, unit=(u.deg), obstime="J2000")

ptg_coords = SkyCoord([str(obshead['RA_PNT'])+" "+str(obshead['DEC_PNT'])], 
	unit=u.deg, obstime='J2000')


data = fits.getdata(obsfile)
x = utils.pixel_to_skycoord(np.arange(data.shape[0]), np.arange(data.shape[1]), WCS(fits.getheader(obsfile)))

bcg1_pix = bcg_coords[0].match_to_catalog_sky(x)[0]  
bcg2_pix = bcg_coords[1].match_to_catalog_sky(x)[0]  

#2 - BCG dist to select snapshots

def distance(dir, ax, label, color):
	potfiles = glob.glob(dir+'/potential/*.fits')
	potfiles.sort()
	
	x = []
	dists = []

	for file in potfiles:
		try: 
			dist=find_peak(file) #kpc
		except TypeError:
			dist = np.nan
		dists.append(dist)
		x.append(int(file.split('_')[-1].split('.fits')[0]))
	
	ax.plot(x, dists, label=label, c=color)
	return ax
	

#3 - shift, rotate, make profiles

def compare_arrays(potfile, xrayfile, tempfile, bcg1_pix=None, bcg2_pix=None, offset=0, minlen_kpc=500):
	islandlist, temp_img, sb_img, xraypeak = make_islands(potfile, xrayfile, tempfile, bcg1_pix, bcg2_pix)
	print ("islands made")

	dx = fits.getheader(potfiles[0])['CDELT1']
	fig, ax = plt.subplots()

	plt.imshow(temp_img, origin = 'lower', cmap = cm.afmhot, norm=colors.Normalize(1,20))

	fig1, ax1 = plt.subplots()
	fig2, ax2 = plt.subplots()

	for isle in islandlist:
		if isle.count() > minlen_kpc/dx:
			print( isle)
			make_profile(temp_img.T,sb_img.T, islandlist, islandlist.index(isle),label=str(islandlist.index(isle)),centre='arc', ax=ax, ax1=ax1, ax2=ax2)

			#but the profiles currently start at r=0. instead, need to focus on the peak


	print("profiles made")

	handles, labels = ax1.get_legend_handles_labels()
	ax1.set_xlim(0,150)
	loc = ax1.get_xticks()
	xlabels = loc*fits.getheader(tempfiles[filenum])['CDELT1']
	ax1.set_xticks(loc[::2])
	ax1.set_xticklabels(['%d' % x for x in xlabels[::2]])
	ax1.set_xlabel('R (kpc)')
	ax1.set_ylabel('kT (keV)')
	fig1.legend(handles, labels)
	fig1.savefig('temperature_profiles_%d.png' % (filenum+offset))

	ax2.set_xlim(0,150)
	loc = ax2.get_xticks()
	xlabels = loc*fits.getheader(tempfiles[filenum])['CDELT1']
	ax2.set_xticks(loc[::2])
	ax2.set_xticklabels(['%d' % x for x in xlabels[::2]])
	ax2.set_xlabel('R (kpc)')
	ax2.set_ylabel(r'Photon emissivity (cts cm$^{-3} s^{-1}$)')
	fig2.legend(handles, labels)
	ax2.set_yscale('log')
	fig2.savefig('xraysb_profiles_%d.png' % (filenum+offset))
	print("formatting complete")
	print(filenum+offset, " done")


def distances(ndirs, colormap):
	fig, ax = plt.subplots()
	for dir in ndirs:
		M = dir.split('_')[0]
		c1 = dir.split('c1=')[1].split('/')[0]
		a1 = dir.split('a1=')[1].split('/')[0]
		c2 = dir.split('c2=')[1].split('/')[0]
		a2 = dir.split('a2=')[1].split('/')[0]
		vrel = dir.split('vrel=')[1].split('/')[0]
		b = dir.split('b=')[1].split('/')[0]
		# label = '%s, c1=%s, a1=%s, c2=%s, a2=%s' % (M, c1, a1, c2, a2)
		label = '%s' % M
		if '2200' in vrel:
			label+=', vrel=2200'
		if '100' in b:
			label+=', b=100'
		ax = distance(dir, ax, label, colormap(float(ndirs.index(dir))/len(ndirs)))

	plt.hlines(470, 9,16,color='k',linestyle='dotted')
	plt.hlines(420, 9,16,color='k',linestyle='dotted')
	plt.xlim(9,16)
	plt.ylim(1e2,1e3)
	xtix = plt.xticks()[0]
	plt.xticks(xtix, ['%0.1f' % x for x in xtix/10.])
	plt.legend(loc=9,ncol=2, fontsize=12)