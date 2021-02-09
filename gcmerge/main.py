################################################################
## Make FITS files, compare them to obs, and save comparisons ##
################################################################

import glob
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import cosmology
import astropy.units as u
from astropy.io import fits
from astropy.wcs import utils, WCS
from astropy.constants import G, c 
from bcg_dist import find_peak
from make_islands import make_islands
from measure_feature import select_feature, standoff_distance, upstream_shock_bcg


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
	
	dists = np.zeros((len(potfiles),2))

	for i in range(len(potfiles)):
		try: 
			dist=find_peak(file) #kpc
		except TypeError:
			dist = np.nan
		dists[i,0] = i
		dists[i,1] = dist
	dists = dists[dists[:,0] > 0]
	np.save('distance', dists)

#3 - shift, rotate, make profiles

def summary_statistics(potfile, xrayfile, tempfile, edgecontrast=None, minlen_kpc=500):
	#selections by contrast 
	#also default xrange yrange select central 400pix ~ 2.4Mpc of image
	islandlist, temp_img, xraypeak = make_islands(xrayfile=xrayfile, tempfile=tempfile, edgecontrast=edgecontrast)
	print ("islands made")

	dx = fits.getheader(tempfile)['CDELT1']
	fig, ax = plt.subplots()

	plt.imshow(temp_img, origin = 'lower', cmap = cm.afmhot, norm=colors.Normalize(1,20))

	for i in range(len(islandlist)):
		features = {}
		#selection by length
		if islandlist[i].count() > minlen_kpc/dx:
			print(isle)
			f, r, c = select_feature(temp_img, sb_img, islandlist, isle)
			plt.scatter(f[:,1], f[:,0], label=i)
			features[i] = f 
	plt.legend()

	if len(features.keys()) > 3:
		p1, p2 = find_peak(potfile, ret_peaks=False) #BCGs of clusters 1 and 2
		
		#so p2 is the bcg we're referring to
		#now what's "left" vs "right" of bcg? remember coordinates are in (y, x)
		left = []
		dist_bcg = np.zeros((2,4))
		i = 0

		if len(features.keys() > 4):  #select four longest features
			length = np.zeros(len(features.keys()))
			for key in features.keys():
				feature = feature[key]
				if feature[:,1].mean() < p2[1]: #left of image
					left.append(key)
				length[key] = len(feature)			
			
			longest = [np.argsort(length)[::-1]][:4]

			for key in longest:
				dist = np.linalg.norm(np.mean(features[key], axis=0) - p2)
				dist_bcg[i] = [key, dist]
				i += 1
		else: #only four features to start with
			for i in range(4):
				dist = np.linalg.norm(np.mean(features[key], axis=0) - p2)
				dist_bcg[i] = [i, dist]

		dist_bcg = dist_bcg[np.argsort(dist_bcg[:,1])] #nearest to farthest

		cf_ind, bs_ind = dist_bcg[:,0][dist_bcg[:,0] in left] #sorted by distance, features on the left
		spur_ind, us_ind = dist_bcg[:,0][dist_bcg[:,0] not in left] #sorted by distance, features on the right

		#ok i think this should work
		standoff 	 = standoff_distance(features[cf_ind], features[bs_index], dx = dx)
		upstream_bcg = upstream_shock_bcg(features[us_ind], bcg_xy=p2, dx=dx)
		
		return fig, ax, cf_ind, bs_ind, us_ind, p2, standoff, upstream_bcg
	else:
		print("Too few features detected at this snapshot")

#ok take a break then come back and test thisn

	