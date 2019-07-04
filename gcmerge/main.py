################################################################
## Make FITS files, compare them to obs, and save comparisons ##
################################################################

import glob, os, yt 
import numpy as np
from read_sim import sciencedir, simfiles, make_fits
from read import open_fits, calibrate
from bcg_dist import bcg_separation
from lensing_array import obs_res, shift_rotate_trim
from compare import shift_rotate_compare
from find_features import 
from fit_arcs import *

potdir = sciencedir+'/fitsfiles/potential/zslice/'

#1 - BCG dist to select snapshots
for filenum in range(len(simfiles)):
	if not os.path.isfile(potdir+'gravitational_potential_z_%s.fits' % filenum):
		make_fits(filenum, 'gravitational_potential', outputdir = potdir, slice=True)

potfiles = glob.glob(outputdir+'*fits')
potfiles.sort()

distmax = 470 #kpc
thetamax = np.deg2rad(20)
distmin = distmin*np.cos(thetamax)

shortlist = []
rot_angle = []
shifts = []
for file in potfiles:
	dist, angle, shift = bcg_separation(file).value #kpc
	if (dist > distmin) & (dist < distmax):
		shortlist.append(file.split('_')[-1].split('.')[0])
		rot_angle.append(angle)
		shifts.append(shift)
#2 - make X-ray emissivity and temperature files
emdir = sciencedir+'/fitsfiles/xray_sb/'

for filenum in shortlist:
	if not os.path.isfile(emdir+'xray_emissivity_0.3_7.0_keV_%s.fits' % filenum):
		make_fits(filenum, 'xray_emissivity_0.3_7.0_keV', outputdir=emdir, proj=True)

tempdir = sciencedir+'/fitsfiles/temp/'
for filenum in shortlist:
	if not os.path.isfile(emdir+'temperature_%s.fits' % filenum):
		make_fits(filenum, 'temperature', outputdir=tempdir, proj=True, weight_field='mazzotta_weighting')

#3 - make lensing files
boxsize = yt.load(simfiles[0]).domain_right_edge.to('Mpc')[0].value
image_res = float(boxsize)/obs_res
rhoprojdir = sciencedir+'/fitsfiles/rhoproj/'
for filenum in shortlist:
	if not os.path.isfile(emdir+'rhoproj_%s.fits' % filenum):
		make_fits(filenum, 'density', outputdir=rhoprojdir, proj=True, image_res=image_res)

#4 - shift, rotate, compare
def compare_arrays(fitsdir, obsfile, errfile):
	simfiles, times, data, errorsq, x, y, xp = init(obsfile, errfile, fitsdir)
	for file in range(len(simfiles)): #be careful here. this works for now because I'm only generating FITS files for shortlist snaps
			#but ideally we should match the snap # so this doesn't collapse if the shortlist changes somehow
		image = fits.getdata(simfiles[file])
		chisq = shift_rotate_compare(data, image, offset, rot_angle[file], shift[file])

"So this is how far i need to get today"
"make all the fits files and compare them"
"let's move to work, where the internet is much faster. plus i get to take a walk"
#5 - find features

