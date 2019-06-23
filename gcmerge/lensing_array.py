#####################################################################
# Convert lensing observations into arrays for comparison with sims #
#####################################################################

import numpy as np
import pandas as pd

def dataFrame(lensing_cat, zobs, bcg_pos_inpix, columns=None):
	arcsecPerPix = 202./4096 #HST-WFC param
	zobs = 0.2323
	
	data = np.genfromtxt(lensing_cat)
	if data.shape[1] == 6:
		columns = ['x','y','g1','g2','wi','sigma']
	elif data.shape[0] == 4:
		columns = ['x','y','g1','g2']
	else:
		print "Error: Unknown column configuration. Find column names from catalog and input as argument."
		return 0

	lensing_data = pd.DataFrame(data=data, columns = columns)
	mpcPerArcsec = cosmology.Planck15.angular_diameter_distance(zobs)*units.arcsec.in_units('radian')

	lensing_data['x'] -= bcga_pos[0]#lensing_data['x'][peak]
	lensing_data['y'] -= bcga_pos[1]#lensing_data['y'][peak]
	lensing_data['x'] *= arcsecPerPix*mpcPerArcsec
	lensing_data['y'] *= arcsecPerPix*mpcPerArcsec
	
	return lensing_data


