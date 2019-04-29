import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
import glob

xraysb = glob.glob('xray_sb/*')
temps = glob.glob('tempproj/*')
xraysb.sort()
temps.sort()
nfiles = min(len(xraysb), len(temps))

#time period where merger active
minfile = 10
maxfile = 50

for file in range(minfile, min(maxfile,nfiles)):
	try:
		sb = fits.getdata(xraysb[file])
		temp = fits.getdata(temps[file])
		pseudo_p = np.sqrt(sb)*temp 
		np.save('pseudopressure/pseusopressure_%d' % file, pseudo_p)
		print("File %d done!" % file)
	except IOError:
		continue