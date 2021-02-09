#make movie of simulation
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
from astropy.io import fits

simname = "1to3b250" #Edit as needed

fig = plt.figure()

files = glob.glob('tempproj/*fits')#('fitsfiles/temp/*fits')
# files = glob.glob('xray_sb/*fits')#('fitsfiles/xray_sb/*fits')
files.sort()

if 'temp' in files[0]:
	norm = plt.Normalize(0, 20)
elif 'xray_sb' in files[0]:
	sbmax = fits.getdata(files[15]).max()
	norm = plt.LogNorm(sbmax/1e4, sbmax)

x = np.arange(fits.getdata(files[0]).shape[0])
y = np.arange(fits.getdata(files[0]).shape[0]).reshape(-1, 1)
ims = []
for add in np.arange(len(files)):
	try:
		ims.append((plt.pcolor(x, y, fits.getdata(files[add]), norm=norm),))
		print add, " finished"
	except (IndexError, IOError, ValueError):
		print add, " failed"
		continue

im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
                                   blit=True)
# To save this second animation with some metadata, use the following command:
im_ani.save('%s.mp4' % simname)

plt.show()