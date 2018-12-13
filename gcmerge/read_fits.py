import numpy as np
import matplotlib.pylab as plt
import matplotlib
from astropy.io import fits
import glob

snaps=glob.glob('xray_emissivity*fits')
snaps.sort()
#goal here is to rescale everything due to mass
#density stays same, R, V scale as M^1/3, T ~ M^2/3

Mobs = 4.15# e14Msun
Msim = 6# e14Msun
scaling = Mobs/Msim

def open_fits(filename, tablename):
	f = fits.open(filename)
	data = f[tablename].data
	header = f[tablename].header
	return data, header

def calibrate(length, header, axis=1):
	refpix = header['CRPIX'+str(axis)]
	refval = header['CRVAL'+str(axis)]
	step = header['CDELT'+str(axis)] 
	values = np.empty(length)
	for ind in xrange(length):
		fitsind = ind+1
		values[ind] = refval+(fitsind-refpix)*step
	return values*scaling**(1./3)

#do this interactively to finetune
def plotmap(data,x,y, vmin, vmax, log = True):#, name, title):
	plt.clf()
	X, Y = np.meshgrid(x, y)
	if log:
		plt.pcolormesh(X, Y,data,norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax))
	else:
		plt.pcolormesh(X, Y,data,norm=matplotlib.colors.Normalize(vmin=vmin,vmax=vmax))
	plt.xticks(np.linspace(x.min(),x.max(),7),['%d' % ind for ind in np.linspace(x.min()/1000,x.max()/1000,7)])
	plt.yticks(np.linspace(y.min(),y.max(),7),['%d' % ind for ind in np.linspace(y.max()/1000, y.min()/1000,7)])
	plt.colorbar()
	plt.xlabel('y (Mpc)')
	plt.ylabel('z (Mpc)')

	# plt.title(title)
	# plt.savefig(name+'.png')
if __name__=="__main__":
	for filename in snaps:
		#general calibration
		xdata, xheader = open_fits(filename, 'XRAY_EMISSIVITY_0.3_7.0_KEV')
		x = calibrate(xdata.shape[1], xheader, axis = 1)
		y = calibrate(xdata.shape[0], xheader, axis = 2)
		plotmap(xdata, x, y, 1e-9, 1e-3) #projecting density adds extra R scaling
		plt.savefig('projy/xray_%s.png' % filename.split('proj_')[1].split('.fits')[0])
		print filename, ' done'
