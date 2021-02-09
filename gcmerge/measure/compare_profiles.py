import numpy as np 
import matplotlib.pylab as plt
from matplotlib import colors, cm
from astropy.io import fits
from astropy.constants import k_B
import glob, os

home = os.getcwd()
kB = k_B.to('keV/K').value

def plot(bowdat, upstreamdat, fitsfile, ax, ax1, rot=False, color=None, res=None, cmap=cm.tab20c):
	if not res:
		res = fits.getheader(fitsfile)['CDELT1']
	b = np.genfromtxt(bowdat)
	u = np.genfromtxt(upstreamdat) #check units tho
	rb = b[:,0]*res
	tb = b[:,1]*(res**2)*kB
	ru = u[:,0]*res
	tu = u[:,1]*(res**2)*kB
	bmax = np.argmax(tb)
	umax = np.argmax(tu)
	if rot:
		theta = int(bowdat.split('eta')[1].split('_')[0])
		phi = int(bowdat.split('phi')[1].split('_')[0])
		lss = {0:'solid', 15:'dashed', 30:'dotted', 60:'dashed'}
		norm = colors.Normalize(vmin=0, vmax=90)
		m = cm.ScalarMappable(norm=norm, cmap=cmap)

		ax.plot(rb[bmax:] - rb[bmax], tb[bmax:], c=m.to_rgba(phi), linestyle = lss[theta], label = r'$\theta$ = %d, $\phi$ = %d' % (theta, phi))
		ax1.plot(ru[umax:]- ru[umax], tu[umax:], c=m.to_rgba(phi), linestyle = lss[theta], label = r'$\theta$ = %d, $\phi$ = %d' % (theta, phi))
	else:
		ax.plot(rb[bmax:] - rb[bmax], tb[bmax:], c=color, label = bowdat.split('_bow')[0])
		ax1.plot(ru[umax:]- ru[umax], tu[umax:], c=color, label = upstreamdat.split('_upstream')[0])

def plotdir(dirname, rot=False, res=None):
	os.chdir(home+dirname)
	fitsfiles = glob.glob('*fits')
	bow = glob.glob('*bow.dat')
	upstream = glob.glob('*upstream.dat')
	fitsfiles.sort(); bow.sort(); upstream.sort()
	# print(fitsfiles, bow, upstream)
	fig, ax = plt.subplots()
	fig1, ax1 = plt.subplots()
	
	if np.mod(len(bow),2):
		cmap = cm.ScalarMappable(cmap=cm.tab10c, norm=colors.Normalize(0,len(bow))) #to avoid a white colour in the middle
	else:
		cmap = cm.ScalarMappable(cmap=cm.tab10c, norm=colors.Normalize(0,len(bow)+1))
	
	for i in range(len(bow)):
		if rot:
			plot(bow[i], upstream[i], fitsfiles[i], ax=ax, ax1=ax1, rot=True, res=res)
		else:
			plot(bow[i], upstream[i], fitsfiles[i], color=cmap.to_rgba(i), ax=ax, ax1=ax1, rot=False,res=res)
	
	for a in (ax, ax1):
		hand, lab = a.get_legend_handles_labels()
		a.legend(hand, lab)
		a.set_xlim(0,60)
		a.set_ylabel('T (keV)')
		a.set_xlabel('R (kpc)')
	fig.savefig('bowshocks_temp.png')
	fig1.savefig('upstreamshocks_temp.png')
