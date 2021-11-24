from basicplots import *

def fig1(xmin=0, xmax=200,ymin=15,ymax=215, tickmax=300):
	snaps = [164, 168, 170, 180]
	dir = 'beta=100/'
	fig, ax = plt.subplots(ncols = 4, nrows = 4, figsize=(8,8))
	fields = ['photon_emissivity', 'temperature', 'magnetic_field', 'beta']
	for i in range(len(fields)):
		for j in range(len(snaps)):
			plot(tfile=get_tfile(dir, snaps[j], fields[i]), ax = ax[i][j], field=fields[i], 
				cmap=cmaps[fields[i]], norm = norms[fields[i]], xmin=0.3, xmax=0.75)
	for a in ax.flatten():
		a.set_xticks([]); a.set_yticks([])
	tfile = get_tfile(dir,snaps[0])
	for i in range(4): 
		axisticks(ax[i][0], tfile=tfile, xmin=ymin, xmax=ymax, x=False, tickmax=tickmax)
	for i in range(4): 
		axisticks(ax[3][i], tfile=tfile, xmin=xmin, xmax=xmax, y=False, ntix=3, tickmax=tickmax)
	plt.tight_layout()
	return fig, ax 


def fig2(xmin=0, xmax=200,ymin=15,ymax=260, tickmax=280):
	times = [0, 170, 174, 180]
	sims = ['beta=200/','beta=100/','beta=50/']
	field = 'rm'
	plt.rcParams.update({'font.family':'Geneva'})
	fig, ax = plt.subplots(nrows=3, ncols=4, sharey=True)
	for i in range(3):
		for j in range(4):
			# print(sims[i], times[j])
			tfile = get_tfile(sims[i], times[j], field)
			# plot(ax[i][j], dir=sims[i], snapnum	= times[j], field=field, xmin=0.2, xmax=0.75,
			#  cmap=cmaps[field], norm=norms[field])
			ax[i, j].text(xmin+10, ymin+10,  r'$\sigma_{RM}$ = %d rad m$^{-2}$' % np.std(fits.getdata(tfile)), fontsize=10) 
			if not i:
				ax[i][j].set_title('%0.2f Gyr' % (times[j]/100.))
	for a in ax.flatten():
		a.set_xticks([]); a.set_yticks([])
		a.set_xlim(xmin, xmax); a.set_ylim(ymin, ymax)
	tfile = get_tfile(sims[0],times[0])
	for i in range(3): 
		axisticks(ax[i][0], tfile=tfile, xmin=ymin, xmax=ymax, x=False)
	for i in range(4): 
		axisticks(ax[2][i], tfile=tfile, xmin=xmin, xmax=xmax, y=False, ntix=3, tickmax=300)
	return fig, ax 

def fig3(xmin=0, xmax=220, ymin=60, ymax=280 ):
	snapnum=180 
	dirs = ['beta=200/', 'beta=100/', 'beta=50/']
	fig, ax = plt.subplots(ncols = 4, nrows = 3, figsize=(8,6))
	fields = ['photon_emissivity', 'temperature', 'magnetic_field', 'beta']
	for i in range(len(dirs)):
		for j in range(len(fields)):
			
			plot(tfile=tfile, ax = ax[i][j], field=fields[j], 
			cmap=cmaps[fields[j]], norm = norms[fields[j]], xmin=0.2, xmax=0.75)
	for a in ax.flatten():
		a.set_xticks([]); a.set_yticks([])
		a.set_xlim(xmin, xmax); a.set_ylim(ymin, ymax)
	tfile = get_tfile(dirs[0],snapnum)
	for i in range(3): 
		axisticks(ax[i][0], tfile=tfile, xmin=xmin, xmax=xmax, x=False)
	for i in range(4): 
		axisticks(ax[2][i], tfile=tfile, xmin=xmin, xmax=xmax, y=False, ntix=3, tickmax=300)
	return fig, ax 

def fig4(xmin=0, xmax=230,ymin=60,ymax=280):
	time = 180
	fig, ax = plt.subplots(nrows=3, ncols=3, sharey=True)
	sims = ['beta=inf/', 'beta=50/', 'beta=50/turnoff_85/']
	fields = ['photon_emissivity', 'ggm', 'temperature']
	for i in range(3):
		for j in range(3):
			print(i,j)
			try: 
				plot(ax[i][j], dir=sims[i], snapnum	= time, field=fields[j], xmin=0.2, xmax=0.75)
			except ValueError:
				continue
	for a in ax.flatten():
		a.set_xticks([]); a.set_yticks([])
		a.set_xlim(xmin, xmax); a.set_ylim(ymin, ymax)
	tfile = get_tfile(sims[0],time)
	for i in range(3): 
		axisticks(ax[i][0], tfile=tfile, xmin=ymin, xmax=ymax, x=False)
	for i in range(3): 
		axisticks(ax[2][i], tfile=tfile, xmin=xmin, xmax=xmax, y=False, ntix=3, tickmax=300)
	return fig, ax 


def fig5():
	run power_spectrum.py

#fig 6 and fig 7 are done in ds9
#make sure to use a snapshot that's in the other figs, probably 180

def fig8(xmin=0,xmax=200, ymin=15,ymax=215):
	fig, ax = plt.subplots(ncols=3,sharey=True)
	sim = 'beta=100/'
	filebase = 'fitsfiles/rm/rm_proj_180_theta'
	ends = ['0phi0.fits', '30phi0.fits', '30phi30.fits']
	for i in range(3):
		a = ax.flatten()
		plot(a, tfile=sim+filebase+ends[i], xmin=0.3, xmax=0.75, field='rm')
		axisticks(a, tfile=tfile, xmin=xmin, xmax=xmax, y=False, ntix=3)

def restest(xmin=0,xmax=230):
	fig, ax = plt.subplots(nrows=3, ncols=3,sharey=True)
	time = 170 #but later make this 176
	sims = ['beta=100/LowRes/','beta=100/','beta=100/HiRes/']
	fields = ['photon_emissivity', 'temperature', 'beta']
	for i in range(3):
		for j in range(3):
			plot(ax=ax[i][j], dir=sims[i], snapnum	= time, field=fields[j], xmin=xmin, xmax=xmax, norm=norms[fields[j]],cmap=cmaps[fields[j]])
	for a in ax.flatten():
		a.set_xticks([]); a.set_yticks([])
	tfile = get_tfile(sims[0],time)
	for i in range(3): 
		axisticks(ax[i][0], tfile=tfile, xmin=xmin, xmax=xmax, x=False)
	for i in range(3): 
		axisticks(ax[2][i], tfile=tfile, xmin=xmin, xmax=xmax, y=False, ntix=3)
	ax[0][0].text(20,180,'max_level = 3', color='w', fontsize=10)  
	ax[1][0].text(20,180,'max_level = 4', color='w', fontsize=10)  
	ax[2][0].text(20,180,'max_level = 5', color='w', fontsize=10)  
	return fig, ax
