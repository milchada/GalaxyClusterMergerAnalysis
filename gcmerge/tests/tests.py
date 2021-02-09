########################################
# Test different parts of the pipeline #
########################################

def test_arcfits(island, showpts=True):
	arcfit = fit_arc(island, 1.5)
	fig, ax = plt.subplots()
	sm = cm.ScalarMappable(cmap=cm.seismic, norm=colors.Normalize(0,1))
	if showpts:
		feature = find_points_above_contrast(island, 0)[:,0]
		xdata = feature[:,1]
		ydata = feature[:,0]
		plt.scatter(ydata, xdata, marker='x', c='k')
	for row in range(len(arcfit)):
		cx, cy, r = arcfit[row, 2:5]
		feature = find_points_above_contrast(island, arcfit[row,0])[:,0]
		xdata = feature[:,1]
		ydata = feature[:,0]
		thetas = np.arctan2((ydata - cy),(xdata - cx))
		xfit = r*np.cos(thetas) + cx
		yfit = r*np.sin(thetas) + cy
		plt.scatter(yfit, xfit, c = sm.to_rgba(arcfit[row,0]), lw=0)
	sm.set_array([])
	fig.colorbar(sm)
	plt.xlim(300,600)
	plt.ylim(300,600)

names = name = {1:'Cold Front', 3: 'Bow', 6: 'Swirl', 9: 'Upstream'}
def arcfits_stability(islandnums=[1,3,6,9],names=names, time=1.5):
	fig1, ax1 = plt.subplots()
	for island in islandnums:
		fig, ax = plt.subplots(ncols= 2)
		arcfit = fit_arc(island, 1.5)
		ax1.plot(arcfit[:,0], arcfit[:,-1], label=names[island])
	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels)
	ax1.set_xlabel('Contrast threshold')
	ax1.set_ylabel(r'$\Sigma d_{min}$')
	return fig1, ax1

def test_peak_finding():
	files = glob.glob('potential/*fits')
	files.sort()
	fig, ax = plt.subplots(nrows=3,ncols=3,sharex=True, sharey=True,start=9,end=18)
	peaks = np.empty([end-start+1, 4])
	for axis, file in zip(ax.flatten(), files[start:end]):
		try:
			find_peak(file, axis)
		except(IOError,IndexError,ValueError):
			continue
		fig.savefig('potential_z.png')