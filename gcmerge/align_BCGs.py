from bcg_dist import bcg_separation
from scipy.ndimage.interpolation import rotate
import math

def angle(pt1, pt2):
	relative_pos = (pt2 - pt1)
	if len(relative_pos) == 1:
		relative_pos = (pt2 - pt1)[0]
	return np.angle(relative_pos[0] + 1j*relative_pos[1]) #output in radians

def rotate_image(image, angle):
	image[np.isnan(image)] = 0 #otherwise rotation gets fucked
	rotimage = rotate(image, angle = angle)
	datasize = image.shape[0]
	rotsize = rotimage.shape[0]

	if rotsize != datasize:
		xstart = (rotsize - datasize)//2
		rotimage = rotimage[xstart:xstart+datasize, xstart:xstart+datasize]
		rotimage[rotimage == 0] = np.nan
	return rotimage

def align_bcgs(potfile, xrayfile):
	peaks = bcg_separation(potfile, ret_peaks=True, xmin=300, xmax=1200, ymin = 300, ymax = 1200)
	potential = fits.getdata(potfile)
	minima1 = potential[peaks[0][0]][peaks[0][1]] 
	minima2 = potential[peaks[1][0]][peaks[1][1]] 

	#this is centered on xray peak
	obsdata, errorsq, x, y, xp = init(obsfile, errfile)

	#currently manual because wtf
	bcg1_pix = xp
	bcg2_pix = (270,230)

	data = fits.getdata(xrayfile)
	header = fits.getheader(xrayfile)
	dx = header['CDELT1']
	sx = calibrate(data.shape[1],header,axis=1)
	sy = calibrate(data.shape[1],header,axis=2)

	if minima1 < minima2:
		peak1 = peaks[0]
		peak2 = peaks[1]
	else:
		peak1 = peaks[1]
		peak2 = peaks[0]

	pos_peak1 = (sx[peak1[0]], sy[peak1[1]])

	#peak1 should = xp = center
	sx -= pos_peak1[0]
	sy -= pos_peak1[1]

	#rotate each about xp counterclockwise by angle 
	bcg_angle_obs = np.degrees(angle(bcg2_pix - xp, bcg1_pix - xp))
	bcg_angle_sim = np.degrees(angle(peak2, peak1))

	f = interp2d(sx, sy, data)
	binned_data = f(x,y)

	rot_image = rotate(binned_data, bcg_angle_sim-bcg_angle_obs, reshape=False)
	
	#good! now output this shifted and rotated array
	return np.flip(rot_image, axis=1) #because obs RA are decreasing, while sim is increasing

#testing in app editing