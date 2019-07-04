from find_peaks import find_peak

def offset_tilt(potfile, xmin=400, xmax=600, ymin = 400, ymax = 600):
	simpeaks =  find_peak(potfile, xmin, xmax, ymin, ymax)	
	#x, y pixel coordinates of 2 potential minima in potfile
	offset = (x, y)
	return offset, rot_angle