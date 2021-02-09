def cut(img, centre, halfwidth = 256, peak_threshold = 0.9):
    #set left and right limits, accounting for possibility that 
    #this may extend beyond box edge

	if (centre[0]-halfwidth) < 0:
		left = 0
		right = 2*halfwidth

	elif (centre[0]+halfwidth) > img.shape[0]:
		right = img.shape[0]
		left = right - 2*halfwidth

	else:
		left, right = centre[0] - halfwidth, centre[0] + halfwidth

    #set top and bottom limits, accounting for possibility that 
    #this may extend beyond box edge
	if (centre[1]-halfwidth) < 0:
		bottom = 0
		top = 2*halfwidth

	elif (centre[1]+halfwidth) > img.shape[0]:
		top = img.shape[0]
		bottom = right - 2*halfwidth

	else:
		bottom, top = centre[1] - halfwidth, centre[1] + halfwidth

	return img[left:right, bottom:top]