####################################################################################
# Given an island of points, select those with a minimum contrast vs neighbourhood #
####################################################################################

import numpy as np
from scipy.ndimage import gaussian_gradient_magnitude

def find_points_above(img,feature, mincontrast=1, gwidth=1, type='sb'):
    points = np.array([feature.lines[0].pixlist[0].rawx, feature.lines[0].pixlist[0].rawy])
    for line in feature.lines:
        for pix in line.pixlist:
            points = np.insert(points, -1, (pix.rawx,pix.rawy))
    points = points[1:-1]
    points = np.reshape(points, (int(len(points)/2), 2))
    if type == 'ggm':
        ggm = gaussian_gradient_magnitude(img, gwidth)
        ggms = []
        for point in points:
            ggms.append(ggm[point[0]][point[1]])
        ggms = np.array(ggms)
        return points[np.argwhere( ggms >= mincontrast*ggms.max())]
    elif type == 'sb':
        sbs = []
        for point in points:
            sbs.append(img[point[0]][point[1]])
        sbs = np.array(sbs)
        return points[np.argwhere( sbs >= mincontrast*sbs.max())]
    else:
        print('Type must be either sb or ggm')