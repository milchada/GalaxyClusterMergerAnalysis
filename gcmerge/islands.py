#!/usr/bin/env python
"""Turn a pixel list into a list of event islands.
"""

class Pixel (object):
    """Single pixel.
    """
    def __init__ (self, pt):
        # read in sorted array of x, y pixel positions of points
        self.pt = pt
        self.rawx = pt[0]
        self.rawy = pt[1]
        self.island = None

class PixLine (object):
    """List of adjacent pixels above threshold ***with the same rawy***.
    Note that pixlist is sorted by x 
    """
    def __init__ (self, pixlist):
        self.pixlist = pixlist
        self.y = pixlist[0].rawy
        self.xmin = pixlist[0].rawx
        self.xmax = pixlist[-1].rawx

    def __str__ (self):
        """Turn point array into list of ('x','y'). 
        But why have them as strings? This might be a vestige of Paul's input file
        """
        s = "["
        n = 0
        for t in self.pixlist:
            if n:
                s += ","
            n += 1
            s += "({},{})".format (t.rawx, t.rawy)
        s +="]"
        return s

        # PixLines with the same y cannot overlap
    def __lt__ (self, other):    
        return (self.y < other.y) or ((self.y == other.y) and (self.xmax < other.xmin))
        
    def __gt__ (self, other):
        return ( self.y > other.y) or ((self.y == other.y) and (self.xmin > other.xmax))
    
    def __eq__ (self, other):
        return not (self > other) and not (self < other)

    def count (self):
        """Number of pixels.
        """
        return len (self.pixlist)

    def setIsland (self, isle):
        """Set the colour of every pixel to the specified value.
        """
        self.island = isle
        for p in self.pixlist:
            p.island = isle

def mkLines (pixlist):
    """Sort pixel list into PixLines.
    """
    # Pixels in ascending order
    
    linelist = []
    line = None
    ylast = -1
    for pix in pixlist:
        if pix.rawy != ylast or pix.rawx != xlast + 1:
            # Beginning of a new line segment
            if line:
                # Save previous line
                linelist.append (PixLine (line))
            # Start a new line
            line = [pix]
            ylast = pix.rawy
        else:
            line.append (pix)
        xlast = pix.rawx
    if line:
        linelist.append (PixLine (line))
    return linelist

class Island (object):
    """List of directly or indirectly adjacent PixLines.
        Takes in list of horizontal line segments
    """
    def __init__ (self, line):
        """Start with a single PixLine.
        """
        self.lines = [line]
        self.ymax = line.y

    def adjoin (self, island):
        """Merge another island into this one.
        """
        # Ordered merge of the lists of line segments
        lnew = []
        u = self.lines
        v = island.lines
        i = 0
        ilen = len (u)
        j = 0
        jlen = len (v)
        while i < ilen and j < jlen:
            """This part sorts line segments by x value"""
            if u[i] < v[j]:
                lnew.append (u[i])
                i += 1
            else:
                lnew.append (v[j])
                j += 1
        if i < ilen:
            lnew += u[i:]
        else:
            lnew += v[j:]
        self.lines = lnew
        if island.ymax > self.ymax:
            self.ymax = island.ymax

    def adjacent (self, seg, maxsep=1):
        """Test if PixLine "overlaps" with any part of this island.
        """
        # Relies on segment being later than the segments in the island
        yseg = seg.y
        xmin = seg.xmin
        xmax = seg.xmax
        for oseg in self.lines [-1::-1]:
        # Search backward to take advantage of ordering
            if oseg.y == yseg:
            # Line segments are disjoint
                continue
            if oseg.y + maxsep < yseg:
                # No earlier PixLine can be adjacent
                return False
            # Here oseg.y + 1 == yseg
            if xmin - maxsep > oseg.xmax:
                # No remaining PixLines can overlap
                return False
            if xmax + maxsep >= oseg.xmin: 
                return True
        return False

    def count (self):
        """Number of pixels in an island.
        """
        n = 0
        for t in self.lines:
            n+= t.count ()
        return n

    def setIsland (self):
        """Record the Island for every PixLine and Pixel in this Island.
        """
        for l in self.lines:
            l.setIsland (self)


def mkIslands (linelist, maxsep = 1):
    """Construct a list of islands from the list of line segments.
    """
    # Start the list of islands
    islands = [Island (linelist[0])]
    # Add line segments to the map in increasing order
    for lineseg in linelist [1:]:
        # Find all current islands adjacent to the new line segment
        yseg = lineseg.y
        adjacents = []
        for isle in islands [-1::-1]:
            # Check if this island is adjacent
            if isle.ymax < yseg - maxsep:
                # Islands ordered by uppermost, rightmost line segment
                break
            if isle.adjacent (lineseg, maxsep):
                # Add adjacent islands to the adjacent list
                adjacents.append (isle)
        newisle = Island (lineseg)
        for isle in adjacents:
            newisle.adjoin (isle)
            islands.remove (isle)
        islands.append (newisle)
    return islands

import numpy as np
from scipy.ndimage import gaussian_gradient_magnitude

def find_points_above(img, islandlist, island, mincontrast=1, gwidth=1, type='sb'):
    feature = islandlist[island]
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