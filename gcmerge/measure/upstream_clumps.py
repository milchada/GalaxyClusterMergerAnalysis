import numpy as np
import matplotlib.pylab as plt
import yt


#overplot velocity vectors
ds = yt.load('Data_000014') #or something

#select a smaller region around shock
p = yt.ProjectionPlot(ds,'z',("gas","kT"), center='c', width=(1, 'Mpc'), 
	cmap='hot', weight_field="mazzotta_weighting")
p.annotate_quiver(("gas","velocity_x"),("gas","velocity_y"), 
	factor=8, plot_args={'color':'white'})