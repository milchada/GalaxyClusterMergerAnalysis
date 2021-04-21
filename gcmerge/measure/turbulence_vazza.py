import numpy as np
import yt
from astropy.convolution import convolve, Box1DKernel

ds = yt.load('Data_000156')
halfwidth = 0.5 #Mpc
_, c = ds.find_min(("gamer","Pote"))
grid = ds.r[(c[0]-halfwidth,"Mpc"):(c[0]+halfwidth,"Mpc"):256j,
            (c[1]-halfwidth,"Mpc"):(c[1]+halfwidth,"Mpc"):256j,
            (c[2]-halfwidth,"Mpc"):(c[2]+halfwidth,"Mpc"):256j]
vx = grid["gas", "velocity_x"].to_value('km/s')				#3D array of velocity		
vy = grid["gas", "velocity_y"].to_value('km/s')
vz = grid["gas", "velocity_z"].to_value('km/s')

#needed parameters and thresholds
n     = 256 						#linear size of the grid
r2    = int(n*0.5-1.) 				#upper limit for L
r1    = 4 							#lower limit
turbo = np.zeros([n,n,n])			#turbulent field
scale = np.zeros([n,n,n]) 			#scale of the flow
sk    = np.ndarray([n,n,n]) 		#skewness
eps   = 0.1 						#tolerance in Eq.5
nk    = 8 							#number of cells to compute skewness
epssk = 1. 							#tolerance for the skewness
drr   = 1. 							#radial step for Eq.5
# preliminary computation of the skewness
meanv 		  = convolve(vel, Box1DKernel(nk), boundary=None)
sc   		  = abs((vx-meanv)/vel)
kernel		  = np.ones([nk,nk,nk])
kernel[0,:,:] = 0.
kernel[:,0,:] = 0.
kernel[:,:,0] = 0.
kernel[nk-1,:,:] = 0.
kernel[:,nk-1,:] = 0.
kernel[:,:,nk-1] = 0.
sk  = convolve((vel-meanv)**2,kernel, boundary=None)
sk  = meanv**3/float(sk**1.5) 				#skewness, Eq.7
sc1 = 0 #..

#.iterations to constrain turbulence
for r in np.arange(r1,r2,drr):
	width = 2.*r+1 											# width of box
	meanv = convolve(vel, Box1DKernel(width), boundary=None)
	#mean local velocity at each scale
	sc  = abs((vel-meanv)/vel) 								#differential change in vel, Eq.5
	skm = convolve(sk,Box1DKernel(width), boundary=None) 	#average skewness within L
	#check of which cells are converged
	############### JOHN PLEASE CHECK THIS ########################
	ibox = np.where((abs(sc-sc1)/float(sc1) > eps) or
			(abs(skm) > epssk) and scalex = 0.)

	if len(ibox) > 0:
		turbo[ibox] = vel[ibox] - meanv[ibox]
		#turbulent velocity in the cell
		scale[ibox] = float(r+0.01) 						#outer scale L
	################################################################	
	sc1 = sc
	############### JOHN PLEASE CHECK THIS ########################
	# the next 3 lines started with a ; so I assume they had to be commented out, but I'm not sure:
	# ;..zc = n*0.5
	# ;tvscl,[vel(*,*,zc),meanv(*,*,zc),
	# turbo(*,*,zc)]
	################################################################	
#saves our final results
np.save("turb.dat", turb)
#outer scale
np.save("scale.dat", scale)
#skewness
np.save("skewness.dat", sk)
