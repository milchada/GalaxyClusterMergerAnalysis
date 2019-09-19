#######################################################
# Applies an unsharp map in the same way as observers #
#######################################################

from astropy.io import fits      
import matplotlib, gc
matplotlib.use("Agg")                                                
import matplotlib.pylab as plt                                                   
from matplotlib import colors, cm               
from astropy import convolution                                                  
import numpy as np                                                              

def convolve(array, sigma= 50/14.):
        kernel = convolution.Gaussian2DKernel(stddev=sigma)
        return convolution.convolve(array, kernel)

def unsharp_mask(file = 'xray_photon_emissivity_0.3_7_keV_proj_14.fits', 
    sigma_arcsec = np.array([1.5, 5,20]), ret = True):
    data = fits.getdata(file)             
    header = fits.getheader(file)        

    sigma_pix = sigma_arcsec/header['CDELT1']                                       

    im0 = convolve(data, sigma_pix[0])
    im1 = convolve(data, sigma_pix[1])                                              
    im2 = convolve(data, sigma_pix[2])                                              

    norm = colors.LogNorm(data.max()/1e4, data.max())                                

    fig, ax = plt.subplots(ncols=2, sharey=True, sharex=True)                                     
    ax1, ax2 = ax.flatten()                                                      

    a1 = ax1.imshow(im0, norm=norm, cmap=cm.afmhot)                               

    for axis in ax: 
        axis.set_xlim(800,1200) 
        axis.set_ylim(950,1100) 
        xtix = axis.get_xticks()
        ytix = axis.get_yticks()

    xla = (xtix - xtix.mean())*header['CDELT1']
    yla = (ytix - ytix.mean())*header['CDELT1']
    plt.xticks(xtix, ['%d' % x for x in xla])
    plt.yticks(ytix, ['%d' % x for x in yla])

    contr = (im1-im2)/(im1+im2)                                                     

    a2 = ax2.imshow(contr, norm=colors.Normalize(-1., 0.4), cmap=cm.afmhot)                                                                             

    outname = file.split('_')[0]+file.split('_')[-1].split('.')[0]+'.png'
    fig.colorbar(a1, ax = ax1)
    fig.colorbar(a2, ax = ax2)
    fig.tight_layout()
    plt.savefig(outname)
    print (outname, 'done!')
    gc.collect()
    if ret:
        return fig, ax

import glob
files = glob.glob('*fits')
files.sort()
for file in files:
    unsharp_mask(file, ret = False)
    plt.close()