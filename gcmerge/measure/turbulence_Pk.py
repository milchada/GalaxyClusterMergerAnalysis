import yt, gc
import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt

def _nd_window(data, filter_function):
    """
    https://stackoverflow.com/questions/27345861/extending-1d-function-across-3-dimensions-for-data-windowing

    Performs an in-place windowing on N-dimensional spatial-domain data.
    This is done to mitigate boundary effects in the FFT.

    Parameters
    ----------
    data : ndarray
           Input data to be windowed, modified in place.
    filter_function : 1D window generation function
           Function should accept one argument: the window length.
           Example: scipy.signal.hamming
    """
    for axis, axis_size in enumerate(data.shape):
        # set up shape for numpy broadcasting
        filter_shape = [1, ] * data.ndim
        filter_shape[axis] = axis_size
        window = filter_function(axis_size).reshape(filter_shape)
        # scale the window intensities to maintain image intensity
        np.power(window, (1.0/data.ndim), out=window)
        data *= window

def image_powerspec(data, Lx, Ly, Lz):
    
    from scipy.signal import hann

    nx, ny, nz = data.shape
    
    dx = Lx/nx
    dy = Ly/ny
    dz = Lz/nz

    # Shift the wavenumbers so that the zero is at the center of the transformed image
    kx = fftpack.fftshift(fftpack.fftfreq(nx, d=dx))
    ky = fftpack.fftshift(fftpack.fftfreq(ny, d=dy))
    kz = fftpack.fftshift(fftpack.fftfreq(nz, d=dz))

    # Compute the 3D grid of wavenumbers

    kx, ky, kz = np.meshgrid(kx, ky, kz, indexing="ij")
    kk = np.sqrt(kx*kx+ky*ky+kz*kz)[:,:,:]

    # Compute the 3D power spectrum

    _nd_window(data, hann)

    P = Lx*Ly*Lz*np.abs(fftpack.fftshift(fftpack.fftn(data)/(nx*ny*nz)))**2

    # Set the maximum and minimum limits on the wavenumber bins

    kmin = 1.0/Lx
    kmax = 1.0/dx

    # Bin up the 3D power spectrum into a 1-D power spectrum

    bins = np.logspace(np.log10(kmin), np.log10(kmax), 100)
    k = np.sqrt(bins[1:]*bins[:-1])
    Pk = np.histogram(kk, bins, weights=P)[0] / np.histogram(kk, bins)[0]

    return k[Pk > 0], Pk[Pk > 0]

def plot(ds, c = [7., 7., 7.], halfwidth = 0.25, ax = None, linestyle='solid', retimg=False, nd=True, color='k', exp=0):
    grid = ds.r[(c[0]-halfwidth,"Mpc"):(c[0]+halfwidth,"Mpc"):256j,
                (c[1]-halfwidth,"Mpc"):(c[1]+halfwidth,"Mpc"):256j,
                (c[2]-halfwidth,"Mpc"):(c[2]+halfwidth,"Mpc"):256j]
    width_kpc = halfwidth*2000.
    if not ax:
        fig, ax = plt.subplots(figsize=(10,10))
    vx = grid["gas", "velocity_x"].to_value('km/s')
    vy = grid["gas", "velocity_y"].to_value('km/s')
    vz = grid["gas", "velocity_z"].to_value('km/s')
    k, vxk = image_powerspec(vx, width_kpc, width_kpc, width_kpc)
    k, vyk = image_powerspec(vy, width_kpc, width_kpc, width_kpc)
    k, vzk = image_powerspec(vz, width_kpc, width_kpc, width_kpc)
    if nd:
        ax.loglog(k, vxk*pow(k,exp), label="v$_x$", linestyle = linestyle, c='tab:blue')
        ax.loglog(k, vyk*pow(k,exp), label="v$_y$", linestyle = linestyle, c='tab:green')
        ax.loglog(k, vzk*pow(k,exp), label="v$_z$", linestyle = linestyle, c='tab:orange')
            
    else:
        #vmag = grid["gas", "velocity_magnitude"].to_value('km/s')
        #k, vk = image_powerspec(vmag, width_kpc, width_kpc, width_kpc)
        ax.loglog(k, (vxk+vyk+vzk)*pow(k,exp), label="v$_{\rm mag}$", c=color)
        
    ax.set_xlabel("k (kpc$^{-1}$)")
    if exp:
        ax.set_ylabel("P(k)$\times$ $k^%d$ (km$^2$ s$^{-2}$ kpc$^{3}$)" % exp)
    else:
        ax.set_ylabel("P(k) (km$^2$ s$^{-2}$ kpc$^{3}$)")
    if retimg:
        return ax
    
def main(dirs, linestyles, figname, c = 'potmin', halfwidth = 0.25, retimg=False, nd=True, exp = 0, ymax=1e12):
    fig, ax = plt.subplots(figsize=(10,10))
    for (dir, ls) in zip(dirs, linestyles):
        ds = yt.load(dir+'/Data_000156') #make sure this file is correctly named in each dir
        if c == 'potmin':
            _, c = ds.find_min(("gamer","Pote"))
            c = c.value
        elif c == 'c':
            c = [7., 7., 7.]
        if nd:
            plot(ds, ax=ax, nd=True, linestyle=ls, c = c, halfwidth = halfwidth, retimg=retimg, exp=exp)
        else:
            plot(ds, ax=ax, nd=False, color=ls, c = c, halfwidth = halfwidth, retimg=retimg,exp=exp)
    ax.set_xlim(3e-3,0.5)
    ax.set_ylim(0.2, ymax)
    fig.savefig('/home/uc24/'+figname)

if __name__ == "__main__":
    colors = ['tab:blue', 'tab:green', 'tab:orange']
    linestyles = ['solid', 'dashed', 'dotted']
    dirs = ['beta=50', 'beta=50/turnoff_at_0_7Gyr', 'beta=inf']
    main(dirs, colors, 'turbulence_comparison_beta50_k4.png', nd=False, exp = 4, ymax=1e5)
    gc.collect(); gc.collect(); gc.collect()
    main(dirs, linestyles, 'turbulence_comparison_3d_beta50_k4.png', nd=True, exp=4, ymax=1e5)
    gc.collect(); gc.collect(); gc.collect()
    
    dirs = ['beta=inf','beta=200', 'beta=100', 'beta=50']
    colors = ['tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    main(dirs, colors, 'turbulence_comparison_allbeta_k4.png', nd=False, exp=4, ymax=1e2)
    gc.collect(); gc.collect(); gc.collect()
    main(dirs, linestyles, 'turbulence_comparison_3d_allbeta_k4.png', nd=True,exp=4, ymax=1e2)
    gc.collect(); gc.collect(); gc.collect()