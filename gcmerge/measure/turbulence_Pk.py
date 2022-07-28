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

def plot(ds, c = [7., 7., 7.], halfwidth = 0.25, ncell=256j, ax = None, linestyle='solid', retimg=False, x=False, y=False, z = False, all=True, color='k', exp=0, field='velocity', units='km/s', label='v'):
    grid = ds.r[(c[0]-halfwidth,"Mpc"):(c[0]+halfwidth,"Mpc"):ncell,
                (c[1]-halfwidth,"Mpc"):(c[1]+halfwidth,"Mpc"):ncell,
                (c[2]-halfwidth,"Mpc"):(c[2]+halfwidth,"Mpc"):ncell]
    width_kpc = halfwidth*2000.
    if not ax:
        fig, ax = plt.subplots(figsize=(10,10))
    if x:
        vx = grid["gas", field+"_x"].to_value(units)
        kx, vxk = image_powerspec(vx, width_kpc, width_kpc, width_kpc)
        ax.loglog(kx, vxk*pow(kx,exp), label=label+"$_x$", linestyle = linestyle, c=x)
    if y:
        vy = grid["gas", field+"_y"].to_value(units)
        ky, vyk = image_powerspec(vy, width_kpc, width_kpc, width_kpc)
        ax.loglog(ky, vyk*pow(ky,exp), label=label+"$_y$", linestyle = linestyle, c=y)
    if z:
        vz = grid["gas", field+"_z"].to_value(units)
        kz, vzk = image_powerspec(vz, width_kpc, width_kpc, width_kpc)
        ax.loglog(kz, vzk*pow(kz,exp), label=label+"$_z$", linestyle = linestyle, c=z)
    if all:
        vmag = grid["gas", field+"_magnitude"].to_value(units)
        k, vk = image_powerspec(vmag, width_kpc, width_kpc, width_kpc)
        ax.loglog(k, (vk)*pow(k,exp), label=label+"$_{\rm mag}$", c=color)
    
    ax.set_xlabel(r"k (kpc$^{-1}$)", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if retimg:
        return fig, ax
    
def overtime(dir, ax, halfwidth = 0.5, ncell=512j, x=False, y=False, z = False, all=True, exp=0, field='velocity', units='km/s', label='v',
    snaps = ['Data_000000', 'Data_000085', 'Data_000170', 'Data_000180', 'Data_000200']):
    colors = ['tab:blue', 'tab:green', 'tab:orange', 'tab:red', 'tab:purple']
    for s, col in zip(snaps, colors[:len(snaps)]):
        ds = yt.load(dir+'/'+s)
        _, c = ds.find_min(("gamer", "Pote"))
        plot(ds, ax=ax, c = c.value, linestyle='solid', color=col, field=field, units=units, label=label)

fig, ax = plt.subplots()
overtime('beta=50',ax=ax, field='magnetic_field', units='G', label='B')
# ax.vlines(1/6.8, 1e3, 1e10, color='k', linestyle='dotted')
ax.set_xlim(.008,.1)
# ax.set_ylim(1e2,1e10)
plt.ylabel(r"P(k)$\times k^3$", fontsize=14)
plt.tight_layout()
fig.savefig('/home/uc24/PkB_t_beta50.png')
ax.cla()
overtime('beta=inf',ax=ax, field='magnetic_field', units='G', label='B')
# ax.vlines(1/6.8, 1e3, 1e10, color='k', linestyle='dotted')
ax.set_xlim(.008,.1)
# ax.set_ylim(1e2,1e10)
plt.ylabel(r"P(k)$\times k^3 $", fontsize=14)
plt.tight_layout()
fig.savefig('/home/uc24/PkB_t_hydro.png')
ax.cla()
overtime('beta=50/turnoff_85',ax=ax, field='magnetic_field', units='G', label='B')
# ax.vlines(1/6.8, 1e3, 1e10, color='k', linestyle='dotted')
ax.set_xlim(.008,.1)
# ax.set_ylim(1e2,1e10)
plt.ylabel(r"P(k)$\times k^3$", fontsize=14)
plt.tight_layout()
fig.savefig('/home/uc24/PkB_t_turb.png')

def main(dirs, linestyles, figname,snapname='/Data_000180', c = 'potmin', halfwidth = 0.25, retimg=False, all=True, exp = 0, ymin=0.2,ymax=1e12, field='velocity', units='km/s', label='v'):
    fig, ax = plt.subplots(figsize=(10,10))
    for (dir, ls) in zip(dirs, linestyles):
        ds = yt.load(dir+snapname) #make sure this file is correctly named in each dir
        if c == 'potmin':
            _, c = ds.find_min(("gamer","Pote"))
            c = c.value
        elif c == 'c':
            c = [7., 7., 7.]
        if nd:
            plot(ds, ax=ax, all=all, linestyle=ls, c = c, halfwidth = halfwidth, retimg=retimg, exp=exp, field=field, units=units, label=label)
        else:
            plot(ds, ax=ax, all=all, color=ls, c = c, halfwidth = halfwidth, retimg=retimg,exp=exp, field=field, units=units, label=label)
    ax.set_xlim(3e-3,0.5)
    ax.set_ylim(ymin, ymax)
    fig.savefig('/home/uc24/'+figname)
    if retimg:
        return fig, ax

if __name__ == "__main__":
    colors = ['tab:blue', 'tab:green', 'tab:orange']
    linestyles = ['solid', 'dashed', 'dotted']
    dirs = ['beta=50', 'beta=50/turnoff_85', 'beta=inf']
    fig, ax = main(dirs, colors, 'Pk_beta50_k4.png', nd=False, exp = 4, ymax=1e2, snapname='/Data_000180')
    gc.collect(); gc.collect(); gc.collect()
    fig, ax = main(dirs[0], linestyles[0], 'beta50_Pk_B.png', nd=True, exp = 1, ymin=1e-16, ymax=1e-8, snapname='/Data_000180', field='magnetic_field', units='G', label='B')
    gc.collect(); gc.collect(); gc.collect()

    dirs = ['beta=200', 'beta=100', 'beta=50'] #'beta=inf',
    colors = ['tab:green', 'tab:orange', 'tab:red'] #'tab:blue', 
    fig, ax = main(dirs, colors, 'Pk_allbeta.png', nd=False, exp=4, ymax=1e2, snapname='/Data_000180')
    gc.collect(); gc.collect(); gc.collect()
    fig, ax = main(dirs, colors, 'allbeta_Pk_B.png', nd=False, exp = 1, ymin=1e-16, ymax=1e-8, snapname='/Data_000180', field='magnetic_field', units='G', label='B')
    gc.collect(); gc.collect(); gc.collect()
