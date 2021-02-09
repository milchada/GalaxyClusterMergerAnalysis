import cluster_generator as cg
import astropy.units as u
from astropy import cosmology
import numpy as np
f = cosmology.default_cosmology
f = f.get()
z = 0.34 # 1Gyr before z = 0.2323
rho_crit = f.critical_density(z).to('Msun kpc**-3')

def rs_NFW(M, c): #kpc 
    return (3*M*u.Msun/(800*np.pi*rho_crit*(c**3)))**(1./3)

def a(rs, c): #Lokas 2001 Eq. 28
    rh = rs * c * (0.6082 - 0.1843 * np.log10(c) - 0.1011 * np.log10(c) ** 2 + 0.03918 * np.log(c) ** 3)
    return rh / 5.478

def rho_s_NFW(M, c): #Msun/kpc**3
    rs = rs_NFW(M, c)
    return M*u.Msun/(4*np.pi*(rs**3)*(np.log(1+c) - c/(1+c)))

def c_sNFW(c):
    return 1.36 + 0.76*c

def main(M200, c, fgas=0.17, npart=1, alpha=1.0, rsfrac=0.5, rc_rs = 0.1, ptcls=True, suffix=''):
    npart *= 20000000
    rsc = rs_NFW(M200, c).value
    a_sNFW = a(rsc, c)
    rvir = rsc * c
    # cs = c_sNFW(c)
    mass = cg.snfw_mass_profile(1, a_sNFW)
    M = M200*(1.0-fgas)/mass(rvir)
    dm = cg.snfw_density_profile(M, a_sNFW)
    dm = cg.snfw_density_profile(M*(1.0-fgas), a_sNFW)
    rsg = rsc * rsfrac
    rcg = rc_rs * rsg
    beta = 2.0 / 3.0
    eps = 3.0
    gas = cg.vikhlinin_density_profile(1.0, rcg, rsg, alpha, beta, eps)
    Mgas = fgas * M200
    
    gas = cg.rescale_profile_by_mass(gas, Mgas, rvir)
    tot = dm + gas
    hse = cg.HydrostaticEquilibrium.from_dens_and_tden(0.1, 30000., gas, tot)
    hse.write_model_to_h5('tNFW_profile_m%0.1e_c%0.1f_a%0.1f%s.h5' % (M200, c, alpha, suffix), in_cgs=True, overwrite=True)
    if ptcls:
        vir = cg.VirialEquilibrium.from_hse_model(hse)
        parts = vir.generate_particles(int(npart))
        parts.write_particles_to_h5('tNFW_%0.1e_particles_c%0.1f_a%0.1f%s.h5' % (M200, c, alpha, suffix), in_cgs=True, overwrite=True)

if __name__=="__main__":
	main(1.6e14, 5.2, alpha=0.0, rsfrac = 1.0, rc_rs = 0.5, ptcls=True,suffix='_rs1.0_rc0.5',npart=0.25)
