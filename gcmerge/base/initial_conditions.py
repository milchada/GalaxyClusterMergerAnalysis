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

def rho_s_NFW(M, c): #Msun/kpc**3
    rs = rs_NFW(M, c)
    return M*u.Msun/(4*np.pi*(rs**3)*(np.log(1+c) - c/(1+c)))

def c_sNFW(c): #Lilley eq ??
    return 1.36 + 0.76*c

def a(rvir, c_sNFW): #Lilley 2018 eq 11 and the line below.
    rs = rvir/c_sNFW
    return 3*rs/2

def main(M200, c, fgas=0.17, npart=1, alpha=1.0, rs_rvir=0.5, rc_rvir = 0.1, ptcls=True, suffix=''):
    npart *= 20000000
    rsc = rs_NFW(M200, c).value
    rvir = rsc * c
    cs = c_sNFW(c)
    a_sNFW = a(rvir, cs)
    mass = cg.snfw_mass_profile(1, a_sNFW)
    M = M200*(1.0-fgas)/mass(rvir)
    dm = cg.snfw_density_profile(M, a_sNFW)
    dm = cg.snfw_density_profile(M*(1.0-fgas), a_sNFW)
    rsg = rvir * rs_rvir
    rcg = rvir * rc_rvir
    beta = 2.0 / 3.0
    eps = 2.0
    gas = cg.vikhlinin_density_profile(1.0, rcg, rsg, alpha, beta, eps)
    Mgas = fgas * M200
    
    gas = cg.rescale_profile_by_mass(gas, Mgas, rvir)
    tot = dm + gas
    hse = cg.HydrostaticEquilibrium.from_dens_and_tden(0.1, 30000., gas, tot)
    hse.write_model_to_h5('sNFW_profile_m%0.1e_c%0.1f_a%0.1f%s.h5' % (M200, c, alpha, suffix), in_cgs=True, overwrite=True)
    if ptcls:
        vir = cg.VirialEquilibrium.from_hse_model(hse)
        parts = vir.generate_particles(int(npart))
        parts.write_particles_to_h5('sNFW_%0.1e_particles_c%0.1f_a%0.1f%s.h5' % (M200, c, alpha, suffix), in_cgs=True, overwrite=True)

if __name__=="__main__":
    masses = [4e14, 5e14, 6e14, 7e14, 8e14, 9e14, 1e15, 1e14, 1.6e15, 2.1e14, 2.4e14, 2.7e14]
    for m200 in masses:
        main(m200, 5, alpha=2, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')
    for c in [4, 5, 6]:
        main(6e14, c, alpha=2, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')
        main(2.1e14, c, alpha=2, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')

    #NCC
    main(6e14, 5, alpha=0, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')
    main(2.1e14, 5, alpha=0, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')
    #WCC
    main(6e14, 5, alpha=1, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')
    main(2.1e14, 5, alpha=1, rs_rvir = 0.6, rc_rvir = 0.1, suffix = 'rs0.6rvir_rc0.1rvir_beta0.67_eps2')
    



