import unyt
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

def hse(M200, c, alpha, fgas, beta, eps, rs_rvir, rc_rvir):
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
    gas = cg.vikhlinin_density_profile(1.0, rcg, rsg, alpha, beta, eps)
    Mgas = fgas * M200
    
    gas = cg.rescale_profile_by_mass(gas, Mgas, rvir)
    tot = dm + gas
    hse = cg.HydrostaticEquilibrium.from_dens_and_tden(0.1, 30000., gas, tot)
    return hse 

def main(M1, M2, c1, c2, fgas=0.17, npart=1.25, alpha1=2.0, alpha2=2.0, rs_rvir=0.6, rc_rvir = 0.1, suffix='', 
    beta = 2./3, eps=2., dclusters=3000., bclusters=100., boxsize=14., vrel=1152.):
    hse1 = hse(M1, c1, alpha1, fgas, beta, eps, rs_rvir, rc_rvir)
    hse2 = hse(M2, c2, alpha2, fgas, beta, eps, rs_rvir, rc_rvir)
    profile1 = 'sNFW_profile_m%0.1e_c%0.1f_a%0.1f%s.h5' % (M1, c1, alpha1, suffix)
    profile2 = 'sNFW_profile_m%0.1e_c%0.1f_a%0.1f%s.h5' % (M2, c2, alpha2, suffix)
    hse1.write_model_to_h5(profile1, in_cgs=True, overwrite=True)
    hse2.write_model_to_h5(profile2, in_cgs=True, overwrite=True)
    
    d = dclusters # distance between the clusters in kpc
    b = bclusters # impact parameter in kpc
    center = np.array([boxsize*500]*3) # center coordinates of the two clusters; converts boxsize in Mpc to midpoint in kpc

    # this function returns two center arrays for the two clusters
    center1, center2 = cg.compute_centers_for_binary(center, d, b)

    # You can define the velocities any way you want, but here’s one way

    # relative velocity
    velocity = np.array([(vrel*unyt.km/unyt.s).to_value("kpc/Myr"), 0.0, 0.0])
    # R is the mass ratio of the clusters
    velocity1 = velocity/(1.0+R)
    velocity2 = -velocity*R/(1.0+R)

    # Now we set up a ClusterICs class that generates the cluster ICs
    basenm = "" # base name of files that the ClusterICs class will write
    profile_files = [profile1, profile2] # filenames of the two HSE profiles
    num_clusters = 2 # number of clusters in the ICs
    num_particles = {"dm": npart*10000000} # the total number of DM particles between the two clusters to use
    r_max = boxsize*500. # maximum radius of particles for each cluster in npc

    ics = cg.ClusterICs(basenm, num_clusters, profile_files,
        [center1, center2], [velocity1, velocity2],
        num_particles=num_particles, r_max=r_max
    )

    # Now we set up GAMER ICs. This function generates the particle files and prints out
    # what you should put in Input__TestProb

    cg.setup_gamer_ics(ics, regenerate_particles=False)

if __name__=="__main__":
    M1 = 5e14
    M2 = 1.6e14
    c1 = c2 = 5
    alpha1 = alpha2 = 2.
    d = 3000
    b = 100
    vrel = 1152.
    boxsize_Mpc = 14.
    main(M1, M2, c1, c2, fgas=0.17, npart=1.25, alpha1=alpha1, alpha2=alpha2, rs_rvir=0.6, rc_rvir = 0.1, suffix='', 
    beta = 2./3, eps=2., dclusters=d., bclusters=b., boxsize=boxsize_Mpc, vrel=vrel)