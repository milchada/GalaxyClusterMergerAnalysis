
##currently using default beta-profile params from vikhlinin, change if well-motivated

from formulas import sNFW_density_profile,     vikhlinin_density_profile,     rescale_profile_by_mass
import yt.units as u
from cluster_generator import HydrostaticEquilibrium, VirialEquilibrium, ClusterModel
from astropy import cosmology
import numpy as np

f = cosmology.default_cosmology
f = f.get()
rho_crit = f.critical_density0.to('Msun kpc**-3')

# The `formulas` package contains a number of astrophysically relevant formulas, including radial profiles. 
# For DM, we'll use the "super-NFW" profile, and for the gas, we'll use Alexey's modifications to the $\beta$-model.
# What these lines do is create the profile formulas only, but doesn't set any values.

def rv(M):
    rv = (4*np.pi*200*rho_crit/(3*M))**(-1./3)
    return rv

def rs(rv, c):
    rs = rv/c
    return rs

def rhos(M, c, rs):
    return M/(4*np.pi*(rs**3)*(np.log(1+c) - c/(1+c)))

def csNFW(cNFW):
    return 0.76*cNFW + 1.36

def a(rv, c):
    rh = rv*(0.6082 - 0.1843*np.log10(c) - 0.1011*np.log10(c)**2 + 0.03918*np.log(c)**3)
    #Lokas+ 2001
    return rh/7.3#5.478 #Lilley+ 2015a(M1, c1)

def main(M1, c1, fgas = 0.17, npartfrac = 1, alpha1 = 1.0):
	npart1 = int(20000000*npartfrac)

	# dm1 = sNFW_density_profile(M="M_1", cs="c_s", rv="r_v") # for the primary cluster
	dm1 = tNFW_density_profile(rho_s="rho_s", r_s="r_s", r_t="r_t")
	"""
	Now we'll decide on some parameters for the DM profiles. We use the yt units package
	to set the actual units. We should always use units of Msun and kpc, since that's what
	cluster_generator understands. I've chosen more or less arbitrary values for now which
	are still physically sensible, but they'll need to be adjusted by you to try to match
	the conditions for A2146.
	"""
	rvir = rv((1-fgas)*M1) #the dark matter mass
	rvir1 = rvir.value * u.kpc
	cs = csNFW(c1)
	rs1 = rs(rvir, c1)
	# a1 = a(rs1, c1)
	rho0 = rhos((1-fgas)*M1, c1, rs1).to('g/cm**3')

	M1 *= u.Msun
	# Now we tell the profiles what the parameter values are
	# dm1.set_param_values(M_1=(1-fgas)*M1, c_s=cs, r_v=rvir)
	dm1.set_param_values(rho_s=rho0, r_s=rs1, r_t=2*rvir1)

	# Now we do all the same stuff for the gas profiles.
	gas1 = vikhlinin_density_profile(rho_0="rho_c1", r_c="r_c1", r_s="r_s1",
	                                 alpha='alpha1', beta="beta1", epsilon='eps1')

	"""
	We won't set the central densities yet, as you'll see why later. Again, picking
	more or less arbitrary values for these parameters...
	"""
	rhoc1=1*u.Msun/(u.kpc**3)
	rc1 = 0.1*rs1
	beta1 = 2./3.
	eps1 = 2.0


	# Setting all parameter values except central densities...
	gas1.set_param_values(r_c1=rc1, r_s1=rs1, alpha1=alpha1, beta1=beta1, eps1=eps1, rho_c1=rhoc1)

	"""
	We didn't set the value of the central densities yet. That's because it makes
	more physical sense to set it by a gas mass fraction of the dark matter. There's
	a function for this:
	"""
	Mgas1 = fgas*M1
	rescale_profile_by_mass(gas1, "rho_c1", Mgas1, rvir1)

	# We get the total density profiles by adding these two together
	tot1 = dm1+gas1

	# cluster_generator doesn't take profiles with units, so we give it the unitless versions
	profiles1 = {"density": gas1.unitless, "total_density": tot1.unitless}

	# now we set up HSE objects, with profiles from 0.1 kpc to 20 Mpc. This may take a while.
	hse1 = HydrostaticEquilibrium.from_scratch("dens_tden", 0.1, 30000.0, profiles1)

	# We now write the hydrostatic profiles to HDF5
	# hse1.write_model_to_h5("ebeling_profile_m%0.1e_c%0.1f_a%0.1f.h5" % (M1, c1, alpha1), in_cgs=True, overwrite=True)
	hse1.write_model_to_h5("tNFW_profile_m%0.1e_c%0.1f_a%0.1f.h5" % (M1, c1, alpha1), in_cgs=True, overwrite=True)
	vir1 = VirialEquilibrium.from_hse_model(hse1)

	# Generate the particles. I've chosen a 3.8/1 ratio in particle number to match the mass ratio
	# This will also take a while.
	parts1 = vir1.generate_particles(npart1)

	# Write the particles
	parts1.write_particles_to_h5("tNFW_%0.1e_particles_c%0.1f_a%0.1f.h5" % (M1, c1, alpha1), in_cgs=True, overwrite=True)

	##vary M1, M2, a1 (=r200/2), a2


if __name__=="__main__":
	M1 = np.array([1.1e15, 9.3e14, 8.1e14])
	M2 = np.array([2.6e14, 2.4e14, 2.1e14])
	ratio = (M2/M1).astype(int)
	for ind in range(3):
		main(M1[ind], 4.1)
		main(M2[ind], 5.2, npartfrac = ratio[ind])

		if ind == 1:
			main(M1[ind], 3.1)
			main(M2[ind], 6.3, npartfrac = ratio[ind])
		
		# main(M2[ind], 6.3, npart1 = 20000000*ratio[ind])
	
