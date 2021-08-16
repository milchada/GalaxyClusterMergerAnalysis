import cluster_generator as cg
import numpy as np
import glob

profs = ['sNFW_profile_m1.6e+14_c5.0_a2.0rs0.6rvir_rc0.1rvir_beta0.67_eps2.h5',
	 'sNFW_profile_m5.0e+14_c5.0_a2.0rs0.6rvir_rc0.1rvir_beta0.67_eps2.h5']
p1 = cg.ClusterModel.from_h5_file(profs[0])
p2 = cg.ClusterModel.from_h5_file(profs[1])

beta = 200.0

p1.set_magnetic_field_from_beta(beta, gaussian=False)
p2.set_magnetic_field_from_beta(beta, gaussian=False)


buffer = 200.0

left_edge = np.array([0.0]*3)-buffer
right_edge = np.array([14.0e3]*3)+buffer
ctr1 = np.array([8552.0,7000.0,7000.0])
ctr2 = np.array([5552.26,7100.0,7000.0])

dims = (256,)*3

bfield = cg.RadialRandomMagneticVectorPotential(left_edge, right_edge, dims, 
                                                10.0, 500.0, ctr1, p1, 
                                                ctr2=ctr2, profile2=p2)

bfield.write_to_h5("B_IC", overwrite=True, length_unit="Mpc",
                   field_unit="sqrt(1.0e14*Msun/Mpc**3)*Mpc/(10*Gyr)")

import h5py
def scale(filename, factor, Bfile):
    b = h5py.open(Bfile)
    bnew = h5py.open(filename, mode='a')
    bnew.create_dataset('x', data=b['x'])
    bnew.create_dataset('y', data=b['x'])
    bnew.create_dataset('z', data=b['x'])
    bnew.create_dataset('magnetic_vector_potential_x', data=b['magnetic_vector_potential_x']*factor)
    bnew.create_dataset('magnetic_vector_potential_y', data=b['magnetic_vector_potential_y']*factor)
    bnew.create_dataset('magnetic_vector_potential_z', data=b['magnetic_vector_potential_z']*factor)
    bnew.flush()
