import cluster_generator as cg
import numpy as np
import glob

with open('Input__TestProb') as f:
    a = f.readlines()

for line in a:
    if 'Merger_Coll_PosX1' in line:
        x1 = float(line.split('PosX1 ')[1].split('\t')[0])
    if 'Merger_Coll_PosY1' in line:
        y1 = float(line.split('PosY1 ')[1].split('\t')[0])
    if 'Merger_Coll_PosX2' in line:
        x2 = float(line.split('PosX2 ')[1].split('\t')[0])
    if 'Merger_Coll_PosY2' in line:
        y2 = float(line.split('PosY2 ')[1].split('\t')[0])
    if 'Merger_File_Prof1' in line:
        profs.append(line.split('       ')[1].split(' \t')[0])
    if 'Merger_File_Prof2' in line:
        profs.append(line.split('       ')[1].split('\t')[0])

p1 = cg.ClusterModel.from_h5_file(profs[0])
p2 = cg.ClusterModel.from_h5_file(profs[1])

beta = 100.0

p1.set_magnetic_field_from_beta(beta, gaussian=False)
p2.set_magnetic_field_from_beta(beta, gaussian=False)


buffer = 200.0
ctr1 = np.array([x1,y1,7000.0])
ctr2 = np.array([x2,y2,7000.0])

with open('Input__Parameter') as f: a = f.readlines()
for line in a:
    if 'BOX_SIZE' in line:
        size = 1000*float(line.split('BOX_SIZE                      ')[1].split(' ')[0])
        #Mpc to kpc

left_edge = np.array([0.0]*3)-buffer
right_edge = np.array([size]*3)+buffer

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
