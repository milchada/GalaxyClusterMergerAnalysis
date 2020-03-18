name = 'gcmerge'
__all__ = ['align_bcgs', 
	'bcg_dist', 
	'bcg_velocities', 
	'compare',
	'concentration',
	'find_features',
	'fit_arcs',
	'islands',
	'make_islands',
	'make_movie',
	'make_profiles',
	'measure_feature',
	'points_above_gradient',
	'read',
	'read_sim']

with open("input_param.txt", "r") as f:
	a = f.readlines()
	inputs = {}
	for line in a:
		key = line.split('\t')[0]
		val = line.split('\t')[-1].split('\n')[0]
		inputs[key] = val 
