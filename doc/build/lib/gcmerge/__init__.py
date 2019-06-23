name = 'gcmerge'
__all__ = ['read', 'compare', 'plot', 'make_profiles', 'islands']

with open("input_param.txt", "r") as f:
	a = f.readlines()
	inputs = {}
	for line in a:
		key = line.split('\t')[0]
		val = line.split('\t')[-1].split('\n')[0]
		inputs[key] = val 
