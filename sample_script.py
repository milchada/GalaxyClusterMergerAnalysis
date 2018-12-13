from gcmerge import *
from multiprocessing.pool import ThreadPool as Pool

if __name__=='__main__':
	simid=int(inputs['simid'])
	sciencedir = a[simid].split('\n')[0]
	read_gamer.make_fits(sciencedir, 'temperature',weight_field='mazzotta_weighting',proj=False,fits=True)