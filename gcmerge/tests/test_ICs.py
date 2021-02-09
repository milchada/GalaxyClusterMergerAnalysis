import numpy as np
import matplotlib.pylab as plt
from astropy import constants
import h5py
from scipy.integrate import quad
from scipy.interpolate import interp1d
from astropy.cosmology import Planck15
import astropy.units as u
import glob

kB = constants.k_B.to('keV/K').value
du = u.cm.to('kpc')  

def entropy(prof):
	m_p = constants.m_p.to('g').value
	t = prof['temperature'].value*kB
	d = prof['density'].value/m_p
	return t*pow(d, -2./3)


def get_files(type):
	if type=='a1':
		files = glob.glob('*prof*5.0*5.1*.0.h5')
		files.sort()		
		param = r'$\alpha$'
	if type=='a2':
		files = glob.glob('*prof*1.6*5.2*.0.h5')
		files.sort()		
		param = r'$\alpha$'
	if type=='c1':
		files = glob.glob('*prof*5.0*a0.0.h5')
		files.sort()
		param = r'$r_s/a$'
	if type=='c2':
		files = glob.glob('*prof*1.6*a0.0.h5')
		files.sort()
		param = r'$r_s/a$'
	return files, param

def compare_ics(type='a2'):
	files, param = get_files(type)

	fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows = 2, ncols=2)
	for file in files:
		if type in ['a1','a2']:
			a = file.split('_a')[1].split('.h5')[0]
		elif type in ['c1','c2']:
			a = file.split('a2.0_')[1].split('rs.h5')[0]
		f = h5py.File(file)['fields']
		ax1.loglog(f['radius'].value*du, f['density'], label = '%s = %s' % (param,a))
		ax2.loglog(f['radius'].value*du, f['total_density'], label = '%s = %s' % (param,a))
		ax3.loglog(f['radius'].value*du, f['temperature'], label = '%s = %s' % (param,a))
		ax4.loglog(f['radius'].value*du, entropy(f), label = '%s = %s' % (param,a))
	han, lab = ax1.get_legend_handles_labels()
	ax1.legend(han, lab)
	ax3.set_xlabel('R (kpc)')
	ax4.set_xlabel('R (kpc)')
	ax1.set_ylabel(r'$\rho_g$')
	ax2.set_ylabel(r'$\rho_{\rm tot}$')
	ax3.set_ylabel('T (K)')
	ax4.set_ylabel(r'$K$')
	ax1.set_ylim(1e-29,1e-21) 
	ax2.set_ylim(1e-29,1e-21) 
	ax3.set_ylim(4e6,3.1e7)
	ax3.set_yscale('linear')
	for a in ax.flatten():
		a.set_xlim(1,1e3)
		a.set_xscale('log')
	plt.tight_layout(h_pad=-1.0)
	fig.savefig('%s_ics_compare.png' % type)

from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

def Temp_integral(f, type = 'proj'):
	m_p = constants.m_p.to('g').value
	n_cc = f['density'].value/m_p
	t_kev = f['temperature'].value*kB 
	r = f['radius'].value*du
	dr = r[1:] - r[:-1]
	mazzotta = (n_cc**2)*pow(t_kev, -3./4)
	# dr = r[1]-r[0]
	# dV = dA*dr
	if type == 'proj':
		d_r = interp1d(r, mazzotta)
		n_r = interp1d(r, mazzotta*t_kev)
		def tproj(i):
			numerator = quad(n_r, r[i], r.max())
			denominator = quad(d_r, r[i], r.max())
			try:
				return numerator[0]/denominator[0]
			except ZeroDivisionError:
				return 0
		t_proj = Parallel(n_jobs=num_cores)(delayed(tproj)(i) for i in range(len(r)))
		return np.array(t_proj)
	elif type == 'average':
		dA = 4*np.pi*(r**2)
		d_r = interp1d(r, mazzotta*dA)
		n_r = interp1d(r, mazzotta*t_kev*dA)
		numerator = quad(n_r, r.min(), r.max())
		denominator = quad(d_r, r.min(), r.max())
		tave = numerator[0]/denominator[0]
		return tave

def rvir(profile, overdensity=180, z = 0.5):
	rho_crit = Planck15.critical_density(z).value 
	mtot = profile['total_mass'].value #g/cm^3
	r = profile['radius'].value #cm
	V = 4*np.pi*(r**3)/3.
	cum_rho = mtot/V
	rvir = r[np.argmin(abs(cum_rho - (overdensity*rho_crit)))]
	return rvir

def compare_TIC(cc_files, ncc_files, param, type):
	plt.clf()
	m = cm.ScalarMappable(cmap = cm.Set1, norm=colors.Normalize(0,len(files)))  
	for i in range(len(cc_files)):
		file = cc_files[i]
		if 'alpha' in param:
			a = file.split('a2')[1].split('_rs')[0]
		elif 'r_c' in param:
			a = str(float(file.split('rs')[1].split('.h5')[0])/10.)
		f = h5py.File(file)['fields']
		tproj = Temp_integral(f, 'proj')
		tave = Temp_integral(f, 'average')
		r180 = rvir(f)
		plt.plot(f['radius'].value/r180, np.array(tproj)/tave, color=m.to_rgba(i), label = '%s = %s' % (param,a))
		file = ncc_files[i]
		f = h5py.File(file)['fields']
		tproj = Temp_integral(f, 'proj')
		tave = Temp_integral(f, 'average')
		plt.plot(f['radius'].value/r180, np.array(tproj)/tave, color=m.to_rgba(i), linestyle='dotted',label = '%s = %s' % (param,a))
	plt.legend()
	plt.xlabel(r'r/$r_{180}$')
	plt.ylabel(r'T/<T>')
	plt.ylim(0,2.5)
	# plt.xscale('log')
	plt.xlim(0, 0.5)
	plt.tight_layout(h_pad=-1.0)
	plt.savefig('%s_ics_compare.png' % type)

def plot(file, ax, label=None, linestyle='solid', rnorm=None, znorm=None):
		f = h5py.File(file)['fields']
		t = f['temperature'].value*kB
		r = f['radius'].value*du
		d = f['density'].value
		k = entropy(f)
		if rnorm:
			r500 = rvir(f, overdensity=500)
			r /= r500
		if znorm:
			rho_crit = Planck15.critical_density(znorm).value 
			d /= rho_crit
		ax1, ax2, ax3 = ax.flatten()
		ax1.plot(r, t, color=m.to_rgba(i), label=label, linestyle=linestyle)
		ax2.loglog(r, d, color=m.to_rgba(i), linestyle=linestyle)
		ax3.loglog(r, k, color=m.to_rgba(i), linestyle=linestyle)

def fig12(cc_files, ncc_files, rnorm=None, znorm=None):
	fig, ax = plt.subplots(ncols=3, sharex=True)
	fig.set_figwidth(8.5)
	m = cm.ScalarMappable(cmap = cm.tab10, norm=colors.Normalize(0,len(cc_files)))  
	for i in range(len(cc_files)):
		file = cc_files[i]
		try:
			a = str(float(file.split('rs')[1].split('.h5')[0])/10.)
		except ValueError:
			a = ''
		plot(file, ax, label = '%s = %s' % (param,a), rnorm=rnorm, znorm=znorm)
		plot(ncc_files[i], ax, linestyle='dotted',rnorm=rnorm,znorm=znorm)
	plt.legend()
	plt.xlabel(r'R')
	ax1.ylabel(r'T (K)')
	ax2.ylabel(r'$\rho (g/cm^3)$')
	ax3.ylabel(r'K (keV cm$^2$)')
	plt.xscale('log')
	if rnorm:
		plt.xlim(.005, 2)
	else:
		plt.xlim(1, 1e3)
	plt.tight_layout()
	plt.savefig('%s_ics_compare.png' % type)
