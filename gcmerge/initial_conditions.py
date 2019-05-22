#masses (Msun) from King et al 2016
M1 = 1.06e15
M2 = 2.6e14

#Diemer & Kravtsov 2016
#scatter in log-c = 0.16 dex
#mean relation eq 9, 10; constants in para below
def c(M):
	f0 = 7.14
	f1 = 1.6
	n0 = 4.1
	n1 = .75
	alpha = 1.4
	beta = .67
	n = -2.3
	cmin = f0 + f1*n 
	vmin = n0 + n1*n 
	#mass in Msun --> v?
	meanc = cmin/2 * ((v/vmin)**(-alpha) + (v/vmin)**beta)
	return 10**(meanc-0.16), 10**meanc, 10**(meanc+0.16)

#likelihood of merger with given mass ratio
#from Fakhouri & Ma 2009 using Millenium+FOF merger trees

#analytical arguments from Sarazin 2001
def v(M1, M2, di):
	#relative velocity in km/s at initial separation di in Mpc
	#M1, M2 in 10^15 solar masses
	return 2930*np.sqrt(M1+M2)*pow(di, -0.5)

def f(M1, M2):
	return (M1 + M2)**3 * (1 - (M1**(5./3) + M2**(5./3))/(M1+M2)**5./3)**(3./2)

def vr(M1, M2, lmd, di):
	#returns tangential velocity in km/s given spin parameter lmd, 
	#halo masses in 10^15 Msun and initial separation in Mpc
	return 93*(lmd/0.05)*np.sqrt(M1+M2)*pow(di, -0.5)*f(M1, M2)/2

def b(M1, M2, lmb, di):
	#characteristic impact parameter in kpc
	#for given halo masses in 10^15 Msun, spin parameter lmb, and initial separation in Mpc
	return 160 * (lmd/0.05) * np.sqrt(d) * f(M1, M2)/2.

def vikhlinin_profile(r, r_c, r_s, alpha, beta, epsilon, gamma=3):
	"""
	rc = core radius
	rs = scale radius
	alpha = inner slope
	beta = intermediate slope
	epsilon = width of transition to outer slope
	gamma = outer slope
	"""
	return (r/r_c)**(-alpha/2.) * \
        (1+(r/r_c)**2)**(-3*beta/2.+alpha/4.) * \
        (1+(r/r_s)**gamma)**(-epsilon/gamma/2.)

