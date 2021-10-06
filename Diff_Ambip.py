import numpy as np
import h5py
from const import *
from configparser import ConfigParser

import matplotlib.pyplot as plt

# Read Diff_Amb.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Diff_Ambip.ini')
Params = params['params']

h5eqatm = Params['h5eqatm']
Output  = Params['Output']

## Functions ##

def Alpha(T0):

	term1 = MIH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MIH)) * SIGMA_IN
	term2 = MEH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MEH)) * SIGMA_EN
	alpha = term1 + term2

	return alpha
	
def Cjbb(Dc0, Dn0, alpha):

	chi_n = Dn0 / (Dc0 + Dn0)
	cjbb = chi_n**2 / (alpha * Dc0 * Dn0) 

	return cjbb

# Obtain temperature from equilibrium file
h5Obj0 = h5py.File(h5eqatm, 'r')
mx, my, mz = h5Obj0.attrs['dims']
my_ghost = h5Obj0.attrs['my_ghost']
Pc0 = np.array(h5Obj0['pe_c' ][:,0,0])
Dc0 = np.array(h5Obj0['rho_c'][:,0,0])
Pn0 = np.array(h5Obj0['pe_n' ][:,0,0])
Dn0 = np.array(h5Obj0['rho_n'][:,0,0])
B0  = np.array(h5Obj0['bx'   ][:,0,0])
dz  = h5Obj0.attrs['dz']
h5Obj0.close()
T0 = Pn0 / (Dn0 * BK) * MH

# Obtain alpha
alpha = Alpha(T0)
# Obtain cjbb coefficient
cjbb = Cjbb(Dc0, Dn0, alpha)
# Obtain the ambipolar diffusivity
nu_A = cjbb * B0**2 #/ MU0

plt.close()
plt.semilogy(np.arange(len(nu_A)) * dz * 1e-6 + 0.5, nu_A)
plt.show()

# Save results 
with h5py.File(Output, 'w') as f:

	f.attrs['dx']   = dz
	f.attrs['dy']   = dz
	f.attrs['dz']   = dz
	f.attrs['dims'] = [mx,my,mz]
	f.attrs['my_ghost'] = my_ghost
	dset = f.create_dataset("nu_A", (len(nu_A), 1, 1),  dtype='f8', data=nu_A, chunks=True, shuffle=True)
	
