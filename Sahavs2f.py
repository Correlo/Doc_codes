import numpy as np
import h5py
import matplotlib.pyplot as plt
from const import *

def Saha_mod(T, n):

	fac1 = 1 / n
	fac2 = 2 * np.pi * ME * BK * T / h**2 
	fac3 = np.exp(- H_ION_ENERGY / (BK * T)) 
	
	S = fac1 * fac2**1.5 *fac3
	
	return S
	
def X(S):

	x = - S / 2 + np.sqrt(S**2 + 4*S) / 2
	return x
	
# Read the densities from 2F equilibrium file
h5eqatm_2F = '../Doc_data/Equilibrium/Strat_B0x_3212.h5'
with h5py.File(h5eqatm_2F, 'r') as h5Obj0:

	_, _, mz = h5Obj0.attrs['dims']
	dz  = h5Obj0.attrs['dz']
	my_ghost = h5Obj0.attrs['my_ghost']
	Dc0 = np.array(h5Obj0['rho_c'][my_ghost:-my_ghost,0,0])
	Dn0 = np.array(h5Obj0['rho_n'][my_ghost:-my_ghost,0,0])
	
# Obtain the density and temperature from the 1F equilibrium file
h5eqatm_1F = '../Doc_data/Equilibrium/Strat_B0x_3208_1f.h5'
with h5py.File(h5eqatm_1F, 'r') as h5Obj0:

	_, _, mz = h5Obj0.attrs['dims']
	my_ghost = h5Obj0.attrs['my_ghost']
	D0 = np.array(h5Obj0['rho_c'][my_ghost:-my_ghost,0,0])
	P0 = np.array(h5Obj0['pe_c'][my_ghost:-my_ghost,0,0])
	
T0 = P0 / (D0 * BK) * MH


Z = np.arange(mz) * dz * 1e-3 + 520 #Km

T00  = T0[0]
Dn00 = Dn0[0]
Dc00 = Dc0[0]

Hn = BK * T00 / (MH * G)
Hc = 2 * BK * T00 / (MH * G)

Dn0_2 = Dn00 * np.exp(- Z * 1e3 / Hn)
Dc0_2 = Dc00 * np.exp(- Z * 1e3 / Hc)

n_2 = Dn0_2 /MH + 2 * Dc0_2 / (2 * MH) 


# Obtain the fraction of neutrals of each model
xi_n_2f = Dn0 / (Dn0 + Dc0)

n = Dn0 /MH + 2 * Dc0 / (2 * MH) 
S = Saha_mod(T00, n_2)
#print(S)
x_i_Saha = X(S)
#print(x_i_Saha)
ne = x_i_Saha * n
nc = 2 * ne
mu_c = 0.5
mu_n = 1
rho_c = MH * mu_c * nc
#print(rho_c)
rho_n = MH * mu_n * (n - ne) 
xi_n_Saha = rho_n / (rho_n + rho_c)
#print(xi_n_Saha)

plt.close()
plt.plot(Z, xi_n_2f  , label = '2F'  )
plt.plot(Z, xi_n_Saha, label = 'Saha')
plt.xlabel('Height (km)')
plt.ylabel(r'$\xi_n$')
plt.legend(frameon = False)
plt.show()


	
	


	