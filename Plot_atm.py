import numpy as np
import matplotlib.pyplot as plt
import h5py
from configparser import ConfigParser
# Modules from the developer
from const import *

# Read Plot_atm.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Plot_atm.ini')
Params = params['params']

# Read equilibrium files 
H5file0_2f = Params['H5file0_2f'] 
H5file0_1f = Params['H5file0_1f']
outname    = Params['outname']


with h5py.File(H5file0_2f, 'r') as H5Obj0:
	dz = H5Obj0.attrs['dz'] 
	my_ghost = H5Obj0.attrs['my_ghost']
	mz  = H5Obj0.attrs['dims'][-1]
	Pn0 = np.array(H5Obj0['pe_n' ][my_ghost:-my_ghost,0,0])
	Pc0 = np.array(H5Obj0['pe_c' ][my_ghost:-my_ghost,0,0])
	Dn0 = np.array(H5Obj0['rho_n'][my_ghost:-my_ghost,0,0])
	Dc0 = np.array(H5Obj0['rho_c'][my_ghost:-my_ghost,0,0])
	
P0 = Pn0 + Pc0
D0 = Dn0 + Dc0

with h5py.File(H5file0_1f, 'r') as H5Obj0:
	my_ghost = H5Obj0.attrs['my_ghost']
	mu_c0 = np.array(H5Obj0['mu_c'][my_ghost:-my_ghost,0,0])
	
# Temperatures
T0  = MH * mu_c0 * P0  / (BK * D0)
Tn0 = MH * Pn0 / (BK * Dn0)
Tc0 = MH * Pc0 / (2 * BK * Dc0)

# Sound speed
cn0 = np.sqrt(gamma * Pn0 / Dn0)
cc0 = np.sqrt(gamma * Pc0 / Dc0)
cs0 = np.sqrt(gamma * P0  / D0 )

# Output figures to characterize solar chromosphere
Z = np.arange(mz) * dz * 1e-3 + 520 # Km

plt.close()
fig, ax = plt.subplots(2, 2, figsize = (15, 13))

ax[0,0].semilogy(Z, Pc0, '-r', label = 'charges', linewidth = 3)
ax[0,0].semilogy(Z, Pn0, '-b', label = 'neutrals', linewidth = 3)
ax[0,0].semilogy(Z, P0 , '--g', label = '1f', linewidth = 3)
ax[0,0].set_xlim(min(Z), max(Z))
ax[0,0].set_xlabel('Height (km)', fontsize = 15)
ax[0,0].set_ylabel('P (Pa)', fontsize = 15)
ax[0,0].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=14, top = True, right = True)
ax[0,0].tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=14, top = True, right = True)
ax[0,0].legend(frameon = False, fontsize = 14)
ax[0,0].minorticks_on()

ax[0,1].semilogy(Z, Dc0, '-r', linewidth = 3)
ax[0,1].semilogy(Z, Dn0, '-b', linewidth = 3)
ax[0,1].semilogy(Z, D0 , '--g', linewidth = 3)
ax[0,1].set_xlim(min(Z), max(Z))
ax[0,1].set_xlabel('Height (km)', fontsize = 15)
ax[0,1].set_ylabel(r'$\rho$ ($kg/m^3$)', fontsize = 15)
ax[0,1].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=14, top = True, right = True)
ax[0,1].tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=14, top = True, right = True)
ax[0,1].minorticks_on()

ax[1,0].plot(Z, Tn0, '-b', linewidth = 3)
ax[1,0].plot(Z, T0, '--g', linewidth = 3)
ax[1,0].set_xlim(min(Z), max(Z))
ax[1,0].set_xlabel('Height (km)', fontsize = 15)
ax[1,0].set_ylabel('T (K)', fontsize = 15)
ax[1,0].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=14, top = True, right = True)
ax[1,0].tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=14, top = True, right = True)
ax[1,0].minorticks_on()

ax[1,1].plot(Z, cc0 * 1e-3, '-r', linewidth = 3)
ax[1,1].plot(Z, cn0 * 1e-3, '-b', linewidth = 3)
ax[1,1].plot(Z, cs0 * 1e-3, '--g', linewidth = 3)
ax[1,1].set_xlim(min(Z), max(Z))
ax[1,1].set_xlabel('Height (km)', fontsize = 15)
ax[1,1].set_ylabel(r'$c_s$ (km/s)', fontsize = 15)
ax[1,1].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=14, top = True, right = True)
ax[1,1].tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=14, top = True, right = True)
ax[1,1].minorticks_on()



plt.savefig(outname)


