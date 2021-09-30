#####################################################################################
#    Code to plot Velocity, Pressure, Temperature, Density and Magnetic results     #
#####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from configparser import ConfigParser
# Modules from the developer
from const import *

# Read iniCond.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/VPTDBplot2.ini')
Params = params['params']

# Parameters of the code
eqfilename = Params['eqfilename']
H5file_2f  = Params['H5file_2f']
H5file_1f  = Params['H5file_1f']
figname    = Params['figname']
div        = float(Params['div'])


with h5py.File(H5file_2f, 'r') as H5obj_2f:

	# Ghost shells
	my_ghost_2f = int(H5obj_2f.attrs['my_ghost'])
	# Obtain dz
	dz = H5obj_2f.attrs['dz']
	# Obtain time
	t = H5obj_2f.attrs['time'] 
	# Obtain data (2f)
	vz_c = np.array(H5obj_2f['vz_c'  ])[my_ghost_2f:-my_ghost_2f,0,0]
	vz_n = np.array(H5obj_2f['vz_n'  ])[my_ghost_2f:-my_ghost_2f,0,0]
	Pc   = np.array(H5obj_2f['pe_c'  ])[my_ghost_2f:-my_ghost_2f,0,0]
	Pn   = np.array(H5obj_2f['pe_n'  ])[my_ghost_2f:-my_ghost_2f,0,0]
	Tc   = np.array(H5obj_2f['temp_c'])[my_ghost_2f:-my_ghost_2f,0,0]
	Tn   = np.array(H5obj_2f['temp_n'])[my_ghost_2f:-my_ghost_2f,0,0]
	Dc   = np.array(H5obj_2f['rho_c' ])[my_ghost_2f:-my_ghost_2f,0,0]
	Dn   = np.array(H5obj_2f['rho_n' ])[my_ghost_2f:-my_ghost_2f,0,0]
	B2f  = np.array(H5obj_2f['bx'    ])[my_ghost_2f:-my_ghost_2f,0,0]
	N    = H5obj_2f.attrs['dims'][2]

with h5py.File(H5file_1f, 'r') as H5obj_1f:

	# Ghost shells
	my_ghost_1f = int(H5obj_1f.attrs['my_ghost'])
	# Obtain data (1f)
	vz  = np.array(H5obj_1f['vz_c'  ])[my_ghost_1f:-my_ghost_1f,0,0]
	Pe  = np.array(H5obj_1f['pe_c'  ])[my_ghost_1f:-my_ghost_1f,0,0]
	T   = np.array(H5obj_1f['temp_c'])[my_ghost_1f:-my_ghost_1f,0,0]
	D   = np.array(H5obj_1f['rho_c' ])[my_ghost_1f:-my_ghost_1f,0,0]
	B1f = np.array(H5obj_1f['bx'    ])[my_ghost_1f:-my_ghost_1f,0,0]
	
with h5py.File(eqfilename, 'r') as H5obj0:

	# Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	# Obtain data (1f)
	Pc0   = np.array(H5obj0['pe_c'  ])[my_ghost0:-my_ghost0,0,0]
	Pn0   = np.array(H5obj0['pe_n'  ])[my_ghost0:-my_ghost0,0,0]
	Dc0   = np.array(H5obj0['rho_c' ])[my_ghost0:-my_ghost0,0,0]
	Dn0   = np.array(H5obj0['rho_n' ])[my_ghost0:-my_ghost0,0,0]

# Temporal
vz2f = ((Dc0 + Dc) * vz_c + (Dn0 + Dn) * vz_n) / (Dc0 + Dc + Dn0 + Dn)
P2f  = Pc + Pn
T2f  = ((Dc0 + Dc) * Tc   + (Dn0 + Dn) * Tn  ) / (Dc0 + Dc + Dn0 + Dn)
D2f  = Dc + Dn

# z-axis
Z  = np.arange(N) * dz * 1e-3 + 520 # Km

# Create a figure with two subplots
plt.close()
fig, ax = plt.subplots(5, 1, sharex = True, figsize = (12, 8))
fig.subplots_adjust(hspace = 0)

ax[0].plot(Z, vz , '-k', label = '1f ambipolar'); ax[0].plot(Z, vz2f, '--r' , label = '2f');
ax[1].plot(Z, Pe , '-k'); ax[1].plot(Z, P2f , '--r')
ax[2].plot(Z, T  , '-k'); ax[2].plot(Z, T2f , '--r')
ax[3].plot(Z, D  , '-k'); ax[3].plot(Z, D2f , '--r')
ax[4].plot(Z, B1f, '-k'); ax[4].plot(Z, B2f , '--r')
for i in range(5): ax[i].axvline(Z[int(N/div)], linestyle = '--', color = 'k')

ax[0].set_ylabel(r'$u$ $(m$ $s^{-1})$', fontsize = 14); ax[1].set_ylabel(r'$P_1$ $(Pa)$', fontsize = 14)
ax[2].set_ylabel(r'$T_1$ $(K)$', fontsize = 14); ax[3].set_ylabel(r'$\rho_1$ $(kg$ $m^{-3})$', fontsize = 14)
ax[4].set_ylabel(r'$B_{1,x}$ (T)', fontsize = 14)
ax[4].set_xlabel('Height (km)', fontsize = 14)

ax[0].legend(fontsize = 13, loc = 'upper left')
ax[4].set_xlim(min(Z), max(Z) - 50)

ax[0].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, top = True, bottom = False)
ax[0].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[1].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, bottom = False)
ax[1].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[2].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, bottom = False)
ax[2].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[3].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, bottom = False)
ax[3].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[4].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True)
ax[4].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)


for i in range(5): ax[i].ticklabel_format(axis = "y", style = "sci", scilimits = (0,0))
for i in range(5): ax[i].minorticks_on()

plt.savefig(figname)







