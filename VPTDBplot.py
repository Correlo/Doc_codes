#####################################################################################
#    Code to plot Velocity, Pressure, Temperature, Density and Magnetic results     #
#####################################################################################

import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

# Obtain the name of hdf5 file
H5file_2f = str(sys.argv[1]) 
H5file_1f = str(sys.argv[2]) 
#H5file_1f_2 = str(sys.argv[2]) 

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
	B    = np.array(H5obj_2f['bx'    ])[my_ghost_2f:-my_ghost_2f,0,0]
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
	N    = H5obj_1f.attrs['dims'][2]
	dz = H5obj_1f.attrs['dz']
# 	# Obtain time
	t = H5obj_1f.attrs['time']
	
# with h5py.File(H5file_1f_2, 'r') as H5obj_1f:
# 
# 	# Ghost shells
# 	my_ghost_1f = int(H5obj_1f.attrs['my_ghost'])
# 	# Obtain data (1f)
# 	vz_2  = np.array(H5obj_1f['vz_c'  ])[my_ghost_1f:-my_ghost_1f,0,0]
# 	Pe_2  = np.array(H5obj_1f['pe_c'  ])[my_ghost_1f:-my_ghost_1f,0,0]
# 	T_2   = np.array(H5obj_1f['temp_c'])[my_ghost_1f:-my_ghost_1f,0,0]
# 	D_2   = np.array(H5obj_1f['rho_c' ])[my_ghost_1f:-my_ghost_1f,0,0]
# 	B1f_2 = np.array(H5obj_1f['bx'    ])[my_ghost_1f:-my_ghost_1f,0,0]

#print(Pe[100], Pc[100], Pn[100])

eqfilename = '../Doc_data/Equilibrium/Strat_B0x_3212.h5'	
with h5py.File(eqfilename, 'r') as H5obj0:

	# Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	# Obtain data (1f)
	Pc0   = np.array(H5obj0['pe_c'  ])[my_ghost0:-my_ghost0,0,0]
	Pn0   = np.array(H5obj0['pe_n'  ])[my_ghost0:-my_ghost0,0,0]
	Dc0   = np.array(H5obj0['rho_c' ])[my_ghost0:-my_ghost0,0,0]
	Dn0   = np.array(H5obj0['rho_n' ])[my_ghost0:-my_ghost0,0,0]


# z-axis
Z  = np.arange(N) * dz * 1e-3 + 520 # Km

# Create a figure with two subplots
plt.close()
fig, ax = plt.subplots(5, 1, sharex = True, figsize = (12, 8))
fig.subplots_adjust(hspace = 0)

ax[0].set_title('Time = %.1f s' % t, fontsize = 15)
ax[0].plot(Z, vz_c, '-r', label = 'charges'); ax[0].plot(Z, vz_n, '-b' , label = 'neutrals');
ax[0].plot(Z, vz  , '-g', label = '1f'     )
#ax[0].plot(Z, vz_2, '-r', label = '1f'     )
ax[1].plot(Z, Pc / Pc0 , '-r'); ax[1].plot(Z, Pn / Pn0 , '-b');
ax[1].plot(Z, Pe / (Pc0 + Pn0), '-g')
#ax[1].plot(Z, Pe_2 / (Pc0 + Pn0), '-r')
ax[2].plot(Z, Tc  , '-r'); ax[2].plot(Z, Tn , '-b');
ax[2].plot(Z, T , '-g')
#ax[2].plot(Z, T_2 , '-r')
ax[3].plot(Z, Dc /Dc0 , '-r'); ax[3].plot(Z, Dn / Dn0, '-b'); 
ax[3].plot(Z, D / (Dn0 + Dc0), '-g')
#ax[3].plot(Z, D_2 / (Dn0 + Dc0), '-r')
ax[4].plot(Z, B   , '-y'); 
ax[4].plot(Z, B1f, '-g')
#ax[4].plot(Z, B1f_2, '-r')

# r'$P_1$ $[Pa]$'
# r'$\rho_1$ $[kg$ $m^{-3}]$'
ax[0].set_ylabel(r'$u$ $[m$ $s^{-1}]$', fontsize = 14); ax[1].set_ylabel(r'$P_1 / P_0$', fontsize = 14)
ax[2].set_ylabel(r'$T_1$ $[K]$', fontsize = 14); ax[3].set_ylabel(r'$\rho_1 / \rho_0$', fontsize = 14)
ax[4].set_ylabel(r'$B_{1,x}$ [T]', fontsize = 14)
ax[4].set_xlabel('Altitude [Km]', fontsize = 14)

ax[0].legend(fontsize = 13, loc = 'upper left')
ax[4].set_xlim(min(Z), max(Z))

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

plt.show()







