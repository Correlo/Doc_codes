import numpy as np 
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from configparser import ConfigParser
# Modules from the developer
from const import *

# Read Energy_terms_1F.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Energy_terms_1F.ini')
Params = params['params']

# Functions
def Alpha(T0):

	term1 = MIH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MIH)) * SIGMA_IN
	term2 = MEH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MEH)) * SIGMA_EN
	alpha = term1 + term2

	return alpha


# Input parameters
H5file    = Params['H5file'   ]
H5file0   = Params['H5file0'  ]
H5dir     = Params['H5dir'    ]
Outname   = Params['Outname'  ] 
nu_A_file = Params['nu_A_file']
n_file  = int(Params['n_file' ])
dt = float(Params['dt']) # s
with_B0x = bool(int(Params['with_B0x']))


# Parameters that do not come from the user
n_dpts = 100

# Read the fields
with h5py.File(H5file, 'r') as H5obj:

	# Ghost shells
	my_ghost = int(H5obj.attrs['my_ghost'])
	# Obtain dz
	dz = H5obj.attrs['dz'][0]
	# Obtain data
	vz = np.array(H5obj['vz_c'  ])[my_ghost:-my_ghost,0,0]
	P  = np.array(H5obj['pe_c'  ])[my_ghost:-my_ghost,0,0]
	T  = np.array(H5obj['temp_c'])[my_ghost:-my_ghost,0,0]
	D  = np.array(H5obj['rho_c' ])[my_ghost:-my_ghost,0,0]
	mz = H5obj.attrs['dims'][2]
	
	if with_B0x:
	
		B    = np.array(H5obj['bx'])[my_ghost:-my_ghost,0,0]
	
	
with h5py.File(H5file0, 'r') as H5obj0:

    # Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	P0 = np.array(H5obj0['pe_c' ])[my_ghost0:-my_ghost0,0,0]
	D0 = np.array(H5obj0['rho_c'])[my_ghost0:-my_ghost0,0,0]
	T0 = MH * P0 / (BK * D0)
	
	if with_B0x:
	
		B0   = np.array(H5obj0['bx'])[my_ghost0:-my_ghost0,0,0]
	
# Total fields
PT = P + P0 
DT = D + D0  
TT = T + T0  
if with_B0x:

	BT = B + B0

# Obtain the time derivative of the kinetic energy
Filenames = sorted(glob.glob(H5dir + '*'))
vz_aux = np.zeros((n_dpts, mz))
P_aux  = np.zeros((n_dpts, mz))
D_aux  = np.zeros((n_dpts, mz))

for i, H5file in enumerate(Filenames[n_file - int(n_dpts/2):n_file + int(n_dpts/2)]):
	with h5py.File(H5file, 'r') as H5obj:
		vz_aux[i]  = np.array(H5obj['vz_c'  ][my_ghost:-my_ghost,0,0])
		P_aux[i]   = np.array(H5obj['pe_c'  ][my_ghost:-my_ghost,0,0])
		D_aux[i]   = np.array(H5obj['rho_c' ][my_ghost:-my_ghost,0,0])

# Obtain the density of internal energy
PTT = P_aux + P0; del(P_aux)
DTT = D_aux + D0; del(D_aux)  

e = PTT / (gamma - 1)  ; del(PTT)
			
# Obtain the time derivative of these quantities
dedt = np.gradient(e, dt, axis = 0)
to_plot = int(n_dpts/2)

# Obtain the Z-axis
Z = np.arange(mz) * dz * 1e-3 + 520 # Km

# Obtain the density of internal energy times velocity divergence
divue = np.gradient(vz * PT / (gamma - 1), dz) 

# Obtain pressure power 
PP = PT * np.gradient(vz, dz)

# Obtain the ambipolar power 
if with_B0x:

	with h5py.File(nu_A_file, 'r') as f:
		my_ghost = f.attrs['my_ghost']
		nu_A = np.array(f['nu_A'][my_ghost:-my_ghost,0,0])
        
	eta_tilde_A = nu_A / B0**2
	dBdz =  np.gradient(BT, dz) 
	ambW = eta_tilde_A / MU0**2 * dBdz**2 * B
	

plt.close()
plt.figure(figsize = (8, 7))
plt.plot(Z, divue / D0, '-g' , label = r'$\frac{\partial u e}{\partial z}$'  )
plt.plot(Z, -PP / D0, '-m' , label = r'$-\frac{P \partial u }{\partial z}$'  )
plt.plot(Z, dedt[to_plot] / D0, '-y' , label = r'$\frac{\partial e}{\partial t}$')
if with_B0x:
	plt.plot(Z, ambW / D0, '-b' , label = r'$W_{amb}$')

plt.legend(frameon = False, fontsize = 14)
plt.xlim(min(Z), max(Z))
plt.xlabel('Height (km)', fontsize = 14)
plt.ylabel(r'I.E. terms per unit of $\rho$ $(Pa$ $m^3$ $s^{-1}$ $kg^{-1} )$', fontsize = 14)

plt.tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
plt.tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
plt.minorticks_on()

plt.savefig('../Doc_data/Figures/' + Outname)
