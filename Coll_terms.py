import numpy as np 
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from configparser import ConfigParser
# Modules from the developer
from const import *


# Functions
def Alpha(T0):

	term1 = MIH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MIH)) * SIGMA_IN
	term2 = MEH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MEH)) * SIGMA_EN
	alpha = term1 + term2

	return alpha

# Read Coll_terms.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Coll_terms.ini')
Params = params['params']

# Input parameters
H5file_1f  = Params['H5file_1f' ]
H5file_2f  = Params['H5file_2f' ]
H5file0_1f = Params['H5file0_1f']
H5file0_2f = Params['H5file0_2f']
Outname    = Params['Outname'   ] 
nu_A_file  = Params['nu_A_file' ]
n_file  = int(Params['n_file' ])
dt = float(Params['dt']) # s


### Single fluid collisional terms ###

# Read the fields
with h5py.File(H5file_1f, 'r') as H5obj:

	# Ghost shells
	my_ghost = int(H5obj.attrs['my_ghost'])
	# Obtain dz
	dz = H5obj.attrs['dz'][0]
	# Obtain data
	B = np.array(H5obj['bx'])[my_ghost:-my_ghost,0,0]
	
	
with h5py.File(H5file0_1f, 'r') as H5obj0:

    # Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	D0 = np.array(H5obj0['rho_c'])[my_ghost0:-my_ghost0,0,0]
	B0 = np.array(H5obj0['bx'])[my_ghost0:-my_ghost0,0,0]
		
		
# Obtain the ambipolar power 
with h5py.File(nu_A_file, 'r') as f:
	my_ghost = f.attrs['my_ghost']
	nu_A = np.array(f['nu_A'][my_ghost:-my_ghost,0,0])
        
# Obtain the total fields
BT = B + B0
        
eta_tilde_A = nu_A / B0**2
dBdz =  np.gradient(BT, dz) 
ambW = eta_tilde_A / MU0**2 * dBdz**2 * B
	

### Two fluid collisional terms ###

# Read the fields
with h5py.File(H5file_2f, 'r') as H5obj:

	# Ghost shells
	my_ghost = int(H5obj.attrs['my_ghost'])
	# Obtain data
	vz_c = np.array(H5obj['vz_c'  ])[my_ghost:-my_ghost,0,0]
	vz_n = np.array(H5obj['vz_n'  ])[my_ghost:-my_ghost,0,0]
	Tc   = np.array(H5obj['temp_c'])[my_ghost:-my_ghost,0,0]
	Tn   = np.array(H5obj['temp_n'])[my_ghost:-my_ghost,0,0]
	mz   = H5obj.attrs['dims'][2]
	
with h5py.File(H5file0_2f, 'r') as H5obj0:

    # Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	Dc0 = np.array(H5obj0['rho_c'])[my_ghost0:-my_ghost0,0,0]
	Dn0 = np.array(H5obj0['rho_n'])[my_ghost0:-my_ghost0,0,0]
	Pc0 = np.array(H5obj0['pe_c' ])[my_ghost0:-my_ghost0,0,0]
	Pn0 = np.array(H5obj0['pe_n' ])[my_ghost0:-my_ghost0,0,0]
	Tc0  = MH * Pc0 / (2 * BK * Dc0)
	Tn0  = MH * Pn0 / (BK * Dn0)

# Obtain the total fields	
TcT = Tc + Tc0  
TnT = Tn + Tn0 

# Obtain alpha coefficient
alpha = Alpha((Tc0 + Tn0) / 2) #Alpha((TcT + TnT) / 2)
# Obtain collisional pressure per unit of time
Pcoll = 0.5 * alpha * Dc0 * Dn0 * (vz_c - vz_n)**2
# Obtain thermal exchange
TE = 1 / (gamma - 1) * BK/MH * (TcT - TnT) * alpha * Dc0 * Dn0 

# Obtain the Z-axis
Z = np.arange(mz) * dz * 1e-3 + 520 # Km

plt.close()
plt.figure(figsize = (8, 7))
plt.plot(Z, TE / Dc0, '-b' , label = r'$TE$')
plt.plot(Z, Pcoll / Dc0, '-r' , label = r'$FH$')
plt.plot(Z, ambW * 1e2 / D0, '-g' , label = r'$W_{amb}$')
plt.legend(frameon = False, fontsize = 14)
plt.xlim(min(Z), max(Z))
plt.xlabel('Height (km)', fontsize = 14)
plt.ylabel(r'I.E. terms per unit of $\rho$ $(Pa$ $m^3$ $s^{-1}$ $kg^{-1} )$', fontsize = 14)

plt.tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
plt.tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
plt.minorticks_on()
plt.tight_layout(pad = 3.0)

plt.savefig('../Doc_data/Figures/' + Outname)


	
	

	