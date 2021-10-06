import numpy as np 
import matplotlib.pyplot as plt
import h5py
import glob
import sys
from configparser import ConfigParser
# Modules from the developer
from const import *

# Read Energy_terms.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Energy_terms.ini')
Params = params['params']

# Functions
def Alpha(T0):

	term1 = MIH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MIH)) * SIGMA_IN
	term2 = MEH / MH**2 * np.sqrt((8 * BK * T0) / (np.pi * MEH)) * SIGMA_EN
	alpha = term1 + term2

	return alpha


# Input parameters
H5file  = Params['H5file' ]
H5file0 = Params['H5file0']
H5dir   = Params['H5dir'  ]
Outname = Params['Outname'] 
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
	vz_c = np.array(H5obj['vz_c'  ])[my_ghost:-my_ghost,0,0]
	vz_n = np.array(H5obj['vz_n'  ])[my_ghost:-my_ghost,0,0]
	Pc   = np.array(H5obj['pe_c'  ])[my_ghost:-my_ghost,0,0]
	Pn   = np.array(H5obj['pe_n'  ])[my_ghost:-my_ghost,0,0]
	Tc   = np.array(H5obj['temp_c'])[my_ghost:-my_ghost,0,0]
	Tn   = np.array(H5obj['temp_n'])[my_ghost:-my_ghost,0,0]
	Dc   = np.array(H5obj['rho_c' ])[my_ghost:-my_ghost,0,0]
	Dn   = np.array(H5obj['rho_n' ])[my_ghost:-my_ghost,0,0]
	mz   = H5obj.attrs['dims'][2]
	
	if with_B0x:
	
		B    = np.array(H5obj['bx'])[my_ghost:-my_ghost,0,0]
	
	
with h5py.File(H5file0, 'r') as H5obj0:

    # Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	Pc0 = np.array(H5obj0['pe_c' ])[my_ghost0:-my_ghost0,0,0]
	Pn0 = np.array(H5obj0['pe_n' ])[my_ghost0:-my_ghost0,0,0]
	Dc0 = np.array(H5obj0['rho_c'])[my_ghost0:-my_ghost0,0,0]
	Dn0 = np.array(H5obj0['rho_n'])[my_ghost0:-my_ghost0,0,0]
	Tc0  = MH * Pc0 / (2 * BK * Dc0)
	Tn0  = MH * Pn0 / (BK * Dn0)
	
	if with_B0x:
	
		B0   = np.array(H5obj0['bx'])[my_ghost0:-my_ghost0,0,0]
	
# Total fields
PcT = Pc + Pc0 
PnT = Pn + Pn0
DcT = Dc + Dc0  
DnT = Dn + Dn0 
TcT = Tc + Tc0  
TnT = Tn + Tn0 
if with_B0x:

	BT = B + B0

# Obtain the time derivative of the kinetic energy
Filenames = sorted(glob.glob(H5dir + '*'))
vz_c_aux = np.zeros((n_dpts, mz))
vz_n_aux = np.zeros((n_dpts, mz))
Pc_aux   = np.zeros((n_dpts, mz))
Pn_aux   = np.zeros((n_dpts, mz))
Dc_aux   = np.zeros((n_dpts, mz))
Dn_aux   = np.zeros((n_dpts, mz))

for i, H5file in enumerate(Filenames[n_file - int(n_dpts/2):n_file +int(n_dpts/2)]):
	with h5py.File(H5file, 'r') as H5obj:
		vz_c_aux[i] = np.array(H5obj['vz_c'  ][my_ghost:-my_ghost,0,0])
		vz_n_aux[i] = np.array(H5obj['vz_n'  ][my_ghost:-my_ghost,0,0])
		Pc_aux[i]   = np.array(H5obj['pe_c'  ][my_ghost:-my_ghost,0,0])
		Pn_aux[i]   = np.array(H5obj['pe_n'  ][my_ghost:-my_ghost,0,0])
		Dc_aux[i]   = np.array(H5obj['rho_c' ][my_ghost:-my_ghost,0,0])
		Dn_aux[i]   = np.array(H5obj['rho_n' ][my_ghost:-my_ghost,0,0])

# Obtain the density of kinetic energy and internal energy
PcTT = Pc_aux + Pc0; del(Pc_aux)
PnTT = Pn_aux + Pn0; del(Pn_aux)
DcTT = Dc_aux + Dc0; del(Dc_aux)  
DnTT = Dn_aux + Dn0; del(Dn_aux)

ec = PcTT / (gamma - 1)  ; del(PcTT)
en = PnTT / (gamma - 1)	 ; del(PnTT)	
			
# Obtain the time derivative of these quantities
decdt = np.gradient(ec, dt, axis = 0)
dendt = np.gradient(en, dt, axis = 0)
to_plot =int(n_dpts/2)

# Obtain the Z-axis
Z = np.arange(mz) * dz * 1e-3 + 520 # Km

plt.close()
plt.plot(np.arange(len(Pc0))*dz*1e-3 + 520, Pc0 + Pn0 / (gamma - 1))
plt.show()

# Obtain the density of internal energy times velocity divergence
divuec = np.gradient(vz_c * PcT / (gamma - 1), dz) 
divuen = np.gradient(vz_n * PnT / (gamma - 1), dz) 

# Obtain pressure power 
PPc = PcT * np.gradient(vz_c, dz)
PPn = PnT * np.gradient(vz_n, dz)

# Obtain alpha coefficient
alpha = Alpha((Tc0 + Tn0) / 2) #Alpha((TcT + TnT) / 2)
# Obtain collisional pressure per unit of time
Pcoll = 0.5 * alpha * Dc0 * Dn0 * (vz_c - vz_n)**2
# Obtain thermal exchange
TE = 1 / (gamma - 1) * BK/MH * (TcT - TnT) * alpha * Dc0 * Dn0 

plt.close()
fig, ax = plt.subplots(1, 2, figsize = (15, 7))
ax[0].plot(Z, divuen / Dn0, '-g' , label = r'$\frac{\partial u_n e_{n}}{\partial z}$'  )
ax[0].plot(Z, -PPn / Dn0, '-m' , label = r'$-\frac{P_{n} \partial u_n }{\partial z}$'  )
ax[0].plot(Z, dendt[to_plot] / Dn0 , '-y' , label = r'$\frac{\partial e_{n}}{\partial t}$'  )
ax[0].plot(Z, - TE / Dn0, '-b' , label = r'- $TE$')
ax[0].plot(Z, Pcoll / Dn0, '-r' , label = r'$FH$')

ax[1].plot(Z, divuec / Dc0, '-g' , label = r'$\frac{\partial u_c e_{c}}{\partial z}$'  )
ax[1].plot(Z, -PPc / Dc0, '-m' , label = r'$-\frac{P_{c} \partial u_c }{\partial z}$'  )
ax[1].plot(Z, decdt[to_plot]  / Dc0 , '-y' , label = r'$\frac{\partial e_{c}}{\partial t}$'  )
ax[1].plot(Z, TE / Dc0, '-b' , label = r'$TE$')
ax[1].plot(Z, Pcoll / Dc0, '-r' , label = r'$FH$')

ax[0].legend(frameon = False, fontsize = 14)
ax[1].legend(frameon = False, fontsize = 14)
ax[0].set_xlim(min(Z), max(Z))
ax[1].set_xlim(min(Z), max(Z))
ax[0].set_xlabel('Height (km)', fontsize = 14)
ax[1].set_xlabel('Height (km)', fontsize = 14)
ax[0].set_ylabel(r'I.E. terms for neutrals per unit of $\rho_n$ $(Pa$ $m^3$ $s^{-1}$ $kg^{-1} )$', fontsize = 14)
ax[1].set_ylabel(r'I.E. terms for charges per unit of $\rho_c$ $(Pa$ $m^3$ $s^{-1}$ $kg^{-1} )$', fontsize = 14)

ax[0].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
ax[0].tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
ax[1].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
ax[1].tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
ax[0].minorticks_on(); ax[1].minorticks_on()
plt.tight_layout(pad = 3.0)

plt.savefig('../Doc_data/Figures/' + Outname)
