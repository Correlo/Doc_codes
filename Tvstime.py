import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob
from configparser import ConfigParser

# function to obtain minimums and maximus
def Peaks(Z, Y):

	# Initialize the lists
	up, down  = [], []
	
	for i in range(1, len(Y) - 1):
		
		sb = np.sign(Y[i] - Y[i - 1])
		sf = np.sign(Y[i] - Y[i + 1])
		if (sb == sf):
			
			if sb ==   1: up.append([Z[i], Y[i]])
			if sb == - 1: down.append([Z[i], Y[i]])
			
	up, down = np.array(up), np.array(down)
	
	return up, down
	
def T_mean_2f_func(snap_i, snap_f, Filenames):

	# Read density and temperature perturbations
	T_mean = np.zeros_like(Z)
	for i in range(snap_i, snap_f):
		with h5py.File(Filenames[i], 'r') as H5obj:
	
			# Ghost shells
			if (i == snap_i): my_ghost = int(H5obj.attrs['my_ghost'])
			Dn = np.array(H5obj['rho_n' ][my_ghost:-my_ghost,0,0])
			Dc = np.array(H5obj['rho_c' ][my_ghost:-my_ghost,0,0])
			Tn = np.array(H5obj['temp_n'][my_ghost:-my_ghost,0,0])
			Tc = np.array(H5obj['temp_c'][my_ghost:-my_ghost,0,0])
		
	
		# Obtain the temperature
		T2f  = ((Dc0 + Dc) * Tc   + (Dn0 + Dn) * Tn  ) / (Dc0 + Dc + Dn0 + Dn)
		up, down = Peaks(Z, T2f)
		up_interp   = np.interp(Z, up[:,0]  , up[:,1]  )
		down_interp = np.interp(Z, down[:,0], down[:,1])
		T = (up_interp + down_interp) / 2
	
		T_mean = T_mean + T

	T_mean = T_mean / 50
	
	return T_mean
	
def T_mean_1f_func(snap_i, snap_f, Filenames):

	# Read density and temperature perturbations
	T_mean = np.zeros_like(Z)
	for i in range(snap_i, snap_f):
		with h5py.File(Filenames[i], 'r') as H5obj:
	
			# Ghost shells
			if (i == snap_i): my_ghost = int(H5obj.attrs['my_ghost'])
			Dc = np.array(H5obj['rho_c' ][my_ghost:-my_ghost,0,0])
			Tc = np.array(H5obj['temp_c'][my_ghost:-my_ghost,0,0])
	
		up, down = Peaks(Z, Tc)
		up_interp   = np.interp(Z, up[:,0]  , up[:,1]  )
		down_interp = np.interp(Z, down[:,0], down[:,1])
		T = (up_interp + down_interp) / 2
	
		T_mean = T_mean + T

	T_mean = T_mean / 50
	
	return T_mean

# Read TvsZ.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Tvstime.ini')
Params = params['params']

# Obtain the filenames
Filenames_2f = sorted(glob.glob(Params['H5dir_2f'] + '*.h5'))
Filenames_1f = sorted(glob.glob(Params['H5dir_1f'] + '*.h5'))
eqfilename   = Params['eqfilename']
figname      = Params['figname']
div          = float(Params['div'])

# Read the equilibrium files
with h5py.File(eqfilename, 'r') as H5obj0:

	# Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	# Obtain dz
	dz  = H5obj0.attrs['dz']
	Dn0 = np.array(H5obj0['rho_n'][my_ghost0:-my_ghost0,0,0])
	Dc0 = np.array(H5obj0['rho_c'][my_ghost0:-my_ghost0,0,0])
	mz  = H5obj0.attrs['dims'][2]
	
Z    = np.arange(mz) * dz * 1e-3 + 520 # Km
	
# Two-fluid temperature perturbation 
T_snap_2f = 50
snap_2f   = np.arange(1300, 2500, T_snap_2f)
T_mean_t_2f = np.zeros(len(snap_2f))
for j, snap_i in enumerate(snap_2f):

	T_mean = T_mean_2f_func(snap_i, snap_i + T_snap_2f, Filenames_2f)
	T_mean_t_2f[j] = T_mean[-140]	# 760
	
# Single-fluid temperature perturbation 
T_snap_1f = 50
snap_1f   = np.arange(1300, 2500, T_snap_1f)
T_mean_t_1f = np.zeros(len(snap_1f))
for j, snap_i in enumerate(snap_1f):

	T_mean = T_mean_1f_func(snap_i, snap_i + T_snap_1f, Filenames_1f)
	T_mean_t_1f[j] = T_mean[-140]
	

plt.close()
plt.figure(figsize = (8,6))
plt.plot(snap_2f/10, T_mean_t_2f, '.k', label = '2f')
plt.plot(snap_1f/10   , T_mean_t_1f, 'dk', label = '1f')
plt.xlabel('t (s)', fontsize = 14)
plt.ylabel(r'$T_1$ (K)', fontsize = 14)
plt.xlim(min(snap_1f/10) - 2, max(snap_1f/10) + 2)
plt.legend(frameon = False, fontsize = 13)
plt.text(152.5, max(max(T_mean_t_2f), max(T_mean_t_1f)) - 0.05,
         r'z ~ 1950 km', fontsize = 13)
plt.minorticks_on()
plt.tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
plt.tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
plt.savefig(figname)
	

