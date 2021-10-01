import numpy as np
import matplotlib.pyplot as plt
import h5py
import glob

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

def T_mean_2f_func(snap_i, snap_f, my_ghost, Filenames):

	# Read density and temperature perturbations
	T_mean = np.zeros_like(Z)
	for i in range(snap_i, snap_f):
		with h5py.File(Filenames[i], 'r') as H5Obj:
	
			Dn = np.array(H5Obj['rho_n' ][my_ghost:-my_ghost,0,0])
			Dc = np.array(H5Obj['rho_c' ][my_ghost:-my_ghost,0,0])
			Tn = np.array(H5Obj['temp_n'][my_ghost:-my_ghost,0,0])
			Tc = np.array(H5Obj['temp_c'][my_ghost:-my_ghost,0,0])
	
		# Obtain the temperature
		T2f  = ((Dc0 + Dc) * Tc   + (Dn0 + Dn) * Tn  ) / (Dc0 + Dc + Dn0 + Dn)
		up, down = Peaks(Z, T2f)
		up_interp   = np.interp(Z, up[:,0]  , up[:,1]  )
		down_interp = np.interp(Z, down[:,0], down[:,1])
		T = (up_interp + down_interp) / 2
	
		T_mean = T_mean + T

	T_mean = T_mean / 50
	
	return T_mean
	
def T_mean_1f_func(snap_i, snap_f, my_ghost, Filenames):

	# Read density and temperature perturbations
	T_mean = np.zeros_like(Z)
	for i in range(snap_i, snap_f):
		with h5py.File(Filenames[i], 'r') as H5Obj:
	
			Dc = np.array(H5Obj['rho_c' ][my_ghost:-my_ghost,0,0])
			Tc = np.array(H5Obj['temp_c'][my_ghost:-my_ghost,0,0])
	
		up, down = Peaks(Z, Tc)
		up_interp   = np.interp(Z, up[:,0]  , up[:,1]  )
		down_interp = np.interp(Z, down[:,0], down[:,1])
		T = (up_interp + down_interp) / 2
	
		T_mean = T_mean + T

	T_mean = T_mean / 50
	
	return T_mean


# Obtain the filenames
Filenames_2f = sorted(glob.glob('../MA_OUT_strat_P5s_ststep/*.h5'   ))
Filenames_1f = sorted(glob.glob('../MA_OUT_strat_P5s_amb_2/*.h5'))

# The number of ghost shells
my_ghost = 3

# Read the equilibrium files
H5file0 = '../Equilibrium/Strat_B0x_3206.h5'
with h5py.File(H5file0, 'r') as H5Obj0:

	dz  = H5Obj0.attrs['dz']
	Dn0 = np.array(H5Obj0['rho_n'][my_ghost:-my_ghost,0,0])
	Dc0 = np.array(H5Obj0['rho_c'][my_ghost:-my_ghost,0,0])
	
Z    = np.arange(len(Dc0)) * dz * 1e-3 + 500 # Km
	
# Two-fluid temperature perturbation 
T_snap_2f = 50
snap_2f   = np.arange(1300, 2000, T_snap_2f)
T_mean_t_2f = np.zeros(len(snap_2f))
for j, snap_i in enumerate(snap_2f):

	T_mean = T_mean_2f_func(snap_i, snap_i + T_snap_2f, my_ghost, Filenames_2f)
	T_mean_t_2f[j] = T_mean[-600]	
	
# Single-fluid temperature perturbation 
T_snap_1f = 5
snap_1f   = np.arange(130, 200, T_snap_1f)
T_mean_t_1f = np.zeros(len(snap_1f))
for j, snap_i in enumerate(snap_1f):

	T_mean = T_mean_1f_func(snap_i, snap_i + T_snap_1f, my_ghost + 1, Filenames_1f)
	T_mean_t_1f[j] = T_mean[-600]
	

plt.close()
plt.figure(figsize = (8,6))
plt.plot(snap_2f/10, T_mean_t_2f, '.k', label = '2f')
plt.plot(snap_1f   , T_mean_t_1f, 'dk', label = '1f')
plt.xlabel('t (s)', fontsize = 14)
plt.ylabel(r'$T_1$ (K)', fontsize = 14)
plt.xlim(min(snap_1f) - 2, max(snap_1f) + 2)
plt.legend(frameon = False, fontsize = 13)
plt.text(142.5, 1.65, r'z ~ 1700 km', fontsize = 13)
plt.minorticks_on()
plt.tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
plt.tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
plt.savefig('../Figures/Tvstime.png')
	

