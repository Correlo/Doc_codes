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
			Tc = np.array(H5obj['temp_c'][my_ghost:-my_ghost,0,0])
	
		up, down = Peaks(Z, Tc)
		up_interp   = np.interp(Z, up[:,0]  , up[:,1]  )
		down_interp = np.interp(Z, down[:,0], down[:,1])
		T = (up_interp + down_interp) / 2
	
		T_mean = T_mean + T

	T_mean = T_mean / 50
	
	return T_mean
	
def T_mean_2f_func2(snap_i, snap_f, Filenames):

	# Read density and temperature perturbations
	T_mean_c = np.zeros_like(Z)
	T_mean_n = np.zeros_like(Z)
	
	for i in range(snap_i, snap_f):
		with h5py.File(Filenames[i], 'r') as H5obj:
	
			# Ghost shells
			if (i == snap_i): my_ghost = int(H5obj.attrs['my_ghost'])
			Tc = np.array(H5obj['temp_c'][my_ghost:-my_ghost,0,0])
			Tn = np.array(H5obj['temp_n'][my_ghost:-my_ghost,0,0])
	
		up_c, down_c  = Peaks(Z, Tc); up_n, down_n = Peaks(Z, Tn)
		up_interp_c   = np.interp(Z, up_c[:,0]  , up_c[:,1]  )
		up_interp_n   = np.interp(Z, up_n[:,0]  , up_n[:,1]  )
		down_interp_c = np.interp(Z, down_c[:,0], down_c[:,1])
		down_interp_n = np.interp(Z, down_n[:,0], down_n[:,1])
		Tc = (up_interp_c + down_interp_c) / 2
		Tn = (up_interp_n + down_interp_n) / 2
	
		T_mean_c = T_mean_c + Tc
		T_mean_n = T_mean_n + Tn

	T_mean_c = T_mean_c / 50
	T_mean_n = T_mean_n / 50
	
	return [T_mean_c, T_mean_n]
	
# Read TvsZ.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/TvsZ.ini')
Params = params['params']

# Obtain the filenames
Filenames_2f = sorted(glob.glob(Params['H5dir_2f'] + '*.h5'))
Filenames_1f = sorted(glob.glob(Params['H5dir_1f'] + '*.h5'))
eqfilename   = Params['eqfilename']
figname      = Params['figname']
div          = float(Params['div'])
temp3        = bool(int(Params['temp3']))
ranges       = Params['ranges'].split(',')
f_step       = int(Params['f_step'])
pmlFraction  = float(Params['pmlFraction'])


# Build the range list
ranges = [int(ranges[0]), int(ranges[1])]
low_r = ranges[0]
ranges_l = []

while low_r < ranges[1]:

	ranges_l.append((low_r, low_r + f_step))
	low_r += f_step
	
print(ranges_l)

with h5py.File(eqfilename, 'r') as H5obj0:

	# Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	# Obtain dz
	dz  = H5obj0.attrs['dz']
	Dn0 = np.array(H5obj0['rho_n'][my_ghost0:-my_ghost0,0,0])
	Dc0 = np.array(H5obj0['rho_c'][my_ghost0:-my_ghost0,0,0])
	mz  = H5obj0.attrs['dims'][2]
	
Z    = np.arange(mz) * dz * 1e-3 + 520 # Km


# Plots
plt.close()
plt.figure(figsize = (8,6))

plt.axvline(Z[int(mz/div)], linestyle = '--', color = 'k')
plt.axvline(Z[- int(pmlFraction * mz)], linestyle = '-', color = 'k')

colormap = plt.cm.get_cmap('gist_rainbow', 256)
colors   = colormap(np.linspace(0, 1, len(ranges_l)))

for i, r in enumerate(ranges_l[::-1]):

	if temp3:

		T_mean_2f_c_0, T_mean_2f_n_0 = T_mean_2f_func2(r[0], r[1], Filenames_2f)
		plt.plot(Z, T_mean_2f_c_0, '-.', color = colors[i])
		plt.plot(Z, T_mean_2f_n_0, '--', color = colors[i])

	else:

		T_mean_2f_0 = T_mean_2f_func(r[0], r[1], Filenames_2f)
		plt.plot(Z, T_mean_2f_0, '--', color = colors[i])
	
	T_mean_1f_0 = T_mean_1f_func(r[0], r[1], Filenames_1f)
	plt.plot(Z, T_mean_1f_0, '-', label = '%d - %d s' % (r[0] / 10, r[1] / 10), color = colors[i])

plt.xlabel('Height (km)', fontsize = 14)
plt.ylabel(r'$T_1$ (K)', fontsize = 14)
plt.xlim(min(Z), max(Z))
plt.legend(frameon = False, fontsize = 13)
plt.minorticks_on()
plt.tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
plt.tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
plt.savefig(figname)
	

