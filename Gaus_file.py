import numpy as np
import h5py
from configparser import ConfigParser

def Gauss(z, mu, sigma, N):

	# N   = 1 / (sigma * np.sqrt(2 * np.pi))
	EXP = np.exp(- (z - mu)**2 / (2 * sigma**2))
	G = N * EXP
	return G

# Read iniCond.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/Gaus_file.ini')
Params = params['params']



# Configuration parameters
eqfilename  = Params['eqfilename']
reffile     = Params['reffile']
mx          = int(Params['mx'])
my          = int(Params['my'])
mz          = int(Params['mz'])
my_ghost    = int(Params['my_ghost'])
pmlFraction = float(Params['pmlFraction'])
sigma_p     = float(Params['sigma_p'])    
N           = float(Params['N'])

# Ghost param
mz_ghost = mz + 2*my_ghost

# Obtain dz 
with h5py.File(reffile, 'r') as f: dz = f.attrs['dz']

# Gauss parameters
z  = np.arange(mz_ghost)
mu = mz_ghost - int(mz_ghost * pmlFraction) - 1
sigma = int(mz_ghost * pmlFraction) * sigma_p

# Guassian function
Gaus = np.zeros(mz_ghost)
Gaus[mu:] = Gauss(z[mu:], mu, sigma, N)




with h5py.File(eqfilename, 'w') as f:
		f.attrs['dx'] = dz
		f.attrs['dy'] = dz
		f.attrs['dz'] = dz
		f.attrs['my_ghost'] = my_ghost
		f.attrs['time'] = 0.0
		f.attrs['dims'] = [mx,my,mz]
		
		dset = f.create_dataset("gaus" , (mz_ghost, my, mx), dtype='f8', data=Gaus, chunks=True, shuffle=True)
		






