import numpy as np
import matplotlib.pyplot as plt
import h5py
from configparser import ConfigParser
# Modules from the developer
from const import *

# Read ne.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/ne.ini')
Params = params['params']

# Read params
eqfilename = Params['eqfilename']
outname    = Params['outname']

# Read the density of charges
with h5py.File(eqfilename, 'r') as h5Obj:

	dz  = h5Obj.attrs['dz']
	mx, my, mz = h5Obj.attrs['dims']
	my_ghost = int(h5Obj.attrs['my_ghost'])
	Dc   = np.array(h5Obj['rho_c' ])[my_ghost:-my_ghost,0,0]
	
# Obtain the density of electrons
ne = Dc / MH

# Save results 
with h5py.File(outname, 'w') as f:

	f.attrs['dx']   = dz
	f.attrs['dy']   = dz
	f.attrs['dz']   = dz
	f.attrs['dims'] = [mx, my, mz]
	f.attrs['my_ghost'] = my_ghost
	dset = f.create_dataset("ne", (len(ne), 1, 1),  dtype='f8', data=ne, chunks=True, shuffle=True)
	

	


