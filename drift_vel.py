import numpy as np
import matplotlib.pyplot as plt
import h5py
from configparser import ConfigParser
# Modules from the developer
from const import *

# Read drift_vel.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/drift_vel.ini')
Params = params['params']

# Read params
H5file  = Params['H5file']
outfile = Params['outfile']

# Read the velocities
with h5py.File(H5file, 'r') as H5obj:

	dz  = H5obj.attrs['dz']
	mx, my, mz = H5obj.attrs['dims']
	my_ghost = int(H5obj.attrs['my_ghost'])
	vz_n   = np.array(H5obj['vz_n' ])[my_ghost:-my_ghost,0,0]
	vz_c   = np.array(H5obj['vz_c' ])[my_ghost:-my_ghost,0,0]
	
# z-axis
Z  = np.arange(mz) * dz * 1e-3 + 520 # km

plt.close()
plt.figure(figsize = (12,6))
plt.plot(Z, vz_n - vz_c, '-k')
plt.xlabel('Height (km)', fontsize = 14)
plt.ylabel('w (m/s)', fontsize = 14)
plt.xlim(min(Z), max(Z))
plt.minorticks_on()
plt.tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, top = True, right = True)
plt.tick_params(axis='both', direction='in', which='major',
                 length=8, width=1, labelsize=13, top = True, right = True)
                 
plt.savefig(outfile)



	
