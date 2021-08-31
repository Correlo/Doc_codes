import numpy as np
import h5py
from scipy.interpolate import interp1d
from configparser import ConfigParser
# Packages from the developer
from const import *

#Read params.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/iniCond.ini')
Params = params['params']

# Make sure filename has h5 extension otherwise visit won't recognize it!
eqfilename = Params['eqfilename']
valc_file  = Params['valc_file'] 

# Configuration parameters
mx = int(Params['mx'])
my = int(Params['my'])
mz = int(Params['mz'])
my_ghost = int(Params['my_ghost'])
B00 = float(Params['B00'])
oneFluid = bool(Params['oneFluid'])
with_B0x = bool(Params['with_B0x'])
div = float(Params['div'])


# Ghost param
mz_ghost = mz + 2*my_ghost

# VALC
zzValc, tempValc, nnValc, neValc, rhoValc = np.loadtxt(valc_file, usecols=(0,3,5,6,9), unpack=True)


## SUP LIMIT at T=9000K (at z approx 2.1 Mm) ##
indValc = (tempValc <= 9000)
# Fields 
zzValc   = zzValc[indValc]
tempValc = tempValc[indValc]
rhoValc  = rhoValc[indValc]
nnValc   = nnValc[indValc]
neValc   = neValc[indValc]

## INF LIMIT  at 0.5 Mm ##
indValc = (zzValc > 500)
# Fileds
zzValc   = zzValc[indValc]
tempValc = tempValc[indValc]
rhoValc  = rhoValc[indValc]
nnValc   = nnValc[indValc]
neValc   = neValc[indValc]

# Convert valc units to MKS
zzValc  = zzValc  * 1e3 #km to m
rhoValc = rhoValc * 1e3
nnValc  = nnValc  * 1e6
neValc  = neValc  * 1e6

print("Z MIN %e , MAX %e" % (np.min(zzValc), np.max(zzValc)))

# Atmosphere altitude
z, dz = np.linspace(np.min(zzValc), np.max(zzValc), mz_ghost, retstep = True)

## Interpolate values from VALC ##

# Density
func = interp1d(zzValc, rhoValc , kind='cubic')
Rho  = func(z)
# Temp
func = interp1d(zzValc, tempValc, kind='cubic')
Temp = func(z)
# nn
func = interp1d(zzValc, nnValc  , kind='cubic')
nn   = func(z)
# ne
func = interp1d(zzValc, neValc  , kind='cubic')
ne   = func(z)


# Get values at the bottom of the atmosphere in the VALC model
nn00 = nn[0]
ne00 = ne[0]
nc00 = 2*ne00

print("nc00 %e" % (nc00))
print("nn00 %e" % (nn00))


# Set charges and neutrals temperature equal to the VALC temperature
Temp_c = Temp
Temp_n = Temp

pe_n00 = nn00 * BK * Temp_n[0]
pe_c00 = nc00 * BK * Temp_c[0]

Temp = 0.5 * (Temp_n + Temp_c)
alpha = MIH / MH**2 * np.sqrt((8 * BK * Temp) / (np.pi * MIH)) * SIGMA_IN + \
        MEH / MH**2 * np.sqrt((8 * BK * Temp) / (np.pi * MEH)) * SIGMA_EN

print("pec00 %e" % (pe_c00))
print("pen00 %e" % (pe_n00))
print("BETA PLASMA %e " % (pe_c00 / (B00**2 / (2 * MU0))))
print("BETA PLASMA TOTAL %e " % ((pe_c00 + pe_n00)/(B00**2 / (2 * MU0))))

pe_n0  = np.zeros((mz_ghost,my,mx))
pe_c0  = np.zeros((mz_ghost,my,mx))
rho_n0 = np.zeros((mz_ghost,my,mx))
rho_c0 = np.zeros((mz_ghost,my,mx))
B0x    = np.zeros((mz_ghost,my,mx))

print("Bx00 %e" % (B00))
print("Va0 %e" % np.sqrt(B00**2/(MU0 * MH* (nn00 + ne00) )))

print("G = ", G)
print("Temp_n = ", np.min(Temp_n), np.max(Temp_n), np.mean(Temp_n))
print("Temp_c = ", np.min(Temp_c), np.max(Temp_c), np.mean(Temp_c))

# Integrate MHS equilibrium

for i in range(mx):
	for j in range(my):
		sumTemp_c = 0.0
		sumTemp_n = 0.0
		for k in range(mz_ghost):
			expArg = (MH * G * dz * sumTemp_c * pe_c00)/(2.0 * BK * (2 * pe_c00 + B00**2/MU0))
			pe_n0[k,j,i] = pe_n00 * np.exp(-(MH * G * dz * sumTemp_n)/BK)
			pe_c0[k,j,i] = pe_c00 * np.exp(-2*expArg)	
			B0x[k,j,i]=B00*np.exp(-expArg)
			#densities
			rho_n0[k,j,i] = pe_n0[k,j,i] * MH / (      BK * Temp_n[k])
			rho_c0[k,j,i] = pe_c0[k,j,i] * MH / (2.0 * BK * Temp_c[k])
			sumTemp_c = sumTemp_c + 1.0/ Temp_c[k]
			sumTemp_n = sumTemp_n + 1.0/ Temp_n[k]

#make pMag = pe_n0[eqP]
eqP = int(mz_ghost / div) 
pmag = 0.5 * B0x**2 / MU0
const = pe_n0[eqP,0,0] - pmag[eqP,0,0]
B0x = np.sqrt(2 * MU0 * (pmag + const))

zero = np.zeros((mz_ghost,my,mx), dtype = np.longdouble)

with h5py.File(eqfilename, 'w') as f:
		f.attrs['dx'] = dz
		f.attrs['dy'] = dz
		f.attrs['dz'] = dz
		f.attrs['my_ghost'] = my_ghost
		f.attrs['time'] = 0.0
		f.attrs['dims'] = [mx,my,mz]
		#!!! lwz  compression does not seem to work, do no use it
		if oneFluid:
			dset = f.create_dataset("pe_c" , (mz_ghost, my, mx), dtype='f8', data=pe_c0+pe_n0  , chunks=True, shuffle=True)
			dset = f.create_dataset("rho_c", (mz_ghost, my, mx), dtype='f8', data=rho_c0+rho_n0, chunks=True, shuffle=True)
		else:
			dset = f.create_dataset("pe_c" , (mz_ghost, my, mx), dtype='f8', data=pe_c0 , chunks=True, shuffle=True)
			dset = f.create_dataset("rho_c", (mz_ghost, my, mx), dtype='f8', data=rho_c0, chunks=True, shuffle=True)
			dset = f.create_dataset("pe_n" , (mz_ghost, my, mx), dtype='f8', data=pe_n0 , chunks=True, shuffle=True)
			dset = f.create_dataset("rho_n", (mz_ghost, my, mx), dtype='f8', data=rho_n0, chunks=True, shuffle=True)
		# Mag fiels
		if with_B0x:
			dset = f.create_dataset("bx", (mz_ghost, my, mx),  dtype='f8', data=B0x, chunks=True, shuffle=True)







