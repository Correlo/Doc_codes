#********************************************************************************************#
#    Code to plot Velocity, Pressure, Temperature and Density results from Acoustic_WaveD    #
#********************************************************************************************#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import glob
from configparser import ConfigParser

# Definitions

def animate(i):

    # Open file
    with h5py.File(Filenames_2f[i], 'r') as H5obj_2f:
	
    	# Attributes
    	dz       = H5obj_2f.attrs['dz']
    	t        = H5obj_2f.attrs['time']
    	my_ghost = H5obj_2f.attrs['my_ghost'][0]
    	mz       = H5obj_2f.attrs['dims'][2]
    	# Obtain data
    	vz_c = np.array(H5obj_2f['vz_c'  ])[my_ghost:-my_ghost, 0, 0]
    	vz_n = np.array(H5obj_2f['vz_n'  ])[my_ghost:-my_ghost, 0, 0]
    	Pc   = np.array(H5obj_2f['pe_c'  ])[my_ghost:-my_ghost, 0, 0]
    	Pn   = np.array(H5obj_2f['pe_n'  ])[my_ghost:-my_ghost, 0, 0]
    	Tc   = np.array(H5obj_2f['temp_c'])[my_ghost:-my_ghost, 0, 0]
    	Tn   = np.array(H5obj_2f['temp_n'])[my_ghost:-my_ghost, 0, 0]
    	Dc   = np.array(H5obj_2f['rho_c' ])[my_ghost:-my_ghost, 0, 0]
    	Dn   = np.array(H5obj_2f['rho_n' ])[my_ghost:-my_ghost, 0, 0]
    	B_2f = np.array(H5obj_2f['bx'    ])[my_ghost:-my_ghost, 0, 0]
    	
    with h5py.File(Filenames_1f[i], 'r') as H5obj_1f:
    
    	# Attributes
    	my_ghost = H5obj_1f.attrs['my_ghost'][0]
    	# Obtain data
    	vz   = np.array(H5obj_1f['vz_c'  ])[my_ghost:-my_ghost, 0, 0]
    	P    = np.array(H5obj_1f['pe_c'  ])[my_ghost:-my_ghost, 0, 0]
    	T    = np.array(H5obj_1f['temp_c'])[my_ghost:-my_ghost, 0, 0]
    	D    = np.array(H5obj_1f['rho_c' ])[my_ghost:-my_ghost, 0, 0]
    	B_1f = np.array(H5obj_1f['bx'    ])[my_ghost:-my_ghost, 0, 0]

    # z-axis
    Z = np.arange(mz) * dz * 1e-3  + 520 # Km

    # Set data
    Line[0 ].set_data(Z, vz_c); Line[1 ].set_data(Z, vz_n); Line[2 ].set_data(Z, vz)
    Line[3 ].set_data(Z, Pc  ); Line[4 ].set_data(Z, Pn  ); Line[5 ].set_data(Z, P )
    Line[6 ].set_data(Z, Tc  ); Line[7 ].set_data(Z, Tn  ); Line[8 ].set_data(Z, T )
    Line[9 ].set_data(Z, Dc  ); Line[10].set_data(Z, Dn  ); Line[11].set_data(Z, D )
    Line[12].set_data(Z, B_2f); Line[13].set_data(Z, B_1f)

    # Update time
    Line[14].set_text('t = %.1f s' % t)

    if i == 0:

        # Xlim
        ax[0].set_xlim(min(Z), max(Z))
        
        # Ylim
        ax[0].set_ylim(-1.2e3, 1.2e3)
        ax[1].set_ylim(-1.6e1, 1.6e1)
        ax[2].set_ylim(-4e2, 8e2)
        ax[3].set_ylim(-2.7e-7, 2.7e-7)
        ax[4].set_ylim(-3e-4, 2e-4)

    return Line

# Read VPTDBgif.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/VPTDBgif.ini')
Params = params['params']

# Read the params
H5dir_2f = Params['H5dir_2f']
H5dir_1f = Params['H5dir_1f']
figname  = Params['figname' ]


# Obtain the name of hdf5 file
Filenames_2f = sorted(glob.glob(H5dir_2f + '*'))
Filenames_1f = sorted(glob.glob(H5dir_1f + '*'))

# create a figure with four subplots and prepare configuration
plt.close()
fig, ax = plt.subplots(5, 1, sharex = True, figsize = (12, 8))
fig.subplots_adjust(hspace = 0)

#Axis labels
ax[0].set_ylabel(r'$u$ $(m$ $s^{-1})$'      , fontsize = 14)
ax[1].set_ylabel(r'$P_1$ $(Pa)$'            , fontsize = 14)
ax[2].set_ylabel(r'$T_1$ $(K)$'             , fontsize = 14) 
ax[3].set_ylabel(r'$\rho_1$ $(kg$ $m^{-3})$', fontsize = 14)
ax[4].set_ylabel('$B_1$ $(Km)$'             , fontsize = 14) 
ax[4].set_xlabel('Altitude (Km)'            , fontsize = 14)

# Sci tick labels
for i in range(5): ax[i].ticklabel_format(axis = "y", style = "sci", scilimits = (0,0))
# title
title = ax[0].set_title('t = %.1f s' % 0.0, fontsize = 15, pad = 14)

# List of input plot objects to animation function
Line = []
line, = ax[0].plot([], [], color='r', label = 'charges')
Line.append(line)
line, = ax[0].plot([], [], color='b', label = 'neutrals')
Line.append(line)
line, = ax[0].plot([], [], color='g', label = '1f ambipolar')
Line.append(line)

# Legend
ax[0].legend(fontsize = 13, loc = 'upper right')

for i in range(1,4):
    line, = ax[i].plot([], [], color='r')
    Line.append(line)
    line, = ax[i].plot([], [], color='b')
    Line.append(line)
    line, = ax[i].plot([], [], color='g')
    Line.append(line)
    
line, = ax[4].plot([], [], color='y')
Line.append(line)
line, = ax[4].plot([], [], color='g')
Line.append(line)

# Add title at the end
Line.append(title)

ax[0].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, top = True, bottom = False)
ax[0].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[1].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, bottom = False)
ax[1].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[2].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, bottom = False)
ax[2].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[3].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True, bottom = False)
ax[3].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)
ax[4].tick_params(axis ='both', direction='inout', which='minor',
                 length=3, width=.5,labelsize=13, right = True)
ax[4].tick_params(axis='both', direction='in', which='major',
                 length=6, width=1, labelsize=13, right = True, top = True)

for i in range(5): ax[i].ticklabel_format(axis = "y", style = "sci", scilimits = (0,0))
for i in range(5): ax[i].minorticks_on()

ani = animation.FuncAnimation(fig, animate, interval = 100, frames = len(Filenames_2f),
                              blit = False, repeat = False)
ani.save(figname, writer = 'imagemagick', fps = 100)
