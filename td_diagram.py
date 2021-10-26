import numpy as np
import matplotlib.pyplot as plt
import h5py 
import glob 
from configparser import ConfigParser

# Read td_diagram.ini
params = ConfigParser()
params.sections()
params.read('../Doc_data/Config/td_diagram.ini')
Params = params['params']

# Function to make the plots
def plotter(F, F_label, figname, ti, tf, mz, dz):

	plt.close()
	fig, ax = plt.subplots(figsize = (8,6))
	im = ax.imshow(F.transpose(), aspect='auto', origin='lower', cmap = plt.cm.bwr_r, 
                   extent = (ti, tf, 520, (mz - 1) * dz * 1e-3 + 520))
	cbar = plt.colorbar(im)
	cbar.set_label(F_label, fontsize = 14)
	cbar.ax.tick_params(labelsize=13) 
	plt.xlabel('t (s)', fontsize = 14)
	plt.ylabel('Height (km)', fontsize = 14)
	plt.tick_params(axis ='both', direction='out', which='minor',
                    length=3, width=.5,labelsize=14, top = True, right = True)
	plt.tick_params(axis='both', direction='out', which='major',
                    length=8, width=1, labelsize=14, top = True, right = True)
	plt.minorticks_on()	
	plt.tight_layout()
	plt.savefig(figname)

	

# Read the parameters
H5dir      = Params['H5dir']
ranges     = Params['ranges'].split(',')
Outname    = Params['Outname']
eqfilename = Params['eqfilename']

ranges = [int(ranges[0]), int(ranges[1])]

# Read files
Filenames = sorted(glob.glob(H5dir + '*'))[ranges[0]:ranges[1]]

# Read the equilibrium file
with h5py.File(eqfilename, 'r') as H5obj0:

	# Ghost shells
	my_ghost0 = int(H5obj0.attrs['my_ghost'])
	# Obtain data (1f)
	Pc0   = np.array(H5obj0['pe_c'  ])[my_ghost0:-my_ghost0,0,0]
	Pn0   = np.array(H5obj0['pe_n'  ])[my_ghost0:-my_ghost0,0,0]
	Dc0   = np.array(H5obj0['rho_c' ])[my_ghost0:-my_ghost0,0,0]
	Dn0   = np.array(H5obj0['rho_n' ])[my_ghost0:-my_ghost0,0,0]

# Read the perturbations
Vn, Vc, Pn, Pc, Tn, Tc, Dn, Dc, B, = [], [], [], [], [], [], [], [], []
for i, H5file in enumerate(Filenames):

	with h5py.File(H5file, 'r') as h5obj:
	
		if i == 0:
		
			dz       = h5obj.attrs['dz'][0]
			my_ghost = h5obj.attrs['my_ghost'][0]
			mz       = h5obj.attrs['dims'][2]
			ti       = h5obj.attrs['time'][0]
			
		elif i == (len(Filenames) - 1):
		
			tf       = h5obj.attrs['time'][0]
	    	
		Vn.append(h5obj['vz_n'  ][my_ghost:-my_ghost,0,0])
		Vc.append(h5obj['vz_c'  ][my_ghost:-my_ghost,0,0])
		Pn.append(h5obj['pe_n'  ][my_ghost:-my_ghost,0,0])
		Pc.append(h5obj['pe_c'  ][my_ghost:-my_ghost,0,0])
		Tn.append(h5obj['temp_n'][my_ghost:-my_ghost,0,0])
		Tc.append(h5obj['temp_c'][my_ghost:-my_ghost,0,0])
		Dn.append(h5obj['rho_n' ][my_ghost:-my_ghost,0,0])
		Dc.append(h5obj['rho_c' ][my_ghost:-my_ghost,0,0])
		B.append( h5obj['bx'    ][my_ghost:-my_ghost,0,0])
	    
# Transform the list into numpy arrays
Vn = np.array(Vn); Vc = np.array(Vc)
Pn = np.array(Pn); Pc = np.array(Pc)
Tn = np.array(Tn); Tc = np.array(Tc)
Dn = np.array(Dn); Dc = np.array(Dc)
B  = np.array(B )

# Plot the results
plotter(Vn      , r'$u_{n}~(m/s)$'    , Outname + '_un.png', ti, tf, mz, dz)
plotter(Vc      , r'$u_{c}~(m/s)$'    , Outname + '_uc.png', ti, tf, mz, dz)
plotter(Pn / Pn0, r'$P_{n} / P_{n,0}$', Outname + '_pn.png', ti, tf, mz, dz)
plotter(Pc / Pc0, r'$P_{c} / P_{c,0}$', Outname + '_pc.png', ti, tf, mz, dz)
plotter(Tn      , r'$T_{n}~(K)$'      , Outname + '_tn.png', ti, tf, mz, dz)
plotter(Tc      , r'$T_{c}~(K)$'      , Outname + '_tc.png', ti, tf, mz, dz)
plotter(Dn / Dn0, r'$D_{n} / D_{n,0}$', Outname + '_dn.png', ti, tf, mz, dz)
plotter(Dc / Dc0, r'$D_{c} / D_{c,0}$', Outname + '_dc.png', ti, tf, mz, dz)
plotter(B       , r'$B~(T)$'          , Outname + '_bx.png', ti, tf, mz, dz)
