from __future__ import division
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import math
from utils import *
from nicer_plots import *
import pyfits
from astroML.plotting import hist
from astropy.coordinates.distances import Distance
from SAMI_Mass_fn import *


data = np.genfromtxt("data/SAMI_Master_Catalogue_20141009.txt", names=['SAMI_ID', 'ra', 'dec', 'z', 'M_st', 'R_e'], usemask=False, dtype="a40,f8,f8, f8,f8, f8")
ids = np.genfromtxt("data/sami_files.txt", dtype="object")

sami_id = data['SAMI_ID']
all_z = data['z']
all_R_e = data['R_e']
all_M_st = data['M_st']

R_e = np.zeros((ids.shape))
z = np.zeros((ids.shape))
distances = np.zeros((ids.shape))
M_st = np.zeros((ids.shape))



for i, sami_id in enumerate(ids):
	name = str(sami_id[:-15])
	ind = np.where(data['SAMI_ID'] == name)
	R_e[i] = (np.round(all_R_e[ind], 1))
	z[i] = all_z[ind]
	distances[i] = Distance(z=z[i]).Mpc
	M_st[i] = all_M_st[ind]

phys_R_e = angular2physical(R_e, z)

fig = plt.figure()
hist(R_e, bins='knuth')
plt.xlabel("$R_e$ [arcsec]")
plt.savefig("img/R_e_hist")

fig = plt.figure()
hist(z, bins='knuth', label='SAMI')
plt.axvspan(0.0005, 0.003, facecolor='r', alpha=0.4, label='CALIFA')
plt.xlabel("z")
plt.legend(loc='best')
plt.savefig("img/z_hist")

fig = plt.figure()
hist(distances, bins='knuth')
plt.xlabel("Distance [Mpc]") 
plt.savefig("img/distances_hist")

fig = plt.figure()
hist(distances, bins='knuth')
plt.xlabel("Distance [Mpc]") 
plt.savefig("img/distances_hist")


fig = plt.figure()
hist(phys_R_e, bins='knuth')
plt.xlabel("Physical $R_e$ (kpc)") 
plt.savefig("img/phys_R_e_hist")


fig = plt.figure()
n, bins, patches = hist(M_st, bins='knuth')
plt.xlabel("$log(M_{st}), [M_{\odot}]$") 
plt.savefig("img/M_st_hist")

M, Gama_mass_fn = get_GAMA_Mass_fn()


fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = np.histogram(M, bins=20)

ax.plot(np.log10(M), Gama_mass_fn)
#M_st = 10**M_st	
bin_centres = (bins[:-1] + bins[1:])/2.

hist(M_st, normed=True) 

ax.set_yscale('log')
#ax.set_xscale('log')
plt.savefig("img/Gama_M_fn")

