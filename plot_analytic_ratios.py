from __future__ import division
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import math
from utils import *
from integrals import *
import matplotlib.pyplot as plt

#nice plots
import prettyplotlib as ppl
# This is "import matplotlib.pyplot as plt" from the prettyplotlib library
from prettyplotlib import plt

# This is "import matplotlib as mpl" from the prettyplotlib library
from prettyplotlib import mpl

# Set the random seed for consistency
np.random.seed(12)

from matplotlib import rcParams
rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize'] = 9
rcParams['ytick.labelsize'] = 9
rcParams['legend.fontsize'] = 9
rcParams['legend.frameon'] = False
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['patch.edgecolor'] = 'none'
rcParams['xtick.direction'] = 'in'     # direction: in or out
rcParams['ytick.direction'] = 'in'
rcParams['axes.linewidth']=1.0


fig = plt.figure(figsize=(12, 10))
radius = np.linspace(1, 10000, 100)
R_d = 1000
V_max = 200
Sigma_0 = 1000
for ratio in range(1, 5):
		R_flat = 1/ratio * R_d
		Mass = get_M_integral(radius, Sigma_0, R_d)
		vel = RC_model([V_max, R_flat], radius)
		Sigma = Sigma_model([Sigma_0, R_d], radius)
		#J integral
		J_int, err_int = integrate_J_in_ring(radius, Sigma_0, R_d, V_max, R_flat)	
		#j integral
		j_int = get_j_int(J_int, Mass)
		
		Mass = np.cumsum(Mass)
		Mass = Mass/np.max(Mass)
		J_int = np.cumsum(J_int)
		J_int = J_int/np.max(J_int)
		j_int = np.cumsum(j_int)
		j_int = j_int/np.max(j_int)
	
		vel = vel/V_max
		print radius
		r = radius/R_d
		Sigma = Sigma/np.max(Sigma)
		ax = fig.add_subplot(4, 1, ratio)
		ax.plot(r, Mass, label="Mass")
		ax.plot(r, J_int, label="J")
		ax.plot(r, j_int, label="j")
		ax.plot(r, vel, label="vel")
		ax.plot(r, Sigma, label="$\Sigma$")
		ax.axvline(R_d/R_d, c='k', label='$R_d$')
		ax.axvline(R_flat/R_d, c='b', label='$R_{flat}$')
		plt.title("$R_d/R_{flat}$ ratio: "+str(ratio))
		plt.legend(loc='best')
		plt.xlabel("Radius, R_d")
plt.tight_layout()		
plt.savefig("img/ratios")	
