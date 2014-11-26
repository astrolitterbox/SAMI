from __future__ import division
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import math
from utils import *
from nicer_plots import *

def get_GAMA_Mass_fn():

	M = np.linspace(10**8, 10**11.5, 20)
	phi1 = 3.96*10**-3
	alpha1 = -0.35
	phi2 = 0.79e-3
	alpha2 = -1.47
	M_star = 10**10.66
	Gama_mass_fn = np.exp((-M/M_star)) * (phi1 *(M/M_star)**alpha1 + phi2*(M/M_star)**alpha2) * M/M_star
	return M, Gama_mass_fn


