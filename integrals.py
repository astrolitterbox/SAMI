from __future__ import division
from scipy import integrate
import math
import numpy as np



def Mass_integral(r, S1, S2, R1, R2):
	Sigma_r = Sigma_model_2exp([S1, S2, R1, R2], r)
	ret = 2*math.pi*Sigma_r * r
	return ret

def get_M_integral(radius, S1, S2, R1, R2):
	M_int = np.zeros((radius.shape))
	for i, r in enumerate(radius):
		integral, err = integrate.quad(Mass_integral, radius[i-1], r, args=(S1, S2, R1, R2))
		M_int[i] = integral
	M_int[0] = integrate.quad(Mass_integral, radius[i-1], r, args=(S1, S2, R1, R2))[0]
	return M_int
	
def J_integral(r, S1, S2, R1, R2, Vmax, Rflat):
	Sigma_r = Sigma_model_2exp([S1, S2, R1, R2], r)
	vel_r = RC_model([Vmax, Rflat], r)
	ret = r**2 * Sigma_r * vel_r
	return 2*math.pi*ret

def Sigma_model(params, radius):
	S, R = params
	return S*np.exp(-radius/R)

def Sigma_model_2exp(params, radius):
	S1, S2, R1, R2 = params
	return S1*np.exp(-radius/R1) + S2*np.exp(-radius/R2)
	
def RC_model(params, radius):
	Vmax, Rflat = params
	ret = Vmax*(1 - np.exp(-1*radius/Rflat))
	return ret

def get_j_int(J_int, M_integral):
	ret = J_int/M_integral
	ret[np.where(M_integral == 0)] = 0
	return ret

def get_j(J, ring_mass):
	ret = J/ring_mass
	ret[np.where(ring_mass == 0)] = 0
	return ret

def integrate_J_in_ring(radius, S1, S2, R1, R2, Vmax, Rflat):
		J_int = np.zeros((radius.shape))
		err_int = np.zeros((radius.shape))
		for i, r in enumerate(radius):
			integral, err = integrate.quad(J_integral, radius[i-1], r, args=(S1, S2, R1, R2, Vmax, Rflat))
			J_int[i] = integral
			Sigma_r = Sigma_model_2exp([S1, S2, R1, R2], r)
			vel_r = RC_model([Vmax, Rflat], r)
			err_int[i] = err
		J_int[0], err_int[0] = integrate.quad(J_integral, 0, radius[0], args=(S1, S2, R1, R2, Vmax, Rflat))
		return J_int, err_int

def get_M_in_ring(radius, Sigma):
		ring_area = np.zeros((radius.shape))
		for i, r in enumerate(radius):
			ring_area[i] = (radius[i]**2 - radius[i-1]**2)
		ring_area[0] = (radius[0]**2 - 0) # the 0th ring
		ring_area = math.pi*ring_area
		ring_mass = ring_area * Sigma
		return ring_mass



def get_J_in_ring(radius, Sigma, vel):
		ring_mass = get_M_in_ring(radius, Sigma)
		J = ring_mass * vel * radius
		return J

def get_j_in_ring(radius, vel):
		j = vel * radius
		return j

def get_model_j(radius, Rflat, Vmax):
	j = (2 + ((radius/Rflat)**2)/(1 + (radius/Rflat) - np.exp(radius/Rflat)))*Rflat*Vmax
	return j

def get_final_model_j(Rflat, Vmax):
	j = (2 + ((1000/Rflat)**2)/(1 + (1000/Rflat) - np.exp(1000/Rflat)))*Rflat*Vmax
	return j
