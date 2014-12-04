from __future__ import division
import numpy as np
import math
from utils import *
from geom import *
import sampler, ensemble
import scipy.stats as stats
import sys
import random

prior = sys.argv[2]
priorType = sys.argv[3]
name = sys.argv[1]
#priorIncl = np.genfromtxt('chains/incl_'+name)



def getPriors(params):
   vc, c, gamma,  pa, incl  = params    
   if (gamma <= 0.8):
       return -np.inf
   elif (gamma > 1.2):
       return -np.inf       
   else:
       return 0

def lnProb(pars, data):
   x, y, vel, evel, r50 = data
   #data = (x, y, r50)  
   priorProb=getPriors(pars)
   vc, c, gamma,  pa, incl  = pars
   m_vel = model2_Courteau(pars, data)
   vel2 = vel.copy()
   vel2 = reshape_array(vel2)
   evel2 = evel.copy()
   evel2 = reshape_array(evel2)
   resid = np.sum((vel2 - m_vel)**2/(evel2)**2)
   lnl = -0.5*((resid))
   return lnl + priorProb



def getRandParams(initParams):
	out = []
	for i, x in enumerate(initParams):
		if i == 0: #vc
			out.append((np.random.normal(loc=1., scale=0.1)*x))
		elif i == 1: #c
			out.append(np.random.normal(loc=1., scale=0.1)*x)
		elif i == 3: #incl
			out.append(np.radians(np.random.normal(loc=1., scale=0.4)*np.degrees(x)))
		else:
			out.append((np.random.normal(loc=1., scale=0.4)*x))
    	return out
    


def getRotCurveFromVelField(params, data):
  name, x, y, vel, evel, r50, GAMA_incl, linewidth, HI_Vc_err = data
  pa, incl, v0 = params  
  x_rot, y_rot = rotateGalaxy(x, y, -pa)  
  y_rot = np.divide(y_rot, np.cos(incl))
  phi = getPhi(y_rot, x_rot)
  #goodAngles = np.where((np.abs(np.cos(phi)) > 0.3))
  #vel = vel[goodAngles]
  #y_rot = y_rot[goodAngles]
  #x_rot = x_rot[goodAngles]
  #phi = phi[goodAngles]
  vel = vel - v0
  vel = np.divide(vel, np.sin(incl))
  vel = np.divide(vel, np.abs(np.cos(phi)))
  radius = np.sqrt((x_rot)**2+(y_rot)**2)
  return np.reshape(radius, (radius.shape[0],)), vel


def model2_Courteau(params, data):
  #takes in a set of x, y coordinates, where y is already transformed by y = y*cos(incl) and both x, y are rotated by -pa
  x, y, vel, evel, r50 = data
  vc, c, gamma, pa, incl = params              # Parameters,  vc, k, gamma, v0, pa, incl 
  #makes a new set of coordinates by rotating the galaxy back
  phi = getPhi(y, x)   
  #deprojects y for the real radius calculation:
  #calculates real, de-projected radius at each point
  radius = np.sqrt(x**2+y**2)
  model_vel_nom = vc*radius*np.cos(phi - pa)*np.sin(incl)*(np.cos(incl))**gamma
  model_vel_denom = (radius**2*((np.sin(phi - pa))**2 + (np.cos(incl))**2*(np.cos(phi - pa))**2) + c**2*(np.cos(incl))**2)**(gamma/2)
  model_vel= model_vel_nom/model_vel_denom 
  return model_vel

def atan_model(params, data):
  #takes in a set of x, y coordinates, where y is already transformed by y = y*cos(incl) and both x, y are rotated by -pa
  x, y, vel, vel_err, r50 = data
  vc, k, pa, incl = params 
  #makes a new set of coordinates by rotating the galaxy back
  x, y = rotateGalaxy(x, y, pa)  #atsukam atgal
  y = y/np.abs(np.cos(incl))     
  phi = getPhi(y, x)   
  radius = np.sqrt(x**2+y**2)
  model_vel = ((2/math.pi) * vc* np.tanh(radius/(k*r50)))*np.sin(incl) 
  model_vel = model_vel * np.cos(phi - pa)
  return np.reshape(model_vel, (model_vel.shape[0], ))

def expModel(params, data):
	x, y, r50 = data #x, y -- transformed (rotated and inclined)
	vmax, c, pa, incl, v0 = params
	x_rot, y_rot = rotateGalaxy(x, y, -pa)
	y_rot = y_rot/np.abs(np.cos(incl))
	phi = getPhi(y_rot, x_rot)
	radius = np.sqrt(x_rot**2 + y_rot**2)
	model_vel = vmax*(1 - np.exp(-1*(radius/(c*r50))))
	model_vel = model_vel*np.cos(phi)*np.sin(incl) + v0	
	#simple_plot(x, y, model_vel, 'img/models/model_test.png')
	return np.reshape(model_vel, (model_vel.shape[0], ))

		

def lnProbAtan(pars, data):
   x, y, vel, evel, r50 = data
   vc, k, pa, incl = pars 
   if k <= 0:  
       return -np.inf
   m_vel = atan_model(pars, data)
   resid = np.sum((vel - m_vel)**2/(evel)**2)
   lnl = -0.5*((resid))
   return lnl


def getTwoLevelPrior(pars, data):
   name, x, y, vel, evel, r50,GAMA_incl, linewidth, HI_Vc_err = data
   vmax, c, pa, incl, v0 = pars
   incl_prior_sampling = np.random.choice(priorIncl)
   LW = linewidth/(2*math.sin(incl_prior_sampling))
   VmaxPrior = getVmaxPriorProb(vmax, LW, HI_Vc_err)
   return VmaxPrior 

def getPriorWidth(HI_Vc_err):
  priorWidth =  5+HI_Vc_err
  return priorWidth
  
def getProbDistribution(linewidth, HI_Vc_err, sigma):
    lower = linewidth - getPriorWidth(HI_Vc_err)
    upper = linewidth + getPriorWidth(HI_Vc_err)
    probDist = stats.truncnorm((lower - linewidth) / sigma, (upper - linewidth) / sigma, loc=linewidth, scale=sigma)
    return probDist
 
def getInclProbDistribution(incl, cutoff, sigma):
    lower = linewidth - cutoff
    upper = linewidth + cutoff
    probDist = stats.truncnorm((lower - incl) / sigma, (upper - incl) / sigma, loc=incl, scale=sigma)
    return probDist


def getInclPrior(GAMA_incl, incl):
    sigma=0.1   
    probDist = getProbDistribution(np.cos(GAMA_incl), 0.5, sigma)
    priorProb = probDist.pdf(np.cos(incl))   
    return np.log(priorProb)

def getVmaxPriorProb(vmax, LW, HI_Vc_err):
    sigma=getPriorWidth(HI_Vc_err)   
    probDist = getProbDistribution(LW, HI_Vc_err, sigma)
    priorProb = probDist.pdf(vmax)   
    return np.log(priorProb)  

def getLinewidthPrior(pars, data):
   name, x, y, vel, evel, r50, GAMA_incl, linewidth, HI_Vc_err = data
   vmax, c, pa, incl, v0 = pars
   v22 = vmax*(1 - np.exp(-1*(200/(c*r50))))
   LW = linewidth/(2*math.sin(GAMA_incl))
   VmaxPrior = getVmaxPriorProb(v22, LW, HI_Vc_err)
   return VmaxPrior


def lnProbExp(pars, data):
   name, x, y, vel, vel_err, r50, GAMA_incl, HI_linewidth, HI_Vc_err = data
   vmax, c, pa, incl, v0 = pars
   model_data = (x, y, r50) 
   m_vel = expModel(pars, model_data)
   resid = np.sum((vel - m_vel)**2/(vel_err)**2)
   lnl = -0.5*((resid))
   if prior =='True':
		if priorType == '2':
			dataPrior = getInclPrior(GAMA_incl, incl)
		elif priorType == '3':
			dataPrior = getInclPrior(GAMA_incl, incl) + getLinewidthPrior(pars, data)
		elif priorType == '4':
			dataPrior = getLinewidthPrior(pars, data)
		return lnl  + dataPrior
   else:
		return lnl

