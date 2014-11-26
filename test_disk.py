from __future__ import division
import numpy as np
import math
from utils import *
from geom import *
from mcmc import *

def make_exp_test_disk(xmax, ymax, modelParams):
	vmax, c, pa, incl = modelParams
	x = []
	y = []
	for i, yi in enumerate(np.linspace(-ymax, ymax, 2*ymax+1)):
		 for j, xi in enumerate(np.linspace(-xmax, xmax, 2*xmax+1)):
			  #print i, j, xi, yi
			  x.append(xi)
			  y.append(yi)
			  #print yi, xi 
	x = np.reshape(x, (len(x), 1))
	y = np.reshape(y, (len(y), 1))
	#cutting out a round disk
	radius = np.where(np.sqrt(x**2 + y**2) < xmax)  
	x = x[radius]
	y = y[radius]
	X = np.empty((x.shape[0], len(incl)))
	Y = X.copy()
	model_vel = X.copy()  
	for l in range(0, len(incl)):
		  print l
		  X[:, l] = x
		  Y[:, l] = y
	for i, inclination in enumerate(incl): 
		  print np.degrees(pa), np.degrees(inclination), 'pa, incl at make disk'		
		  Y[:, i] = Y[:, i]*np.cos(inclination) #projection
		  #rotating by PA, projecting
		  X[:, i], Y[:, i] = rotateGalaxy(X[:, i], Y[:, i], pa) 
		  data = (X[:, i], Y[:, i], 5) #x, y, vel, vel_err, r50 = data
		  model_vel[:, i] = expModel(modelParams, data) 
	#simple_plot(X, Y, model_vel, 'img/models/exp_test.png')
	return X, Y, model_vel, 0.1*model_vel	

def make_test_disk(xmax, ymax, modelParams):
  x = []
  y = []
  vc, c, gamma, pa, incl = modelParams
  for i, yi in enumerate(np.linspace(-ymax, ymax, 2*ymax+1)):
     for j, xi in enumerate(np.linspace(-xmax, xmax, 2*xmax+1)):
          #print i, j, xi, yi
          x.append(xi)
          y.append(yi)
          #print yi, xi 
  x = np.reshape(x, (len(x), 1))
  y = np.reshape(y, (len(y), 1))
      #cutting out a round disk
  radius = np.where(np.sqrt(x**2 + y**2) < xmax)  
  x = x[radius]
  y = y[radius]
  
  X = np.empty((x.shape[0], len(incl)))
  Y = X.copy()
  model_vel = X.copy()  
  for l in range(0, len(incl)):
      print l
      X[:, l] = x
      Y[:, l] = y
  for i, inclination in enumerate(incl):  
      Y[:, i] = Y[:, i]*np.cos(inclination) #projection
      #rotating by PA, projecting
      X[:, i], Y[:, i] = rotateGalaxy(X[:, i], Y[:, i], pa)  
      data = (X[:, i], Y[:, i], np.zeros((X.shape)), np.zeros((Y.shape)), 5) #5 -- dummy R_e parameter
      params = (vc, c, gamma, pa, inclination)
      model_vel[:, i] = model2_Courteau(params, data)
      model_vel[:, i] = model_vel[:, i] 
  simple_plot(X, Y, model_vel, 'img/models/test.png')
  return X, Y, model_vel, 10*np.ones((model_vel.shape))	

