import math
import numpy as np
from scipy import stats


def fix_geometry(incl, pa, vc):
	incl = wrapAngle(incl)
	incl, vc = foldIncl(incl, vc)
	pa = wrapAngle(pa)
	pa, vc = foldPa(pa, vc)
	pa_std = np.std(pa)
	circularMeanPA = getCircularMean(pa)
	circularMeanIncl = getCircularMean(incl)
	if np.mean(vc) >=0:
		folded_vc = np.mean(np.abs(vc), axis=None)
	else:
		folded_vc = -1*np.mean(np.abs(vc), axis=None)
	return circularMeanIncl, circularMeanPA, folded_vc


def getInclVec(ba):
  q = 0.2 
  ret = np.zeros((ba.shape))
  ret[np.where(ba <= 0.2)] = math.pi
  ret[np.where(ba > 0.2)] = np.arccos(np.sqrt((ba[np.where(ba > 0.2)]**2 - q**2)/(1 - q**2)))
  
  return ret


def getIncl(ba):
  q = 0.2 
  if ba < 0.2:
    return math.pi/2
  else:
    return np.arccos(np.sqrt((ba**2 - q**2)/(1 - q**2)))

def getPhi(y, x):
   return np.arctan2(y, x)


def rotateGalaxy(x, y, ang):
   x_rot = (x)*np.cos(ang) - (y)*np.sin(ang)
   y_rot = (x)*np.sin(ang) + (y)*np.cos(ang)
   return (x_rot, y_rot)

def foldPa(pa, vc):
    vc[np.where(pa > math.pi)] *=-1
    pa[np.where(pa >= math.pi)] -= math.pi
    pa[np.where(pa <=-math.pi)] += math.pi
    pa[np.where(pa < 0)] += math.pi 
    if np.any(np.abs(pa) >= math.pi):
		pa, vc = foldPa(pa, vc)    
    return pa, vc	

def wrapAngle(ang):
  ang = np.mod(ang, 2*math.pi)
  return ang

def reshape_array(array):
	return np.reshape(array, (array.shape[0], 1))

def getCircularMean(angles):
  n = len(angles)
  sineMean = np.divide(np.sum(np.sin(angles)), n)
  cosineMean = np.divide(np.sum(np.cos(angles)), n)
  vectorMean = math.atan2(sineMean, cosineMean)
  print n, np.degrees(sineMean), np.degrees(cosineMean), np.degrees(vectorMean)
  return vectorMean

def getQuadrant(angle):
	ret = -1*np.ones((angle.shape[0], ), dtype=int)
	ret[(angle >= 0) & (angle < math.pi/2)] = 1
	ret[(angle >= math.pi/2) & (angle < math.pi)] = 2
	ret[(angle >= math.pi) & (angle < math.pi+math.pi/2)] = 3
	ret[(angle >= math.pi+math.pi/2) & (angle < 2*math.pi)] = 4
	if (ret == -1).any():
	 print	'wrong angle value!'
	 exit()
	return ret

def foldIncl(incl, vc):
	if type(incl) is np.float64:
	  incl = np.asarray([incl])
	  vc = np.asarray([vc])
	  conv = "1"
	else:
	  conv = "0"
	Quad = getQuadrant(incl)
	#2nd quadrant
	incl[np.where(Quad == 2)] = -1*(incl[np.where(Quad == 2)] - math.pi) 
	#vc[np.where(Quad == 2)] *= -1
	
	#3rd quadrant
	incl[np.where(Quad == 3)] -= math.pi #TODO: check if vc should be flipped
	vc[np.where(Quad == 3)] *= -1
	
	#4th quadrant
	incl[np.where(Quad == 4)] = -1 * (incl[np.where(Quad == 4)] - 2*math.pi)
	vc[np.where(Quad == 4)] *= -1
	if conv == "1":
	  incl = incl[0]
	  vc = vc[0]
   	return incl, vc


#print pa
#pa = wrapAngle(pa)
#print pa
#pa = foldPa(pa, vc)
#print pa
#print vc
#print getQuadrant(np.asarray([0.2, math.pi/2+0.1, math.pi+0.1, 2*math.pi-0.1]))

#incl = np.asarray([(np.radians(1)), np.radians(120.), np.radians(220), np.radians(330)])
#vc = np.asarray([300, 300, 300, 300])
#print np.degrees(incl)
#incl, vc = foldIncl(incl, vc)
#print np.degrees(incl)
#print vc

#print np.degrees(getCircularMean([0, np.pi/2, -np.pi/2]))
