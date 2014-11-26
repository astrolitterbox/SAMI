from __future__ import division
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import math
from utils import *
from integrals import *
from scipy.odr import Data, Model, ODR, RealData, odr_stop
from nicer_plots import *


filenames = np.genfromtxt("Simona/THINGS_files.txt", dtype='object')
J_HLR = []
for filename in filenames:
	name = filename[:-9]
	data = np.genfromtxt("Simona/"+filename, delimiter = ",")
	mask = np.ones(data.shape[0], dtype=bool)
	mask[np.where(data[:, 3] == 0)] = False
	data = data[mask]
	radius = data[:, 0]
	Sigma = data[:, 1]
	Sigma_err = data[:, 2]# + 1/Sigma
	vel = data[:, 3]
	vel_err = 10 + data[:, 4] * 1/Sigma
	#Half-light radius from data
	hlr = np.where(np.round(np.cumsum(Sigma)/np.max(np.cumsum(Sigma)), 1) > 0.5)[0][0]
	print hlr, 'hlr', np.round(np.cumsum(Sigma)/np.max(np.cumsum(Sigma)), 1)
	HLR = radius[hlr]
	#removing velocity range above 1 HLR:
	
	vel2 = np.ma.masked_where((radius > HLR) & (radius < radius[-1]) , vel)
	print vel.shape
	vel_err2 = np.ma.masked_where((radius > HLR) & (radius < radius[-1]) , vel_err)
	radius2 = np.ma.masked_where((radius > HLR) & (radius < radius[-1]) , radius)
	print radius2
	vel_err2[-1] = 40
	#vel2[-1] = Vmax

	#Fitting rotation curve with ODR:
	RC = Model(RC_model)
	mydata = RealData(radius2, vel2, sy=vel_err2/np.std(vel2))
	myodr = ODR(mydata, RC, beta0=[150, 10000])
	myoutput = myodr.run()
	myoutput.pprint()
	Vmax, Rflat = myoutput.beta
	stopReason = myoutput.stopreason
	stopInt = myoutput.info
	
	#Fitting Sigma with ODR:
	exponential = Model(Sigma_model_2exp)
	mydata = RealData(radius, Sigma, sy=Sigma_err/np.std(Sigma))
	myodr = ODR(mydata, exponential, beta0=[10000, 500, 100, 1000])
	myoutput = myodr.run()
	myoutput.pprint()
	S1, S2, R1, R2 = myoutput.beta
	
	
	#J calculation from data 
	J = get_J_in_ring(radius, Sigma, vel)
	ring_mass = get_M_in_ring(radius, Sigma)
	j = get_j(J, ring_mass) 

	#J from integrals
	J_int, err_int = integrate_J_in_ring(radius, S1, S2, R1, R2, Vmax, Rflat)	
	
	J_HLR.append(np.cumsum(J)[hlr]/np.max(np.cumsum(J)))
		
	#j integral
	M_integral = get_M_integral(radius, S1, S2, R1, R2)
	j_int = get_j_int(J_int, M_integral)
	
	fig = plt.figure(figsize=(12, 10))
	ratio_string = str(np.round(R2/Rflat, 2))+r" $R_{s2}/R_{flat}$"
	if stopInt <> 1:
		plt.title(stopReason[0], color='r')
	else:
		plt.title(stopReason[0], color='k', fontsize=11)
	ax = fig.add_subplot(611)
	ax.plot(radius, RC_model([Vmax, Rflat], radius), c='r', label='Fit')
	#ax.errorbar(radius[radius < hlr], vel[radius < hlr], yerr=vel_err[radius < hlr], c='k', linestyle='none')
	ax.axvline(Rflat, c='r', label="$R_{flat}$")
	ax.axvline(HLR, c='k', label="$R_{0.5}$")
	ax.plot(radius, vel, c='k')
	ax.axvspan(0, HLR, facecolor='g', alpha=0.4)
	ax.axvspan(np.max(radius)-1, np.max(radius), facecolor='g', alpha=0.4)
	ax.axis([0, radius[-1]+300, 0, np.max(vel)+20])
	plt.xlabel("radius, pc")
	plt.ylabel("velocity, km/s")

	ax = fig.add_subplot(612)
	ax.errorbar(radius, Sigma, yerr=Sigma_err, c='k', linestyle='none')
	ax.plot(radius, Sigma, c='k')
	ax.plot(radius, Sigma_model_2exp([S1, S2, R1, R2], radius), c='r', label='Fit')
	ax.axvline(HLR, c='k', label="$R_{0.5}$")
	plt.xlabel("radius, pc")
	plt.ylabel(r"$\Sigma, M_{\odot}/pc^2$")
	plt.yscale('log')


	ax = fig.add_subplot(613)
	ax.plot(radius, J_int, c='r', label='Int')
	ax.errorbar(radius, J_int, yerr=err_int, c='r', linestyle='none')
	ax.plot(radius, J, c='k', label='Data')
	ax.axvline(HLR, c='k', label="$R_{0.5}$")
	ax.legend(loc='best')
	plt.xlabel("radius, pc")
	plt.ylabel(r" $J_{ring}$")
	

	ax = fig.add_subplot(614)
	ax.plot(radius, np.cumsum(J), c='k', label='Data')
	ax.plot(radius, np.cumsum(J_int), c='r', label='Int')
	ax.axvline(HLR, c='k', label="$R_{0.5}$")
	ax.legend(loc='best')
	plt.xlabel("radius, pc")
	plt.ylabel(r"Cumulative $J$")
	
	ax = fig.add_subplot(615)
	ax.plot(radius, np.cumsum(M_integral), c='r', label='Int')
	ax.plot(radius, np.cumsum(ring_mass), c='k', label='Data')
	ax.axvline(HLR, c='k', label="$R_{0.5}$")
	ax.legend(loc='best')
	plt.xlabel("radius, pc")
	plt.ylabel(r"Cumulative $M$")
	
	ax = fig.add_subplot(616)
	ax.plot(radius, np.cumsum(j_int), c='r', label='Int')
	ax.plot(radius, np.cumsum(j), c='k', label='Data, '+ratio_string)
	ax.axvline(HLR, c='k', label="$R_{0.5}$= "+str(round(HLR/1000, 1)))
	ax.legend(loc='best')
	plt.xlabel("radius, pc")
	plt.ylabel(r"Cumulative $j$")	
	plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	plt.savefig("img/test_"+name, bbox_inches='tight')

	#f = open('fit_values.csv', 'a')
	#f.write(name+", "+str(odr_S)+", "+str(odr_R_s)+", "+str(Vmax)+", "+str(Rflat)+"\n")
	#f.close	


