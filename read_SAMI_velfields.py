from __future__ import division
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import math
from utils import *
from nicer_plots import *
import pyfits
import db
from mcmc import *
from geom import *
from test_disk import *
import sys
import itertools
from scipy import stats
from astroML.plotting import hist
from integrals import get_final_model_j


#MCMC settings:
Nwalkers = 400 # number of walkers
Nburn = 200 # number of burn-in samples
Nmcmc = 400 # number of MCMC samples
Nthreads = 4 # number of threads to run
fit = 'exp'
test=False
prior = sys.argv[2]
priorType = sys.argv[3]
baryonType = sys.argv[4]



def plotPDF(params, names, filename):
	fig = plt.figure(figsize=(10, 10))
	nSubplots = int(math.factorial(len(params))/(math.factorial(2) * math.factorial(len(params) - 2)))
	print nSubplots, len(params)
	indices = itertools.combinations(range(0, len(params)), 2)
	for i, ind in enumerate(indices):
		print i, ind
		ax = plt.subplot(nSubplots, 1, i)
		plt.scatter(params[ind[0]], params[ind[1]], s=1, edgecolor='none', c='k')
		plt.xlabel(names[ind[0]])
		plt.ylabel(names[ind[1]])
	plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
	plt.savefig(filename, bbox_inches='tight') 
	

		
def fit_galaxy(name, x, y, vel, vel_err, r50, GAMA_incl, HI_linewidth, HI_Vc_err, initParams, fit, test):	
	#get data
	data = (name, x, y, vel, vel_err, r50, GAMA_incl, HI_linewidth, HI_Vc_err)
	Npar=len(initParams)
	if fit == 'model2':
		lnprob0 = lnProb(initParams, data)
		p0 = [getRandParams(initParams) for i in xrange(Nwalkers)]
		sampl = ensemble.EnsembleSampler(Nwalkers, Npar, lnProb, args=[data], threads=Nthreads)
	elif fit == 'exp':
		lnprob0 = lnProbExp(initParams, data)
		p0 = [getRandParams(initParams) for i in xrange(Nwalkers)]
		sampl = ensemble.EnsembleSampler(Nwalkers, Npar, lnProbExp, args=[data], threads=Nthreads)
	
	pos, prob, state = sampl.run_mcmc(p0, Nburn)
	sampl.reset()
	sampl.run_mcmc(pos, Nmcmc, rstate0=state)
	acceptanceFraction = np.mean(sampl.acceptance_fraction)
	
	if fit == 'model2':
		vc = sampl.flatchain[:,0]
		c = sampl.flatchain[:,1]
		gamma = sampl.flatchain[:,2]
		pa = sampl.flatchain[:,3]
		incl = sampl.flatchain[:,4]
		incl, pa, vc = fix_geometry(incl, pa, vc)
		vc_mod, c_mod, gamma_mod, pa_mod, incl_mod, v0_mod = (vc, stats.mode(c), stats.mode(gamma), pa, incl, stats.mode(v0))
		print 'vc', vc_mod, 'c', c_mod, 'gamma', gamma_mod, 'pa', np.degrees(pa_mod), 'incl', np.degrees(incl_mod), 'v0', v0_mod
		modelParams = (vc_mod, c_mod, gamma_mod, pa_mod, incl_mod)
		modelVelField = model2_Courteau(modelParams, data)      
		model_radius = np.linspace(0, 25, 100)
		
		model_RC = (model_radius*np.abs(vc_mod))/((c_mod**2 + model_radius**2)**(gamma_mod/2))

	elif fit == 'exp':
		if prior =='True':
			PDFfilename = 'img/pdf/'+priorType+'_prior_'+name+baryonType
			mcmcOutFile = priorType+'_prior_mcmc.csv'
			mcmcImgName = 'img/models/good/'+priorType+baryonType+'_prior_'+fit+'_RC_'+name
			histName = 'img/models/hist/'+priorType+baryonType+'_prior_incl_'+name
			incl_chainFilename = 'chains/'+priorType+baryonType+'_prior_incl_'+name
			vc_chainFilename = 'chains/'+priorType+baryonType+'_prior_vc_'+name
			c_chainFilename = 'chains/'+priorType+baryonType+'_prior_c_'+name
		else:
			PDFfilename = 'img/pdf/'+name+baryonType
			mcmcOutFile = 'mcmc'+baryonType+'.csv'
			mcmcImgName = 'img/models/good/'+baryonType+fit+'_RC_'+name
			histName = 'img/models/hist/incl_'+baryonType+name
			incl_chainFilename = 'chains/incl_'+baryonType+name		
			vc_chainFilename = 'chains/vc_'+baryonType+name
			c_chainFilename = 'chains/c_'+baryonType+name		
		
		vc = sampl.flatchain[:,0]
		c = sampl.flatchain[:,1]
		pa = sampl.flatchain[:,2]
		incl = sampl.flatchain[:,3]
		v0 = sampl.flatchain[:,4]
		chain_incl = foldIncl(wrapAngle(incl), vc)[0]
		np.savetxt(incl_chainFilename, chain_incl)
		np.savetxt(vc_chainFilename, vc)
		np.savetxt(c_chainFilename, c)
		
		plot_hist(np.degrees(chain_incl), histName)

		incl_mod, pa_mod, vc_mod = fix_geometry(incl, pa, vc)
		c_mod, v0_mod = (np.mean(c), np.mean(v0))
		print 'vc', vc_mod, 'c', c_mod, 'pa', np.degrees(pa_mod), 'incl', np.degrees(incl_mod), 'v0', v0_mod		
		modelParams = (vc_mod, c_mod, pa_mod, [incl_mod], v0_mod)
		modelVelField = expModel(modelParams, (x, y, r50))      
		model_radius = np.linspace(0, 25, 1000)
		model_RC = vc_mod*(1 - np.exp(-1*(model_radius/(c_mod*r50)))) #((2/math.pi) * vc_mod* np.tanh(model_radius/(c_mod*r50))) 
		model_linewidth = HI_linewidth/(2*np.sin(GAMA_incl))
		modelV_rot = vc_mod*(1 - np.exp(-1*(80/(c_mod*r50))))
		
	delta_coords = get_delta_coords(name)
	delta_z = get_delta_z(name)
	plotPDF((vc, c, np.degrees(foldPa(pa, vc)[0]), np.degrees(foldIncl(wrapAngle(incl), vc)[0])), ('Vmax', 'c', 'pa', 'incl'), PDFfilename) #check number of parameters -- adjustable plot
	fit_radius, fit_vel = getRotCurveFromVelField((pa_mod, incl_mod, v0_mod), data)	  
	fig = plt.figure()
	fig.add_subplot(221)
	#plt.title("$\Delta$ ra, $\Delta$ dec: "+str(delta_coords)+" $\Delta$z: "+str(delta_z)	)
	plt.title("SAMI "+str(name)+", v.08 kinematics")
	plt.axhline(c='k')
	plt.axvline(c='k')
	cd = plt.scatter(x, y, c=vel, marker='s', edgecolor='none', vmin=-150, vmax=150)
	plt.colorbar(cd, label='vel, km/s')


	fig.add_subplot(222)
	#plt.title("Fit rotation curve")
	plt.title(' Incl:' +str(round(np.degrees(incl_mod), 1))+' GAMA incl:' +str(round(np.degrees(GAMA_incl), 1)))
	cd = plt.scatter(model_radius, np.abs(model_RC), c='k', s=4, edgecolor='none', label='Model RC', zorder=100)
	d = plt.scatter(fit_radius, np.abs(fit_vel), c=np.abs(vel), s=5, edgecolor='none', label='Fit RC')
	plt.axhline(model_linewidth, c='r', linewidth=2, linestyle='dotted', label='W50/(2sin(i))')	
	plt.errorbar(fit_radius, np.abs(fit_vel), yerr=vel_err, c='k', linestyle='none', zorder=-200) 
	plt.axis([0, 30, 0, 230])
	plt.legend(loc='best')
	fig.add_subplot(223)
	cd = plt.scatter(x, y, c=modelVelField, marker='s', edgecolor='none', vmin=-150, vmax=150)
	plt.colorbar(cd, label='Model velocity, km/s')
	plt.axhline(c='k')
	plt.axvline(c='k')
	plt.title('Vc:' + str(np.round(modelV_rot, 1))+'km/s, HI linewidth: '+str(round(model_linewidth))+" km/s, W50: "+str(round(HI_linewidth)))
	fig.add_subplot(224)
	cd = plt.scatter(x, y, c= vel - modelVelField, edgecolor='none', marker='s', vmin=-150, vmax=150)
	plt.colorbar(cd, label='Residuals, km/s')
	plt.title("Residuals")
	plt.axhline(c='k')
	plt.axvline(c='k')

	plt.savefig(mcmcImgName)
	plt.close()
	
	#f = open(mcmcOutFile, 'a')
	#f.write(str(filename[:-15])+", "+str(vc_mod)+", "+str(np.std(sampl.flatchain[:,0]))+", "+str(c_mod)+", "+str(v0_mod)+", "+str(pa_mod)+", "+str(incl_mod)+", "+str(model_linewidth)+", "+str(HI_Vc_err)+"\n")
	#f.close()
	
	
	if prior == 'True':
		#j PDF figure
		#reading no prior fit values from previous chains
		c_no_prior = np.genfromtxt('chains/c_'+name)
		vc_no_prior = np.genfromtxt('chains/vc_'+name)
		
		j_no_prior = get_final_model_j(r50*c_no_prior, vc_no_prior)
		
		j = get_final_model_j(r50*c, vc)

		n, bins, patches = hist(j, bins='scott')
		n_no_prior, bins_no_prior, patches = hist(j_no_prior, bins='scott')
		
		fig = plt.figure()	
		ax_Pxy = plt.axes((0.2, 0.34, 0.27, 0.52))
		ax_Py = plt.axes((0.474, 0.34, 0.1, 0.52))


		plt.title("Fit rotation curve and j PDF")
		cd = ax_Pxy.scatter(model_radius, np.abs(model_RC), c='k', s=4, edgecolor='none', label='Model RC', zorder=100)
		d = ax_Pxy.scatter(fit_radius, np.abs(fit_vel), c=np.abs(vel), s=5, edgecolor='none', label='Fit RC')
		ax_Pxy.axhline(model_linewidth, c='r', linewidth=2, linestyle='dotted', label='W50/(2sin(i))')	
		ax_Pxy.errorbar(fit_radius, np.abs(fit_vel), yerr=vel_err, c='k', linestyle='none', zorder=-200) 
		ax_Pxy.axis([0, 30, 0, 230])
		plt.xlabel("radius, arcsec")
		plt.ylabel(r"$v_{rot}$, km/s")
		ax_Pxy.legend(loc='best')

		bin_centers = (bins[:-1] + bins[1:])/2.
		bin_centers_no_prior = (bins_no_prior[:-1] + bins_no_prior[1:])/2.
		
		ax_Py.plot(n, bin_centers/np.max(bin_centers), c='r', label='j PDF, HI/inclination prior')
		ax_Py.plot(n_no_prior, bin_centers_no_prior/np.max(bin_centers), c='teal', label='j PDF, no prior')
		plt.gca().xaxis.set_major_locator(plt.NullLocator())
		plt.gca().yaxis.set_major_locator(plt.NullLocator())
		ax_Py.set_ylim([-1.03,2])
		ax_Py.legend(loc='best')
		ax_Py.patch.set_visible(False)
		ax_Py.axis('off')
		plt.savefig(mcmcImgName+"_"+baryonType+"_marginal_j")

if test:
	name='test'
	testParams = 150, 0.2, np.radians(30), np.radians([60])
	x, y, vel, vel_err = make_exp_test_disk(10, 10, testParams)
	r50 = 5 #test
	#init fit
	vc_init = 120
	c_init = 0.5
	gamma_init = 1 #corresponds to arctan fit, see Courteau 1997
	pa_init = np.radians(40)
	incl_init = np.radians(75)
	if fit == 'exp':
		#Vmax, c,  pa, incl
		initParams = [vc_init, c_init, pa_init, incl_init]
	elif fit == 'model2':
		#Vc, c, gamma,  pa, incl
		initParams = [150, 1, 1.0, np.radians(45), np.radians(45)]
	fit_galaxy(name, x, y, vel, vel_err, r50, initParams, fit, test)	


else:
	if sys.argv[1] == 'all':
		names = db.dbUtils.getFromDB('sami_id', 'db/SAMI.sqlite', 'ALFALFA_Xmatch ')
	else:
		names = [sys.argv[1]]
	for name in names:
		print name, 'NAME'
		#filename = 'data/'+name+'_1_comp.fits.gz'
		filename = 'data/stars/'+str(name)+'_kinematics.fits'
		try:
			#x, y, vel, vel_err = get_gas_velfield(filename)
			x, y, vel, vel_err = get_stellar_velfield(filename)
			simple_plot(x, y, vel, 'img/test/stars/'+name)
			r50, W50, W50_err = get_SAMI_data(name)
			print W50, W50_err
			GAMA_incl = get_GAMA_incl(name)
		except IOError:
			print 'No such file!'
			continue
		if fit == 'exp':
			print 'FITTING ', filename[:-15]
			#Vmax, c,  pa, incl
			initParams = [150, 2.2, np.radians(45), GAMA_incl, -5]
		elif fit == 'model2':
			#Vc, c, gamma,  pa, incl
			initParams = [150, 1, 1.0, np.radians(45), GAMA_incl]
		fit_galaxy(name, x, y, vel, vel_err, r50, GAMA_incl, W50, W50_err, initParams, fit, test)
