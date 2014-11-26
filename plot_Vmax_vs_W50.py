from __future__ import division
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import math
from utils import *
from nicer_plots import *
import db


def plot_Vmax_vs_W50():
	sami_ids = db.dbUtils.getFromDB('sami_id', 'db/SAMI.sqlite', 'mcmc ')
	Vmax = np.abs(db.dbUtils.getFromDB('Vmax', 'db/SAMI.sqlite', 'mcmc '))
	Vmax_err = db.dbUtils.getFromDB('Vmax_err', 'db/SAMI.sqlite', 'mcmc ')
	W50 = np.abs(db.dbUtils.getFromDB('HI_Vc', 'db/SAMI.sqlite', 'mcmc '))
	W50_err = np.abs(db.dbUtils.getFromDB('HI_Vc_err', 'db/SAMI.sqlite', 'mcmc '))
	fig = plt.figure()
	plt.scatter(W50, Vmax, c='k')
	plt.plot([0, 300], [0, 300])
	plt.axis([0, 400, 0, 400])
	plt.ylabel("Vmax, km/s")
	plt.title(str(len(W50))+" galaxies")
	plt.xlabel("Inclination-corrected W50/2")
	plt.errorbar(W50, Vmax, linestyle='none', xerr=W50_err, yerr=Vmax_err, c='k', zorder=-100)
	plt.savefig("img/Vmax_vs_HI_Vc")


def plot_Vmax_vs_W50_prior():
	sami_ids = db.dbUtils.getFromDB('sami_id', 'db/SAMI.sqlite', 'mcmc_prior ')
	Vmax = np.abs(db.dbUtils.getFromDB('Vmax', 'db/SAMI.sqlite', 'mcmc_prior '))
	Vmax_err = db.dbUtils.getFromDB('Vmax_err', 'db/SAMI.sqlite', 'mcmc_prior ')
	W50 = np.abs(db.dbUtils.getFromDB('HI_Vc', 'db/SAMI.sqlite', 'mcmc_prior '))
	W50_err = np.abs(db.dbUtils.getFromDB('HI_Vc_err', 'db/SAMI.sqlite', 'mcmc_prior '))
	fig = plt.figure()
	plt.scatter(W50, Vmax, c='k')
	plt.plot([0, 300], [0, 300])
	plt.axis([0, 400, 0, 400])
	plt.ylabel("Vmax, km/s")
	plt.title(str(len(W50))+" galaxies")
	plt.xlabel("Inclination-corrected W50/2")
	plt.errorbar(W50, Vmax, linestyle='none', xerr=W50_err, yerr=Vmax_err, c='k', zorder=-100)
	plt.savefig("img/Vmax_vs_HI_Vc_prior")
	
	
plot_Vmax_vs_W50_prior()	
plot_Vmax_vs_W50()	
