from __future__ import division
import numpy as np
from astropy.coordinates.distances import Distance
import matplotlib.pyplot as plt
import pyfits
import db
from string import *
from astroML.plotting import hist


def simple_plot(x, y, vel, filename):
  fig = plt.figure(figsize=(10, 10))
  ax = fig.add_subplot(111, aspect='equal',autoscale_on=False, xlim=[-25,25], ylim=[-25,25])
  cb = ax.scatter(x, y, c=vel, edgecolor="none", vmin=-120, vmax=120)
  plt.colorbar(cb)
  ax.axhline(c='k')
  ax.axvline(c='k')
  plt.savefig(filename)

def plot_hist(x, filename):
  fig = plt.figure(figsize=(10, 10))
  hist(x, bins='blocks')
  plt.savefig(filename)
  plt.close()	

def get_SAMI_data(sami_id):
	r50 = db.dbUtils.getFromDB('R_e', 'db/SAMI.sqlite', 'SAMI_Master ', ' where sami_id = '+ str(sami_id))[0]
	W50 = db.dbUtils.getFromDB('W50', 'db/SAMI.sqlite', 'ALFALFA_Xmatch ', ' where sami_id = '+ str(sami_id))[0]
	W50_err = db.dbUtils.getFromDB('W50_err', 'db/SAMI.sqlite', 'ALFALFA_Xmatch ', ' where sami_id = '+ str(sami_id))[0]
	return r50, W50, W50_err

def get_SAMI_coords(sami_id):
	ra = db.dbUtils.getFromDB('ra', 'db/SAMI.sqlite', 'SAMI_Master ', ' where sami_id = '+ str(sami_id))[0]
	dec = db.dbUtils.getFromDB('dec', 'db/SAMI.sqlite', 'SAMI_Master ', ' where sami_id = '+ str(sami_id))[0]
	return ra, dec

def get_delta_z(sami_id):
	z = db.dbUtils.getFromDB('z', 'db/SAMI.sqlite', 'SAMI_Master ', ' where sami_id = '+ str(sami_id))[0]
	alfalfa_id = db.dbUtils.getFromDB('ALFALFA_id', 'db/SAMI.sqlite', 'ALFALFA_Xmatch ', ' where sami_id = '+ str(sami_id))[0]
	alfalfa_z = db.dbUtils.getFromDB('V_P', 'db/SAMI.sqlite', 'ALFALFA ', ' where object ='+"'"+lstrip(str(alfalfa_id))+"'")/300000
	return round(np.abs(float(z) - float(alfalfa_z)), 6)


def	get_delta_coords(sami_id):
	ra, dec = get_SAMI_coords(sami_id)
	alfalfa_id = db.dbUtils.getFromDB('ALFALFA_id', 'db/SAMI.sqlite', 'ALFALFA_Xmatch ', ' where sami_id = '+ str(sami_id))[0]
	print 'where object ='+"'"+lstrip(str(alfalfa_id))+"'"
	alfalfa_ra = 15*db.dbUtils.getFromDB('raopt', 'db/SAMI.sqlite', 'ALFALFA ', ' where object ='+"'"+lstrip(str(alfalfa_id))+"'")
	alfalfa_dec = db.dbUtils.getFromDB('decopt', 'db/SAMI.sqlite', 'ALFALFA ', ' where object ='+"'"+lstrip(str(alfalfa_id))+"'")
	return (round(np.abs(float(alfalfa_ra) - float(ra)), 6), round(np.abs(float(alfalfa_dec) - float(dec)), 6))
	
def get_ALFALFA_W50(ra, dec):
	alfalfa_ra = db.dbUtils.getFromDB('raopt', 'db/SAMI.sqlite', 'ALFALFA')
	alfalfa_dec = db.dbUtils.getFromDB('decopt', 'db/SAMI.sqlite', 'ALFALFA')
	#print np.round(alfalfa_ra*15, 1), np.round(alfalfa_dec, 1)
	W50 = db.dbUtils.getFromDB('W50', 'db/SAMI.sqlite', 'ALFALFA ', ' where round(15*raopt, 1) = '+ str(np.round(ra, 1))+' and round(decopt, 1) = '+str(np.round(dec, 1)))[0]
	return W50



def get_velfield(filename):
		print filename
		all_vel = pyfits.getdata(filename, extname='V', header=False)[1]
		all_vel_err = pyfits.getdata(filename, extname='V_ERR', header=False)[1]
		#good = np.where(all_vel_err) < 100
		#vel = all_vel[good]	
		#vel_err = all_vel_err[good]
		all_vel = np.ma.masked_invalid(all_vel)
		all_vel_err = np.ma.masked_invalid(all_vel_err)
		mask = np.where(all_vel_err < 20)
		
		#getting indices, i.e. y and x:
		ind = np.column_stack(mask) - 25
		x, y = np.asarray(zip(*ind))	
		vel_err = all_vel_err.filled()[mask]
		vel = all_vel.filled()[mask]
		#print 'HI', HI_linewidth
		return x, y, vel, vel_err

def angular2physical(arcsec, z): #return physical effective diameter of the galaxy in kpc
    return (np.radians(arcsec/3600) *Distance(z=z).kpc / (1 + z)**2)


def sqlify(arr):
  strings = ''
  for i in arr:
     if type(i) == type(tuple()):
        i = i[0]   
     strings = strings+","+'"'+strip(str(i))+'"'
  strings = '('+strings[1:]+')'
  return strings
  
  
def convert_pc_to_meters(pc):
 return pc*3.0857*10e16

def decodeU(query_output):
  output = []
  for u in query_output:
    u = str(u)
    output.append(u)
  return output



def get_ALFALFA_data():

	ra = db.dbUtils.getFromDB('ra', 'db/SAMI.sqlite', 'SAMI_Master ')
	dec = db.dbUtils.getFromDB('dec', 'db/SAMI.sqlite', 'SAMI_Master ')
	SAMI_all_ids = db.dbUtils.getFromDB('sami_id', 'db/SAMI.sqlite', 'SAMI_Master ')
	raopt = db.dbUtils.getFromDB('raopt', 'db/SAMI.sqlite', 'ALFALFA')
	decopt = db.dbUtils.getFromDB('decopt', 'db/SAMI.sqlite', 'ALFALFA')
	#print np.round(alfalfa_ra*15, 1), np.round(alfalfa_dec, 1)
	ALFALFA_ids = []
	sami_ids = []
	for sami_ra, sami_dec, sami_id in zip(ra, dec, SAMI_all_ids):
		obj = (db.dbUtils.getFromDB('Object', 'db/SAMI.sqlite', 'ALFALFA ', ' where round(15*raopt, 1) = '+ str(np.round(sami_ra, 1))+' and round(decopt, 1) = '+str(np.round(sami_dec, 1))))
		if len(obj) == 1:
			obj = decodeU(obj)[0]
			ALFALFA_ids.append(obj)
			sami_ids.append(sami_id)
	W50 = db.dbUtils.getFromDB('W50', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	W50_err = db.dbUtils.getFromDB('Werr', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	int_flux = db.dbUtils.getFromDB('sintmap', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	ra =  db.dbUtils.getFromDB('raopt', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	dec =  db.dbUtils.getFromDB('decopt', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	SN =  db.dbUtils.getFromDB('SN', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	rms =  db.dbUtils.getFromDB('rms', 'db/SAMI.sqlite', 'ALFALFA ', ' where Object in'+sqlify(ALFALFA_ids))
	f = open('db/ALFALFA_Xmatch.csv', 'a')
	for i, s in enumerate(W50): #not galaxies have flux measurements	
		print i, s
		f.write(str(sami_ids[i])+", "+ str(ALFALFA_ids[i])+", "+str(ra[i])+", "+str(dec[i])+", "+str(W50[i])+", "+str(W50_err[i])+", "+str(int_flux[i])+", "+str(SN[i])+", "+str(rms[i])+"\n")
	f.close()	

