import matplotlib.pyplot as plt
from astropy.table import Table, Column
import numpy as np
import importlib
import sys
import time

#sub packages
from setup import prefix, Folder, HPMS_eDR3_file, GCNS_cat_file, \
	GCNS_reject_file, HPMS_limit, save_table_process, make_plots, form_out
import plt_setup as ps
ff = []

def re():
	import good_HPMS as GH
	importlib.reload(GH)

def broadcast():
	global prefix, Folder, HPMS_eDR3_file, GCNS_cat_file, \
		GCNS_reject_file, HPMS_limit, save_table_process, make_plots, form_out
	from setup import prefix, Folder, HPMS_eDR3_file, GCNS_cat_file, \
		GCNS_reject_file, HPMS_limit, save_table_process, make_plots, form_out

def load_Data():
	print('HPMS load Data')
	global HPMS
	HPMS = Table.read(Folder + 'Data/' + HPMS_eDR3_file)
def load_GCNS():
	print('HPMS load GCNS')
	global GCNS_cat, GCNS_reject, GCNS_reject_good, GCNS_reject_bad
	GCNS_cat = Table.read(Folder + 'Data/' + GCNS_cat_file)
	GCNS_reject = Table.read(Folder + 'Data/' + GCNS_reject_file)
	GCNS_reject_good = GCNS_reject[np.where(GCNS_reject['GCNS_PROB'] > 0.38)]
	GCNS_reject_bad = GCNS_reject[np.where(GCNS_reject['GCNS_PROB'] < 0.38)]


def match_GCNS():
	print('HPMS match GCNS')
	global index_GCNS_cat, index_GCNS_reject_good, index_GCNS_reject_bad
	global HPMS_in_GCNS_cat, HPMS_in_GCNS_reject_good, HPMS_in_GCNS_reject_bad
	index_GCNS_cat = np.isin(HPMS['source_id'], GCNS_cat['SOURCE_ID'])
	index_GCNS_reject_good = \
			np.isin(HPMS['source_id'], GCNS_reject_good['SOURCE_ID'])
	index_GCNS_reject_bad = \
			np.isin(HPMS['source_id'], GCNS_reject_bad['SOURCE_ID'])
	HPMS_in_GCNS_cat = HPMS[index_GCNS_cat]
	HPMS_in_GCNS_reject_good = HPMS[index_GCNS_reject_good]
	HPMS_in_GCNS_reject_bad = HPMS[index_GCNS_reject_bad]
	

def HPMS_check_parallax(limit = HPMS_limit['px']):
	# check if parallax significance
	out = HPMS['parallax_over_error'] > limit
	print('HPMS px:', np.sum(out), '/', len(out))
	return out


def HPMS_check_ruwe(limit = HPMS_limit['ruwe']):
	# check if ruwe below limit
	out = HPMS['ruwe'] < limit
	print('HPMS Ruwe:', np.sum(out), '/', len(out))
	if make_plots:
		
		fig = plt.figure("ruwe")
		fig.clf()
		ax = plt.axes()
		ps.semilogy(HPMS['phot_g_mean_mag'],\
			HPMS['ruwe'], color =ps.red, \
			ms = 0.5, label = "excluded HPMS")
		ps.semilogy(HPMS[out]['phot_g_mean_mag'], \
			HPMS[out]['ruwe'], color =ps.blue,\
			ms = 0.5, label = "HPMS")
		xlim=plt.xlim()
		plt.plot(xlim,[2,2], color = 'orange', label = 'ruwe = 2')
		# plt.xlim(xlim)
		plt.ylabel("Ruwe")
		plt.xlabel("G [mag]")
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		ax.yaxis.set_major_formatter(ps.formatter)

		plt.pause(0.05)
	return out

	return out

def HPMS_n_obs_sig_g_flux(power = HPMS_limit['n_obs_sig_g_flux_power'],\
	limit = HPMS_limit['n_obs_sig_g_flux']):
	 

	# check if ruwe below limit
	a = HPMS['phot_g_mean_flux_over_error'].data * 1e0
	b =  np.power(HPMS['phot_g_n_obs'].data *1e0, power)
	out =a*b> limit 
	f2 = out == False
	print('HPMS NvsF:', np.sum(out), '/', len(out))
	if make_plots:
		fig = plt.figure("n_vs_sig_g_flux")
		fig.clf()
		ax = plt.axes()
		ps.loglog(HPMS[f2]['phot_g_n_obs'],\
			HPMS[f2]['phot_g_mean_flux_over_error'], color =ps.red, \
			ms = 0.5, label = 'excluded HPMS')
		ps.loglog(HPMS[out]['phot_g_n_obs'], \
			HPMS[out]['phot_g_mean_flux_over_error'], color =ps.blue,\
			ms = 0.5, label = "HPMS")
		plt.loglog(HPMS_in_GCNS_reject_bad['phot_g_n_obs'], \
			HPMS_in_GCNS_reject_bad['phot_g_mean_flux_over_error'], "s", \
			color ="k" , ms = 2, label = "GCNS probability < 0.38")
		plt.ylabel(r"$G_{flux}/\sigma_{G_{flux}}$")
		plt.xlabel(r"$n_{obs}$")
		xlim = plt.xlim()
		ylim = plt.ylim()

		s1 = "%i"%power if power%1==0 else "%.1f"%power
		if limit > 1000: 
			kk = "%.2e"%(limit)
			if kk[2:4]=="00":
				s2 = '%s \cdot 10^{%i}'%(kk[0],int(kk.split('e')[-1]))
			else: 
				s2 = '%s \cdot 10^{%i}'%(kk[0:4],int(kk.split('e')[-1]))
		else:
			s2 = "%i"%limit 
		

		plt.plot([np.power(limit/ylim[0],1/power),np.power(limit/ylim[1],1/power)],ylim,\
			color = ps.limitcolor,label = \
			r"$G_{flux}/\sigma_{G_{flux}}\cdot n_{obs}^{%s} < %s$"%(s1,s2)) 
		# plt.xlim(xlim)
		# plt.ylim(ylim)

		leg = plt.legend()
		for lh in leg.legendHandles:  
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		ax.xaxis.set_major_formatter(ps.formatter)
		ax.yaxis.set_major_formatter(ps.formatter)
		plt.pause(0.05)
	return out


def HPMS_check_GCNS():
	# check if GNCS_prob > 0.38	
	out = index_GCNS_cat | index_GCNS_reject_good
	print('HPMS GCNS good:', np.sum(out), '/', len(out),\
		np.sum(index_GCNS_cat))
	return out


def HPMS_check_GCNS_bad():
	# check if GNCS_prob < 0.38	
	out = index_GCNS_reject_bad
	print('HPMS GCNS bad:', np.sum(out), '/', len(out))
	return out

def HPMS_check_phot(cat = None,limit = HPMS_limit['mag']):
	# check if phot_G < 21	
	global HPMS
	if cat is None: cat = HPMS
	out = cat['phot_g_mean_mag'] < limit
	print('HPMS photG:', np.sum(out), '/', len(out))
	return out

def main():
	broadcast()
	tt = []
	tt.append(time.time())
	if 'HPMS' not in globals():
		load_Data()
	tt.append(time.time())
	if 'GCNS_cat' not in globals():
		load_GCNS()
	tt.append(time.time())
	if 'index_GCNS_cat' not in globals():
		match_GCNS()
	tt.append(time.time())
	px = HPMS_check_parallax()
	rw = HPMS_check_ruwe()
	ns_good = HPMS_check_GCNS()
	ns_bad = HPMS_check_GCNS_bad()
	phot = HPMS_check_phot()
	NvsF = HPMS_n_obs_sig_g_flux()
	good = (px & rw) & (ns_bad == False) & phot & NvsF
	print('HPMS good:', np.sum(good), '/', len(HPMS))
	tt.append(time.time())
	if save_table_process:
		print('save good HPMS')
		print(Folder + 'Results/HPMS.good' + prefix +form_out[0])
		HPMS[good].write(Folder + 'Results/HPMS.good'+prefix + form_out[0],\
			format = form_out[1], overwrite=True)
		HPMS_bad = HPMS[good==False]
		HPMS_bad["out_px"]= Column(px[good==False] == False)
		HPMS_bad["out_ruwe"]= Column(rw[good==False] == False)
		HPMS_bad["out_phot"]= Column(phot[good==False] == False)
		HPMS_bad["out_N_vs_SigFlux"]= Column(NvsF[good==False] == False)
		HPMS_bad["out_ns_bad"]= Column(ns_bad[good==False])
		HPMS_bad.write(Folder + 'Results/HPMS.bad' \
			+ prefix + form_out[0], format = form_out[1], overwrite = True)
	tt.append(time.time())
	tt = np.array(tt)
	cpt = tt[1:]-tt[:-1]
	#print('HPMS:', *(cpt))


	if make_plots:
		fig = plt.figure("HPMS")
		fig.clf()
		pmtot = np.sqrt(HPMS['pmra'] * HPMS['pmra'] \
			+ HPMS['pmdec']*HPMS['pmdec'])
		ax = plt.subplot(111)


		plt.loglog(pmtot[good==False],\
		 	HPMS[good==False]['parallax'],'o', color =ps.red, ms = 0.5, label = 'excluded HPMS')
		ps.loglog(pmtot[good], \
			HPMS[good]['parallax'], color =ps.blue,	ms = 0.5, label = "HPMS")
		plt.xlabel(r"$\mu_{tot}$ [mas/y]")
		plt.ylabel(r"$\varpi$ [mas]")
		leg = plt.legend()
		for lh in leg.legendHandles: 
		 	lh._legmarker.set_alpha(1)
		 	lh._legmarker.set_markersize(5)
		plt.pause(0.05)
	return HPMS[good]['source_id']









