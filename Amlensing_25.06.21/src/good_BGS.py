
#import matplotlib.pyplot as plt
from astropy.table import Table,unique, Column
import numpy as np
import importlib
import sys
import os
import matplotlib.pyplot as plt
import time

#sub packages
import plt_setup as ps
from setup import Folder, BGS_eDR3_file, DR2_BGS_file, BGS_limit,\
	DR2_limit, zeropoint,save_table_process,make_plots, form_out, prefix,\
	random_sample_file, dr2_random_file
from utils import cosdeg

plt.ion()
ff = [None,None,None,None,None,None]

def re():

	import good_BGS as GB
	importlib.reload(ps)
	importlib.reload(GB)

def broadcast():
	global Folder, BGS_eDR3_file, DR2_BGS_file, BGS_limit,\
		DR2_limit, zeropoint,save_table_process,make_plots, form_out, prefix,\
		random_sample_file, dr2_random_file
	from setup import Folder, BGS_eDR3_file, DR2_BGS_file, BGS_limit,\
		DR2_limit, zeropoint,save_table_process,make_plots, form_out, prefix,\
		random_sample_file, dr2_random_file

#for development


def BGS_load_Data(limit = BGS_limit['mag']):
	# load BGS eDR3 table
	print('BGS load Data')
	global BGS
	BGS = Table.read(Folder + 'Data/' + BGS_eDR3_file)
	BGS = BGS[BGS['phot_g_mean_mag'] < limit]
	print('BGS:', len(BGS))

def BGS_ruwe(limit = BGS_limit['ruwe']):
	#check if ruwe below limit
	#only for source with a five parameter solution
	two_parm = np.isnan(BGS['ruwe']) # only two parameter solution
	good_ruwe = BGS['ruwe'] < limit
	print('BGS Ruwe:', np.sum(good_ruwe), '/', len(BGS) - np.sum(two_parm))
	out = two_parm | good_ruwe
	return out

def BGS_gof(limit = BGS_limit['ruwe']):
	# check if gof/n_good_obs below corresponding ruwe limit
	# aplicable also for sources with only a two parameter solution
	
	# translate ruwe limit to GoF limit
	f = lambda x: -4.80459 + 0.000520143 * np.sqrt(4.964e7 * x + 3.57727e7)
	#f(2) = 1.24; f(1.4) = 0.532
	
	out = BGS['astrometric_gof_al'] \
		/ np.sqrt(BGS['astrometric_n_good_obs_al']) < f(limit)

	if make_plots: 
		fig = plt.figure("gof")
		fig.clf()
		ps.loglog(BGS['astrometric_n_good_obs_al'],\
			BGS['astrometric_gof_al'],color = ps.red, \
			ms = 0.5 , label = "excluded sources" )
		ps.loglog(BGS['astrometric_n_good_obs_al'][out],\
			BGS['astrometric_gof_al'][out],color = ps.blue, ms = 0.5, \
			label = "Sou_GoF/sqrt(Sou_N) < %.2f"%f(limit))
		good_ruwe = BGS['ruwe'] < limit
		ps.loglog(BGS['astrometric_n_good_obs_al'][good_ruwe],\
			BGS['astrometric_gof_al'][good_ruwe],color = ps.green, \
			ms = 0.5,  label = "RUWE < %.1f"%limit)		
		xlim = np.array(plt.xlim())
		plt.loglog(xlim, np.sqrt(xlim) * f(limit),color = ps.limitcolor,\
			linewidth = 1)
		# plt.xlim(xlim)
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		
		plt.ylabel('Sou_GoF')
		plt.xlabel('Sou_N')
		plt.pause(0.1)
		ff[0] = fig
	print('BGS GoF:', np.sum(out), '/', len(out))
	return out

def BGS_px(limit = BGS_limit['px']):
	# check if parallax is not signigicant negative
	#only for source with a five parameter solution
	two_parm = np.isnan(BGS['parallax'])
	good_px = BGS['parallax'] > limit * BGS['parallax_error'] + zeropoint
	print('BGS px:', np.sum(good_px), '/', len(good_px)-np.sum(two_parm))
	out = two_parm | good_px
	if make_plots: 
		fig = plt.figure("px")
		fig.clf()
		_,b,_ = plt.hist(BGS['parallax'], color= ps.red, bins = 1000,rwidth = 0.8)
		plt.hist(BGS[out]['parallax'],  color= ps.blue, bins = b, rwidth = 0.5)
		#plt.yscale("log")
		plt.xlim([-5,5])
		plt.ylabel('#')
		plt.xlabel(r'$SOU\_\varpi$')
		plt.pause(0.1)
		ff[1] = fig
	return out

def BGS_mc_Gill():
	# check if PSI value is above 1
	# only for sources with G < 18 
	# see P. McGill et al. 2020
	# not used to filter data, since most are true Gaia_eDR3 sources 
	phot = BGS['phot_g_mean_mag'] > 18
	gamma = np.maximum(pow(10, 0.2 * (BGS['phot_g_mean_mag'] - 18)), 1)
	psi = BGS['astrometric_sigma5d_max'] / (1.2 * gamma)
	good_psi = psi < 1
	out = good_psi | phot
	print('mc_Gill:', np.sum(out) - np.sum(phot), '/', len(out) - np.sum(phot))
	if make_plots: 
		if "random_sample" not in globals():
			if os.path.isfile(Folder + 'Data/' + random_sample_file):
				print("BGS load random_sample")
				global random_sample
				random_sample = \
					Table.read(Folder + 'Data/' + random_sample_file)

		fig = plt.figure("psi")
		fig.clf()

		ps.semilogy(BGS['phot_g_mean_mag'][out], psi[out],\
			color = ps.blue,ms = 0.5,zorder = 1, label = "All BGS")
		ps.semilogy(BGS['phot_g_mean_mag'][out == False], psi[out == False], \
			color = ps.red, ms = 0.5, label = r"$\Psi > 1$ & G < 18 mag", \
			zorder = 0)
		if "random_sample" in globals():
			rgamma = np.maximum(pow(10, \
				0.2 * (random_sample['phot_g_mean_mag'] - 18)), 1)
			rpsi = random_sample['astrometric_sigma5d_max'] / (1.2 * rgamma)
			ps.semilogy(random_sample['phot_g_mean_mag'], rpsi,\
				color = ps.grey, ms = 0.05,  zorder = 3, \
				label = "random sample")
			c = np.arange(5.5,22,0.1)
			d = np.array([np.percentile(rpsi[np.abs(random_sample['phot_g_mean_mag']-i) < 1],90) for i in c])
			d2 = np.array([np.percentile(BGS[np.abs(BGS['phot_g_mean_mag']-i) < 1]["psi"],90) for i in c])
			# plt.plot(c,d2, color = ps.blue, label ="90th percentile BGS",\
			#zorder = 3)
			plt.plot(c,d, color = ps.limitcolor, label ="90th percentile",\
				zorder = 100)
		plt.ylabel(r"$\Psi$")
		plt.xlabel(r"G [mag]")
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		
		if "random_sample" in globals():
			fig = plt.figure("psi_result")
			fig.clf()
			rgamma = np.maximum(pow(10, \
				0.2 * (random_sample['phot_g_mean_mag'] - 18)), 1)
			rpsi = random_sample['astrometric_sigma5d_max'] / (1.2 * rgamma)
			ps.semilogy(random_sample['phot_g_mean_mag'], rpsi,\
				color = ps.grey, ms = 0.05,  zorder = 0, \
				label = "random sample")
			c = np.arange(5.5,22,0.1)
			d = np.array([np.percentile(rpsi[np.abs(random_sample['phot_g_mean_mag']-i) < 1],90) for i in c])
			d2 = np.array([np.percentile(BGS[np.abs(BGS['phot_g_mean_mag']-i) < 1]["psi"],90) for i in c])
			# plt.plot(c,d2, color = ps.blue, label ="90th percentile BGS",\
			#zorder = 3)
			plt.plot(c,d, color = ps.limitcolor, label ="90th percentile",\
				zorder = 100)
			plt.ylabel(r"$\Psi$")
			plt.xlabel(r"G [mag]")
	
		plt.pause(0.1)

	return out 

def BGS_load_dr2(dist_limit = DR2_limit['dist'], mag_limit = DR2_limit['mag']):
	# load Gaia best DR2 neighbour 
	# validate if it is a true match
	print('BGS load Data DR2')
	global DR2_BGS, good_DR2_BGS, bad_DR2_BGS	
	DR2_BGS = Table.read(Folder + 'Data/' + DR2_BGS_file)
	DR2_BGS.sort('angular_distance')
	DR2_BGS = unique(DR2_BGS, keys = 'dr3_source_id')

	good_mag = DR2_BGS['magnitude_difference']**2 < mag_limit**2 \
		* DR2_BGS['angular_distance']**0.4
	good_dist = DR2_BGS['angular_distance'] < dist_limit
	good = good_mag	& good_dist	
	good_DR2_BGS = DR2_BGS[good]
	excluded = (good_mag == False) & good_dist
	bad_DR2_BGS = DR2_BGS[excluded]
	if make_plots:
		fig = plt.figure("DR2")
		fig.clf()
		if "dr2_random" not in globals():
			if os.path.isfile(Folder + 'Data/' + dr2_random_file):
				print("BGS load random_sample")
				global random_sample
				dr2_random = \
					Table.read(Folder + 'Data/' + dr2_random_file)
		if "dr2_random" not in globals():
			ps.semilogx(dr2_random["angular_distance"], \
				np.abs(dr2_random["magnitude_difference"]), ms = 0.1,\
				color = ps.grey,label = "random sample")
		ps.loglog(DR2_BGS[good]["angular_distance"], \
			np.abs(DR2_BGS[good]["magnitude_difference"]), ms = 0.5, \
			color = ps.blue,alpha = 0.5, label = "good matches")
		ps.loglog(DR2_BGS[good == False]["angular_distance"], \
			np.abs(DR2_BGS[good == False]["magnitude_difference"]), \
			ms = 0.5, color = ps.green, label = "mismatches")
		ps.loglog(DR2_BGS[excluded]["angular_distance"], \
			np.abs(DR2_BGS[excluded]["magnitude_difference"]), ms = 0.5, \
			color = ps.red, label = "excluded sources")
		xlim = np.array([plt.xlim()[0]*5,dist_limit])
		f =lambda x: 0.3 * pow(x,0.2)
		plt.loglog(xlim, f(xlim),ps.limitcolor, linewidth = 1)
		plt.loglog([400,400],[2e-4,10],ps.limitcolor, linewidth = 1)

		plt.ylim([1e-4,30])
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		plt.xlabel("$\Delta\phi$ [mas]")
		plt.ylabel("$\Delta G$ [mag]")

		plt.pause(0.1)

	print('BGS in DR2:', len(DR2_BGS), len(good_DR2_BGS), len(bad_DR2_BGS))

def BGS_dr2_dr3_propermotion(pm_limit = BGS_limit['pm'], \
	pm_limit_bad = DR2_limit['pm_bad']): 
	if 'BGS' not in globals():
		BGS_load_Data()
	if 'good_DR2_BGS' not in globals():
		BGS_load_dr2()
	dr2 = good_DR2_BGS[np.isin(good_DR2_BGS['dr3_source_id'],BGS['source_id'])]
	dr3 = BGS[np.isin(BGS['source_id'], good_DR2_BGS['dr3_source_id'])]
	vv = np.isin(BGS['source_id'], good_DR2_BGS['dr3_source_id'])
	dr2.sort('dr3_source_id')
	nn = np.argsort(dr3['source_id'])
	out = np.array([False] * len(BGS))
	pm = np.zeros([len(BGS),2])
	sid = np.zeros(len(BGS))
	qq = np.where(vv)[0][nn]
	out[qq] = True
	T = 0.5
	pm[qq,0] = (dr3[nn]['ra'] - dr2['ra']) * cosdeg(dr2['dec']) / T * 3.6e6
	pm[qq,1] = (dr3[nn]['dec'] - dr2['dec']) / T * 3.6e6


	good = pm[:,0]**2 + pm[:,1]**2 < pm_limit**2
	bad = (pm[:,0]**2 + pm[:,1]**2 < pm_limit_bad**2) \
		& (pm[:,0]**2 + pm[:,1]**2 > pm_limit**2)
	miss_match = (pm[:,0]**2 + pm[:,1]**2 > pm_limit_bad**2)

	five_parm = np.isnan(dr3[nn]['parallax']) == False
	two_par = np.isnan(dr3[nn]['parallax'])
	print('BGS good in DR2:',np.sum(good), '/', len(DR2_BGS))
	print('BGS bad in DR2:',np.sum(bad), '/', len(DR2_BGS))
	if make_plots: 
		fig = plt.figure("DR2_pm", figsize = [6.4,3])
		fig.clf()
		fig.subplots_adjust(right = 0.87, bottom = 0.17)
		ax1 = plt.subplot(121)
		d3 = dr3[nn]
		five = np.isnan(d3['parallax']) == False
		fivedr2 = np.isnan(dr2['parallax']) == False
		twodr2 = np.isnan(dr2['parallax'])

		pp = pm[qq]

		ps.plot(d3[five & twodr2]["pmra"], pp[five & twodr2][:,0], \
			color = ps.blue, label = "2par in dr2")
		ps.plot(d3[fivedr2]["pmra"], pp[fivedr2][:,0], color = ps.green,\
			label = "5par in dr2")


		ax2 = plt.subplot(122)
		ps.plot(d3[five & twodr2]["pmdec"], pp[five & twodr2][:,1], \
			color = ps.blue ,label = "2par in dr2")
		ps.plot(d3[fivedr2]["pmdec"], pp[fivedr2][:,1], color = ps.green,\
			label = "5par in dr2")
		ax2.yaxis.tick_right()
		ax2.yaxis.set_label_position("right")
		ax1.set_xlabel(r'$\mu_{\alpha^{\star},DR3}$ [mas/yr]')
		ax1.set_ylabel(r'$\Delta \phi_{\rm \times2,\,\alpha^{\star}}$ [mas]')
		ax2.set_xlabel(r'$\mu_{\delta,DR3}$ [mas/yr]')
		ax2.set_ylabel(r'$\Delta \phi_{\rm \times2,\,\delta}$ [mas]')
		ax1.set_xlim([-1700,1700])
		ax1.set_ylim([-1700,1700])
		ax2.set_xlim([-1700,1700])
		ax2.set_ylim([-1700,1700])

		leg = ax1.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		# leg = ax2.legend()
		# for lh in leg.legendHandles: 
			# lh._legmarker.set_alpha(1)
			# lh._legmarker.set_markersize(5)
		#plt.tight_layout(pad=0.05)
		plt.pause(0.1)

		ff[0] = fig
	
	pm[miss_match] = [0.,0.]
	pm_dict = dict(zip(BGS['source_id'],pm))
	return out, good, bad, pm_dict


def BGS_pos_error(limit = BGS_limit['pos_err']):
	# Positional error better than 100
	F_pos = BGS['ra_error'] * BGS['ra_error'] \
		+ BGS['dec_error'] * BGS['dec_error'] < limit * limit
	
	if make_plots: 
		fig = plt.figure("pos_err")
		fig.clf()
		ax = plt.axes()
		ps.semilogy(BGS['phot_g_mean_mag'][F_pos==False], \
			np.sqrt(BGS['ra_error'][F_pos==False] \
			* BGS['ra_error'][F_pos==False] \
			+ BGS['dec_error'][F_pos==False] \
			* BGS['dec_error'][F_pos==False]), ms = 0.5, color = ps.red, \
			label = "excluded BGS" )
		ps.semilogy(BGS['phot_g_mean_mag'][F_pos], \
			np.sqrt(BGS['ra_error'][F_pos] * BGS['ra_error'][F_pos] + \
			BGS['dec_error'][F_pos] * BGS['dec_error'][F_pos]), ms = 0.5, \
			color = ps.blue, label = \
			r"$\sqrt{\sigma^{2}_{\alpha} + \sigma^{2}_{\delta}} < $" \
			+ "%i mas"%limit)
		xlim=plt.xlim()
		plt.plot(xlim,[10,10], color = ps.limitcolor, label = r'$\sigma_{\rm pos} = 10$',\
			zorder = 5)
		# plt.xlim(xlim)
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		plt.xlabel("G [mag]")
		plt.ylabel(r"$\sqrt{\sigma_{\alpha}^{2} + \sigma_{\delta}^{2}}$ [mas]")
		ax.yaxis.set_major_formatter(ps.formatter)
		plt.pause(0.01)
		ff[4] = fig
	return F_pos


def main():
	tt = []
	tt.append(time.time())
	if 'BGS' not in globals():
		BGS_load_Data()
	if make_plots:
		BGS_load_dr2()
	tt.append(time.time())
	good_ruwe = BGS_ruwe()
	good_gof = BGS_gof()
	good_px	 = BGS_px()
	good_pos = BGS_pos_error()
	_ = BGS_mc_Gill()

	tt.append(time.time())
	inDR2,goodDR2,bad_DR2, pmDR2 = BGS_dr2_dr3_propermotion()
	badDR2 = np.isin(BGS["source_id"], bad_DR2_BGS['dr3_source_id'])
	tt.append(time.time())
	good = good_ruwe & good_gof & good_px & good_pos & (badDR2 == False)
	good_source_ID = BGS['source_id'][good]

	tt.append(time.time())
	print('BGS good:',np.sum(good),'/', len(BGS))
	if save_table_process:
		print('save good BGS')
		print(Folder + 'Results/BGS.good' +prefix +  form_out[0])
		BGS[good].write(Folder+'Results/BGS.good' + prefix + form_out[0], \
			format = form_out[1], overwrite = True)
		print('save bad DR2 matches')
		print(Folder + 'Results/BGS.bad' +prefix+ form_out[0])
		BGS_BAD = BGS[good==False]
		BGS_BAD["out_ruwe"]= Column(good_ruwe[good==False] == False)
		BGS_BAD["out_gof"]= Column(good_gof[good==False] == False)
		BGS_BAD["out_px"]= Column(good_px[good==False] == False)
		BGS_BAD["dr2_out"]= Column(badDR2[good==False] == False)
		BGS_BAD["out_pos"]= Column(good_pos[good==False] == False)

		BGS_BAD.write(Folder + 'Results/BGS.bad' \
			+ prefix + form_out[0], format = form_out[1], overwrite = True)
	tt.append(time.time())
	tt = np.array(tt)
	cpt = tt[1:] - tt[:-1]
	#print('BGS:', *(cpt))
	return good_source_ID, bad_DR2, pmDR2



