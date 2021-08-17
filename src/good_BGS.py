
from astropy.table import Table,unique, Column
import numpy as np

import sys
import os
import time

#sub packages
import plot_functions as pf

from setup import Folder, BGS_eDR3_file, DR2_BGS_file, BGS_limit,\
	DR2_limit, zeropoint,save_table_process,make_plots, form_out, prefix,\
	random_sample_file, dr2_random_file
from utils import cosdeg




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
	gamma = np.maximum(pow(10, 0.2 * (BGS['phot_g_mean_mag'] - 18)), 1)
	BGS['psi'] = BGS['astrometric_sigma5d_max'] / (1.2 * gamma)
	print('BGS:', len(BGS))

def BGS_ruwe(limit = BGS_limit['ruwe']):
	#check if ruwe below limit
	#only for source with a five parameter solution

	# only two parameter solution
	two_parm = np.isnan(BGS['ruwe']) | (BGS['ruwe'] > 1e10)  

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
		pf.plot_gof(BGS,out,limit)
	print('BGS GoF:', np.sum(out), '/', len(out))
	return out

def BGS_px(limit = BGS_limit['px']):
	# check if parallax is not signigicant negative
	#only for source with a five parameter solution
	two_parm = np.isnan(BGS['parallax']) | (BGS['parallax'] > 1e10)
	good_px = (BGS['parallax'] > limit * BGS['parallax_error'] + zeropoint) 
	print('BGS px:', np.sum(good_px& (two_parm==False)), '/', 
		len(good_px)-np.sum(two_parm))
	out = two_parm | good_px
	if make_plots: 
		pf.plot_px(BGS[two_parm==False],out[two_parm==False])
	return out

def BGS_mc_Gill():
	# check if PSI value is above 1
	# only for sources with G < 18 
	# see P. McGill et al. 2020
	# not used to filter data, since most are true Gaia_eDR3 sources 
	phot = BGS['phot_g_mean_mag'] > 18
	psi = BGS['psi'] 
	good_psi = psi < 1
	out = good_psi | phot
	print('BGS PSI:', np.sum(out) - np.sum(phot), '/', len(out) - np.sum(phot))
	if make_plots: 
		if "random_sample" not in globals():
			if os.path.isfile(Folder + 'Data/' + random_sample_file):
				print("BGS load random_sample")
				global random_sample
				random_sample = \
					Table.read(Folder + 'Data/' + random_sample_file)
			else: sys.exit()
		if "random_sample" in globals():
			pf.plot_psi_result_part_1(random_sample)
			pf.plot_psi(BGS,psi,out,random_sample)
		else: pf.plot_psi(BGS,psi,out, None)

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
		if "dr2_random" not in globals():
			if os.path.isfile(Folder + 'Data/' + dr2_random_file):
				print("BGS load random_sample")
				global dr2_random
				dr2_random = \
					Table.read(Folder + 'Data/' + dr2_random_file)
		if 'dr2_random' in globals():			
			pf.plot_DR2_match(DR2_BGS,good,excluded,dist_limit, dr2_random)
		else: 
			pf.plot_DR2_match(DR2_BGS,good,excluded,dist_limit,None)

	print('BGS in DR2:')
	print('\t-found matches:', len(DR2_BGS))
	print('\t-good matches: ', len(good_DR2_BGS))
	print('\t-bad matches:  ',len(bad_DR2_BGS))
	print('\t-miss matches: ',np.sum(good_dist == False))

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

	five_parm = (np.isnan(dr3['parallax']) == False) \
			& (dr3['parallax'] < 1e10)
	two_par = np.isnan(dr3['parallax']) | (dr3['parallax'] < 1e10)
	print(np.sum(good) + np.sum(bad)+np.sum(miss_match), len(DR2_BGS))

	print('BGS proper motion good in DR2:',np.sum(good), '/', len(DR2_BGS))
	print('BGS proper motion bad in DR2:',np.sum(bad), '/', len(DR2_BGS))
	if make_plots: 
		pf.plot_DR2_pm(dr3[nn],dr2,pm[qq])
	
	pm[miss_match] = [0.,0.]
	pm_dict = dict(zip(BGS['source_id'],pm))
	return out, good, bad, pm_dict


def BGS_pos_error(limit = BGS_limit['pos_err']):
	# Positional error better than 100
	F_pos = BGS['ra_error'] * BGS['ra_error'] \
		+ BGS['dec_error'] * BGS['dec_error'] < limit * limit

	print('BGS pos:',np.sum(F_pos),'/', len(BGS))
	
	if make_plots:
		pf.plot_pos_err(BGS,F_pos,limit)
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
		BGS_BAD["dr2_out"]= Column(badDR2[good==False])
		BGS_BAD["out_pos"]= Column(good_pos[good==False] == False)

		BGS_BAD.write(Folder + 'Results/BGS.bad' \
			+ prefix + form_out[0], format = form_out[1], overwrite = True)
	tt.append(time.time())
	tt = np.array(tt)
	cpt = tt[1:] - tt[:-1]
	#print('BGS:', *(cpt))
	return good_source_ID, bad_DR2, pmDR2



