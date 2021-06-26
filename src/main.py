import numpy as np # version >= 1.19.5
import matplotlib.pyplot as plt # version >= 3.3.3
import os
import sys
import time
from astropy.table import vstack, Table , MaskedColumn # version >= 4.2

import raw_data
import good_BGS
import mass
import microlensing
import calculation 
import find_closest 
import setup
import plt_setup as ps

Folder = setup.Folder
print(Folder)

def broadcast():
	raw_data.broadcast()
def determine_mass(cand):
	# calculate Mass (approximate_mass.py)
	#mass,mass_err,stellar_type = calculate_Mass(good_raw_cands)
	print('calc assumed mass ')
	cMass = mass.approx_Mass(cand)
	cand['star_type'] = cMass[2]
	cand['mass'] = cMass[0]
	cand['mass_error'] = cMass[1]
	return cand

def determine_einstein_radii(cand):

	#replace masked values
	# calculate Einstein Radii (see microlensing.py)
	cEinstein = microlensing.calculate_Einsteinradius(cand)
	cand['ThetaE'] = cEinstein[0]
	cand['ThetaE_error'] = cEinstein[1]
	print('-------------------------')
	return cand
def approx(cand):
	# calculate_approx_date (see find_closest.py) 
	print('calc approx shift ')

	cApprox = find_closest.calc_approx_date(cand)
	cand['approx_tca'] = cApprox[0]
	cand['approx_dist'] = cApprox[1]
	cand['t_aml'] = MaskedColumn(2/np.sqrt((cand['pmra']-cand['ob_pmra'])**2 \
		+(cand['pmdec']-cand['ob_pmdec'])**2)\
		*np.sqrt(cand['ThetaE'] ** 4 /0.1**2+cand['approx_dist']*2),\
		dtype = 'float64',unit = 'years', \
		description = 'Expected duration of the event (shift > ~0.1mas')
	# calculate approximated effect (see microlensing.py)
	approx_shift = microlensing.calculate_Effect(cand, approx = 'True')
	if setup.save_table_process: # save step to file
		prefix = setup.prefix
		form_out = setup.form_out
		print('save approx')
		cand['approx_shift'] = approx_shift
		print(Folder + 'Results/amlensing.approx'+ prefix + form_out[0])
		cand.write(Folder + 'Results/amlensing.approx' + prefix+ \
			form_out[0], format = form_out[1], overwrite = True)
	# filter n approximated effect
	cand = cand[(approx_shift > 0.03) \
		| (cand['approx_dist'] < 2 * cand['parallax'])]
	cand = cand[(cand['approx_tca'] > 2010) \
		& (cand['approx_tca'] < 2070)]

	print('-------------------------')
	return cand

def calc(cand, n_core = 4):
	# split on multiple cores
	if n_core > len(cand): n_core = len(cand)
	if n_core > 1: 
		from joblib import Parallel, delayed # version >= 1.0.0
		import pickle	
		print('split table onto cores')

		#split raw_data into multiple sub folder for parallel computing
		t = (time.time() * 100) % 10000
		instring = ['%s.TMPc%ip%i.pkl' % (Folder, t, core) \
			for core in range(n_core)]
		cands_per_core = len(cand) / n_core	
		for core in range(n_core):
			f = open(instring[core], 'wb')
			raw_cands_parallel = cand[round(cands_per_core * core) \
					: round(cands_per_core * (core + 1))]
			pickle.dump(raw_cands_parallel, f)
			f.close()	

		#parallel computing 
		#(see calculation.py, forwars to microlensing.py and find_closest.py)
		res_par = Parallel(n_jobs = n_core)(delayed(calculation.parallel)(i) \
				for i in instring)

		# get results
		table_out=[]
		for string in instring:
			f = open(string, 'rb') # load temporary files
			table_out.append(pickle.load(f))
			f.close()
			os.remove(string) # delete temporary files
		print('stack table')
		table_out = vstack(table_out)

	else:
		#single computing 
		#(see calculation.py, forwars to microlensing.py and find_closest.py)
		table_out = calculation.single(cand)
	print('-------------------------')
	return(table_out)

def filter_events(table_out):	
	print('filter events')
	prefix = setup.prefix
	form_out = setup.form_out

	if setup.save_table_process: # save step to file
		

		print('save unfiltered results')
		# save all events
		print(Folder + 'Results/amlensing.unfiltered'+ prefix + form_out[0])
		table_out.write(Folder + 'Results/amlensing.unfiltered' + prefix \
			+ form_out[0], format = form_out[1], overwrite = True)
	# filter results 
	shift_01 = table_out['shift_plus'] > 0.1
	shift_L2_01 = table_out['L2_shift_plus'] > 0.1
	mag = table_out['magnification'] > 0.001
	FF = (mag | shift_01) | shift_L2_01


	# compare with Blacklist

	if os.path.isfile(setup.Blacklist_file) \
		or os.path.isfile(Folder + 'Data/'+setup.Blacklist_file):
		
		# load blacklist
		if os.path.isfile(setup.Blacklist_file):
			BL = Table.read(setup.Blacklist_file, data_start = 0)
		elif os.path.isfile(Folder + 'Data/'+setup.Blacklist_file) :
			BL = Table.read(Folder + 'Data/'+setup.Blacklist_file,\
			 data_start = 0)
		BL1 = []
		BL2 = []
		psi = 0
		for i in BL[BL.colnames[0]]:
			print(i)
			if i[0] == '#':
				if i[:5] == '# NO ':
					psi = 1
			elif psi: BL2.append(int(i))
			else: BL1.append(int(i))
		BL1 = np.isin(table_out['ob_source_id'],BL1) \
			| np.isin(table_out['source_id'],BL1)
		BL2 = np.isin(table_out['ob_source_id'],BL2) \
			| np.isin(table_out['source_id'],BL2)

		aa = np.sum(FF)
		bb = np.sum(FF & (BL1 == False))
		cc = np.sum(FF & (BL2 == False))
		print(aa-bb)
		print(aa-cc)

		#does blacklist include ob_source_id?
		BL = [int(i) for i in BL[BL.colnames[0]] if i[0] != '#']
		blacklist = np.isin(table_out['ob_source_id'],BL) \
			| np.isin(table_out['source_id'],BL)

		if setup.save_table_process: # save step to file
			print(Folder + 'Results/amlensing' \
			+ prefix + ".no_blacklist" + form_out[0])
			table_out[FF].write(Folder + 'Results/amlensing' \
			+ prefix + ".no_blacklist" + form_out[0],\
			format = form_out[1], overwrite = True)
		FF = FF & (blacklist == False)
		if setup.make_plots:
			if 'psi_result' in plt.get_figlabels():
				plt.figure('psi_result')
				BGS = good_BGS.BGS
				bl = np.isin(BGS['source_id'],\
					table_out[blacklist]['ob_source_id'])
				print(1)
				bl2= bl & (BGS['dec'] < -30)

				res = np.isin(BGS['source_id'],table_out[FF]['ob_source_id'])
				g = BGS['phot_g_mean_mag'][res]
				g_BL = BGS['phot_g_mean_mag'][bl]
				g_BL2 = BGS['phot_g_mean_mag'][bl2]

				if 'psi' not in BGS.colnames:				
					psi = BGS['psi'][res]
					psi_BL = BGS['psi'][bl]
					psi_BL2 = BGS['psi'][bl2]
				else:
					gamma = np.maximum(pow(10, 0.2 * (BGS['phot_g_mean_mag'] - 18)), 1)
					psi_bgs = BGS['astrometric_sigma5d_max'] / (1.2 * gamma)
					psi = psi_bgs[res]
					psi_BL = psi_bgs[bl]
					psi_BL2 = psi_bgs[bl2]
				plt.plot(g,psi,'o',ms = 1, color = ps.blue,\
					label = r'predicted events', zorder = 1)	
				plt.plot(g_BL,psi_BL,'s',ms = 1, color = ps.red, zorder = 2,\
					label = r'removed candidates  $(\delta >-30^{\circ})$')		
				plt.plot(g_BL2,psi_BL2,'s',ms = 1, color = 'm', zorder = 3,\
					label = r'removed candidates  $(\delta <-30^{\circ})$')
				leg = plt.legend()
				for lh in leg.legendHandles: 
					lh._legmarker.set_alpha(1)
					lh._legmarker.set_markersize(5)
			else: print('No plot found for "psi_result"')

	else: 
		print("No Blacklist file found")
	# save results 
	print('save filtered results')
	print(Folder + 'Results/amlensing' + prefix +form_out[0])
	if form_out[1] == 'fits':
		for j,i in enumerate(table_out.keys().copy()):
			table_out.meta['TCOMM%i'%(j+1)] = table_out[i].description
	table_out[FF].write(Folder + 'Results/amlensing' + prefix +form_out[0],\
		format = form_out[1], overwrite = True)
	print('-------------------------')

	return table_out[FF]

def plot_results(result):
	
	for i in ['All','WD','RG','MS','BD']:
		for j in [[2010,2070],[2014,2032]]:
			for k in [20,6,3]:
				if k == 20 : add =''
				else: add = '_DG_'+str(k)
				if i == 'All':
					if j[1] == 2070:
						fig = plt.figure('Results_'+i+add,figsize= [12.8,4.8])
						plt.subplots_adjust(\
							left = plt.rcParams['figure.subplot.left']/2,
							right= plt.rcParams['figure.subplot.right']/2+0.5)
					else: fig = plt.figure('Results_2030_'+i+add)
					which = np.ones(len(result),bool)
				elif i == 'MS':
					if j[1] == 2070:
						fig = plt.figure('Results_'+i+add, figsize= [12.8,4.8])
						plt.subplots_adjust(\
						left = plt.rcParams['figure.subplot.left']/2,
						right= plt.rcParams['figure.subplot.right']/2+0.5)
					else:fig = plt.figure('Results_2030_'+i+add)
					which = result['star_type']==i
				else:
					if j[1] == 2070: fig = plt.figure('Results_'+i+add)
					else: fig = plt.figure('Results_2030_'+i+add)
					which = result['star_type']==i
				which = which & ((result['ob_phot_g_mean_mag']-result['phot_g_mean_mag']) < k)
				ax = plt.gca()
				twopar = which & (result['ob_parallax'] == 0)
				fivepar = which & (result['ob_parallax'] != 0)
				plt.semilogy(result[twopar]['TCA'],result[twopar]['shift_plus'],'.', color = ps.grey, zorder = -5,
					label = '2 parameter solution')
				plt.errorbar(result[fivepar]['TCA'],result[fivepar]['shift_plus'],fmt = 'o',\
					xerr = result[fivepar]['TCA_error'],
					yerr = [-result[fivepar]['shift_plus_error_m'],\
					result[fivepar]['shift_plus_error_p']],\
					color = 'blue', ms = 2,lw = 0.5, markeredgewidth = 0.2,\
					markeredgecolor = 'k', label = '5 parameter solution')
				plt.ylim([0.05,50])
				plt.yticks([0.1,0.2,0.5,1,2,5,10,20,50],\
					['0.10','0.20','0.50','1.0','2.0','5.0','10','20','50'])
				plt.xlim(j)

				leg = plt.legend()
				for lh in leg.legendHandles: 
					try: 
						lh._legmarker.set_alpha(1)
						lh._legmarker.set_markersize(5)
					except: pass
				plt.ylabel(r'$\delta\theta_{+}$ [mas]')
				plt.xlabel(r'$T_{CA}$ [yr]')
				plt.pause(0.01)
				if k == 3: continue
				if k == 6: 
					z = 1
					c = ps.blue
					label = r'$\Delta G < 6$ mag'
				if k == 20: 
					z = 0
					c = ps.grey 
					label = 'All events'
				if j[1] == 2070:
					if i == 'All' or i == 'MS':
						fig = plt.figure('Results_Hist_'+ i, \
							figsize = [12.8,2.4])
						plt.subplots_adjust(\
							left = plt.rcParams['figure.subplot.left']/2,
							right= plt.rcParams['figure.subplot.right']/2+0.5,
							bottom = plt.rcParams['figure.subplot.bottom']*2)
					else: 
						fig = plt.figure('Results_Hist_' + i,\
						figsize = [6.4,2.4])
						plt.subplots_adjust(\
							bottom = plt.rcParams['figure.subplot.bottom']*2)
					plt.hist(result['TCA'][which],color = c, bins = 60, \
						range= j, rwidth = 0.9, zorder = z, label = label)
				else:
					print(k,z,c)
					fig = plt.figure('Results_Hist_2030_' + i, \
						figsize = [6.4,2.4])
					plt.subplots_adjust(\
							bottom = plt.rcParams['figure.subplot.bottom']*2)
					plt.hist(result['TCA'][which], color = c, range= j,\
						bins = (j[1]-j[0])*2,rwidth = 0.9,\
						weights = np.ones(len(result[which]))*2, zorder = z, label = label)
				leg = plt.legend()

				plt.xlim(j)
				plt.ylabel(r'Events per year')
				plt.xlabel(r'$T_{CA}$ [yr]')

	# plot magnification 
	for i in ['All','WD','RG','MS','BD']:	
		if i == 'All':
			which = np.ones(len(result),bool)
		else:
			which = result['star_type']==i

		fig = plt.figure('Magnification_'+i)
		which = which & (result['magnification']> 0.001)
		ax = plt.gca()
		twopar = which & (result['ob_parallax'] == 0)
		fivepar = which & (result['ob_parallax'] != 0)
		plt.semilogy(result[twopar]['TCA'],result[twopar]['magnification'],'.', color = 'grey',zorder = -5,
			label = '2 parameter solution')
		plt.errorbar(result[fivepar]['TCA'],result[fivepar]['magnification'],fmt = 'o',\
			xerr = result[fivepar]['TCA_error'],
			yerr = [np.minimum(-result[fivepar]['magnification_error_m'],\
			result[fivepar]['magnification']-0.0001),\
			result[fivepar]['magnification_error_p']],\
			color = 'blue', ms = 2,lw = 0.5, markeredgewidth = 0.2,\
			markeredgecolor = 'k', label = '5 parameter solution')
		plt.ylim([0.001,5])
		plt.yticks([0.001,0.005,0.01,0.05,0.1,0.5,1,5],\
			['0.001', '0.005','0.01','0.05','0.1','0.5','1','5'])
		plt.xlim(2010,2070)
		leg = plt.legend()
		for lh in leg.legendHandles: 
			try: 
				lh._legmarker.set_alpha(1)
				lh._legmarker.set_markersize(5)
			except: pass
		plt.ylabel(r'$\Delta m$ [mag]')
		plt.xlabel(r'$T_{CA}$ [yr]')
		plt.pause(0.01)



def main(n = setup.n_core):
	tt1 = time.time()
	if not os.path.isdir(Folder + 'Results'):
		os.mkdir(Folder + 'Results')
	
	broadcast()
	# load and filter Raw Candidates (see raw_data.py, good_BGS.py, good_HPMS.py & Good_Pairs.py)
	good_raw_cands = raw_data.main()

	good_raw_cands = determine_mass(good_raw_cands)

	#replace masked values
	good_raw_cands = raw_data.fill_value(good_raw_cands) 

	good_raw_cands = determine_einstein_radii(good_raw_cands)
	
	#determine aproximated parameters
	good_raw_cands = approx(good_raw_cands)

	# do precise calculations 
	table_out = calc(good_raw_cands, n)

	#filter events
	table_out = filter_events(table_out)

	tt2 = time.time()
	print('Duration: %im:%is' % ((tt2-tt1) // 60,(tt1-tt2) % 60))
	print('Number of events:', len(table_out))
	print('-------------------------')

	if setup.make_plots: # make plots
		print('Update Plots')
		plot_results(table_out)

		#plot results into the CMD created in mass.approx_Mass()
		plt.figure('CMD')
		#select only events with full photometry
		gg = np.where(table_out['phot_bp_mean_mag'] != 0) 
	
		plt.plot(table_out['phot_g_mean_mag'] - table_out['phot_rp_mean_mag'],\
			table_out['phot_bp_mean_mag'] + 5*np.log10(table_out['parallax']/100),\
			'o', color = ps.blue, ms = 0.5, label = 'predicted events', \
			zorder = 2)

		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		plt.pause(0.01)

		#save figures into Images folder 
		if os.path.isdir(Folder+'Images') == False:
			os.mkdir(Folder+'Images')
		for im_name in plt.get_figlabels(): # loop over every image
			f = plt.figure(im_name)
			plt.pause(0.01) 
			prefix = setup.prefix
			print('Save Figure: '+Folder+'Images/'+im_name+prefix+'.png')
			f.savefig(Folder+'Images/'+im_name+prefix+'.png', transparent =False)
			#print('Save Figure: '+Folder+'Images/'+im_name+prefix+'.pdf')
			#f.savefig(Folder+'Images/'+im_name+prefix+'.pdf', transparent =False)

	print('DONE')



	return(table_out)

