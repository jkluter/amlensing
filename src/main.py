import numpy as np # version >= 1.19.5
import matplotlib.pyplot as plt # version >= 3.3.3
from matplotlib.ticker import LogLocator
import os
import sys
import time
from astropy.table import vstack, Table , MaskedColumn # version >= 4.2

import raw_data
import good_BGS as GB
import mass
import microlensing
import calculation 
import find_closest 
import setup
import plot_functions as pf
from plot_functions import plt 
Folder = setup.Folder
print(Folder)

def broadcast():
	raw_data.broadcast()
	
def determine_mass(cand):
	# calculate Mass (see mass.py)
	print('calc assumed mass ')
	cMass = mass.approx_Mass(cand)
	cand['star_type'] = cMass[2]
	cand['mass'] = cMass[0]
	cand['mass_error'] = cMass[1]
	return cand

def determine_einstein_radii(cand):
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
	print('Cand:', len(cand))
	print('-------------------------')
	return cand

def calc(cand, n_core = 4):
	# split on multiple cores
	if n_core > len(cand): n_core = len(cand)
	if n_core > 1: 
		from joblib import Parallel, delayed # version >= 1.0.0
		import pickle	
		print('Split process on %i Cores'%n_core)

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
		print('\nCombine Results')
		table_out = vstack(table_out)

	else:
		#single computing 
		#(see calculation.py, forwars to microlensing.py and find_closest.py)
		table_out = calculation.single(cand)
	print('-------------------------')
	return(table_out)

def filter_events(table_out):	
	print('Filter events')
	prefix = setup.prefix
	form_out = setup.form_out

	if setup.save_table_process: # save step to file
		

		print('save unfiltered results')
		# save all events
		print(Folder + 'Results/amlensing.unfiltered'+ prefix + form_out[0])
		table_out.write(Folder + 'Results/amlensing.unfiltered' + prefix \
			+ form_out[0], format = form_out[1], overwrite = True)
	# filter results 
	shift_01 = table_out['shift_plus']+table_out['shift_plus_error_m'] > 0.1
	shift_L2_01 = table_out['L2_shift_plus']+table_out['shift_plus_error_m']>0.1
	shift_01 = table_out['shift_plus']> 0.1
	shift_L2_01 = table_out['L2_shift_plus']>0.1
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
		for line in BL[BL.colnames[0]]:
			if line[0] == '#':
				if line[:5] == '# NO ':
					psi = 1
			elif psi: BL2.append(int(line))
			else: BL1.append(int(line))
		BL1 = np.isin(table_out['ob_source_id'],BL1) \
			| np.isin(table_out['source_id'],BL1)
		BL2 = np.isin(table_out['ob_source_id'],BL2) \
			| np.isin(table_out['source_id'],BL2)

		aa = np.sum(FF)
		bb = np.sum(FF & (BL1 == False))
		cc = np.sum(FF & (BL2 == False))
		print('Excluded due to large Psi:', aa-bb)
		print('Excluded due to no proper motion:', aa-cc)

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
			BGS = GB.BGS
			pf.plot_psi_result_part_2(BGS,table_out,blacklist,FF)

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

def main(n = setup.n_core):
	plt.ioff()

	tt1 = time.time()
	if not os.path.isdir(Folder + 'Results'):
		os.mkdir(Folder + 'Results')
	
	broadcast()
	# load and filter Raw Candidates 
	# (see raw_data.py, good_BGS.py, good_HPMS.py & Good_Pairs.py)
	good_raw_cands = raw_data.main()

	good_raw_cands = determine_mass(good_raw_cands)

	# replace masked values
	good_raw_cands = raw_data.fill_value(good_raw_cands) 

	good_raw_cands = determine_einstein_radii(good_raw_cands)
	
	# determine aproximated parameters
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
		pf.plot_results(table_out)
		#plot results into the CMD created in mass.approx_Mass()
		pf.plot_CMD(table_out)
		if os.path.isdir(Folder+'Images') == False:
			os.mkdir(Folder+'Images')
		#save figures into Images folder 
		pf.save_figures(Folder)

	print('DONE')



	return table_out 

