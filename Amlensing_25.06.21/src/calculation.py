import os
import time
import numpy as np

#sub packages
from find_closest import estimate_Closest_parallax,\
	estimate_errors_parallax
from astropy.table import MaskedColumn
from microlensing import calculate_Effect
from setup import error_percentile as ep

def parallel(string):
	import pickle
	# load input table
	part = int(string.split('p')[-2][:-1])
	f = open(string, 'rb')
	raw_cands = pickle.load(f)
	f.close()
	os.remove(string) # delete temporary file

	print('Find closet approches')
	result_closest = []
	cc = np.zeros(6)
	tt=time.time()
	# loop over row to finde closest approach
	for i, row in enumerate(raw_cands):
		if (i+1)%100 == 0 : # print progress
			print('Part%i - %i/%i: %im:%is '%(part, i+1,len(raw_cands),\
				(time.time()-tt) // 60,(time.time()-tt) % 60))
		data,cpt = estimate_Closest_parallax(row)
		cc += cpt
		result_closest.append(data)
	result_closest = np.array(result_closest)

	#print("parallel:", *cc)
	data = Process_data(raw_cands,result_closest)
	# store output table temporary
	f = open(string, 'wb')
	pickle.dump(data, f)
	f.close()

def single(raw_cands):
	#raw_cands=raw_cands[:500]
	result_closest = []
	print('Find closet approches')
	tt=time.time()
	cc = np.zeros(6)
	for i, row in enumerate(raw_cands):
		if (i+1)%100 == 0 : # print progress
			print('%i/%i: %is '%(i+1,len(raw_cands),time.time()-tt))
		data,cpt = estimate_Closest_parallax(row)
		cc += cpt
		result_closest.append(data)
	result_closest = np.array(result_closest)
	
	output = Process_data(raw_cands,result_closest)
	return output

def Process_data(raw_cands, result_closest):
	print('calculate_Effect')
	#add closest approach at L2 to table	
	raw_cands['TCA'] = MaskedColumn(result_closest[:,0,0],dtype = 'float64',\
		unit = 'year', description = 'Epoch of the closest approch')
	error_earth = estimate_errors_parallax(raw_cands)
	raw_cands['TCA_error'] = error_earth[0]
	raw_cands['dist'] = MaskedColumn(result_closest[:,1,0],dtype = 'float64',\
		unit = 'mas', description = 'Closest distance')
	raw_cands['dist_error'] = error_earth[1]

	# calculate effect
	c_effect = calculate_Effect(raw_cands, error_percentile = ep)
	# add to table
	if len(c_effect) == 20:
		# add to table
		raw_cands['u'] = c_effect[0]
		raw_cands['u_error'] = c_effect[1]
		raw_cands['u_error_p'] = c_effect[2]
		raw_cands['u_error_m'] = c_effect[3]
		raw_cands['shift'] = c_effect[4]
		raw_cands['shift_error'] = c_effect[5]
		raw_cands['shift_error_p'] = c_effect[6]
		raw_cands['shift_error_m'] = c_effect[7]
		raw_cands['shift_plus'] = c_effect[8]
		raw_cands['shift_plus_error'] = c_effect[9]
		raw_cands['shift_plus_error_p'] = c_effect[10]
		raw_cands['shift_plus_error_m'] = c_effect[11]
		raw_cands['shift_lum'] = c_effect[12]
		raw_cands['shift_lum_error'] = c_effect[13]
		raw_cands['shift_lum_error_p'] = c_effect[14]
		raw_cands['shift_lum_error_m'] = c_effect[15]
		raw_cands['magnification'] = c_effect[16]
		raw_cands['magnification_error'] = c_effect[17]
		raw_cands['magnification_error_p'] = c_effect[18]
		raw_cands['magnification_error_m'] = c_effect[19]
	else:
		raw_cands['u'] = c_effect[0]
		raw_cands['u_error'] = c_effect[1]
		raw_cands['shift'] = c_effect[2]
		raw_cands['shift_error'] = c_effect[3]
		raw_cands['shift_plus'] = c_effect[4]
		raw_cands['shift_plus_error'] = c_effect[5]
		raw_cands['shift_lum'] = c_effect[6]
		raw_cands['shift_lum_error'] = c_effect[7]
		raw_cands['magnification'] = c_effect[8]
		raw_cands['magnification_error'] = c_effect[9]

	# add closest approach at L2 to table
	raw_cands['L2_TCA'] = MaskedColumn(result_closest[:,2,0], \
		dtype = 'float64', unit = 'year', description = \
		'Epoch of the closest approch for Lagrange Point L2')
	error_earth = estimate_errors_parallax(raw_cands,gaia = True, )
	raw_cands['L2_TCA_error'] = error_earth[0]
	raw_cands['L2_dist'] = MaskedColumn(result_closest[:,3,0], \
		dtype = 'float64', unit = 'mas', \
		description = 'Closest distance for Lagrange Point L2')
	raw_cands['L2_dist_error'] = error_earth[1]

	# calculate effect at L2
	c_effect = calculate_Effect(raw_cands, gaia= True, error_percentile = ep)
	if len(c_effect) == 20:
		# add to table
		raw_cands['L2_u'] = c_effect[0]
		raw_cands['L2_u_error'] = c_effect[1]
		raw_cands['L2_u_error_p'] = c_effect[2]
		raw_cands['L2_u_error_m'] = c_effect[3]
		raw_cands['L2_shift'] = c_effect[4]
		raw_cands['L2_shift_error'] = c_effect[5]
		raw_cands['L2_shift_error_p'] = c_effect[6]
		raw_cands['L2_shift_error_m'] = c_effect[7]
		raw_cands['L2_shift_plus'] = c_effect[8]
		raw_cands['L2_shift_plus_error'] = c_effect[9]
		raw_cands['L2_shift_plus_error_p'] = c_effect[10]
		raw_cands['L2_shift_plus_error_m'] = c_effect[11]
		raw_cands['L2_shift_lum'] = c_effect[12]
		raw_cands['L2_shift_lum_error'] = c_effect[13]
		raw_cands['L2_shift_lum_error_p'] = c_effect[14]
		raw_cands['L2_shift_lum_error_m'] = c_effect[15]
		raw_cands['L2_magnification'] = c_effect[16]
		raw_cands['L2_magnification_error'] = c_effect[17]
		raw_cands['L2_magnification_error_p'] = c_effect[18]
		raw_cands['L2_magnification_error_m'] = c_effect[19]
	else:
		# add to table
		raw_cands['L2_u'] = c_effect[0]
		raw_cands['L2_u_error'] = c_effect[1]
		raw_cands['L2_shift'] = c_effect[2]
		raw_cands['L2_shift_error'] = c_effect[3]
		raw_cands['L2_shift_plus'] = c_effect[4]
		raw_cands['L2_shift_plus_error'] = c_effect[5]
		raw_cands['L2_shift_lum'] = c_effect[6]
		raw_cands['L2_shift_lum_error'] = c_effect[7]
		raw_cands['L2_magnification'] = c_effect[8]
		raw_cands['L2_magnification_error'] = c_effect[9]
	return(raw_cands)




