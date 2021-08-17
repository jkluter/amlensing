'''
load a list of all High proper motion stars
must contains the Gaia source ID = 'hpms_source_id'
contains a boolean if the HPMS can be used (hpms_good), 
if not all HPMS in the list are used 

load a list of all background stars, 
must contains the Gaia source ID = 'source_id'
contains a boolean if the HPMS can be used (bgs_good), 
if not all bgs in the list are used
'''


import numpy as np
from astropy.table import Table, MaskedColumn, Column

#sub packages
from setup import Folder, hpms_file,bgs_file,raw_cands_table,\
	save_table_process, make_plots, mask_values,form_out, prefix, do_filter
import  good_Pairs as GP 
import good_HPMS as GH
import good_BGS as GB
import plot_functions as pf
import time
import os
tt = [0]



def broadcast():
	global hpms_file,bgs_file,raw_cands_table,\
		save_table_process, make_plots, mask_values,form_out, prefix, do_filter
	from setup import Folder, hpms_file,bgs_file,raw_cands_table,\
		save_table_process, make_plots, mask_values,form_out, prefix, do_filter
	GB.broadcast()
	GH.broadcast()
	GP.broadcast()
def fill_value(tab):
	# replace None values with default value
	# see setup for values
	for colname in mask_values.keys():
		eq1e20 = tab[colname]== 1e20

		tab[colname].fill_value = mask_values[colname]
		tab[colname][eq1e20] = mask_values[colname]
	return tab.filled()

def load_raw_pairs():
	global raw_cands
	print('Load raw pairs')
	if '.vot' in raw_cands_table:
		raw_cands = Table.read(Folder + 'Data/' + raw_cands_table, \
			format = 'votable')
	if '.fits'in raw_cands_table:	
		raw_cands = Table.read(Folder + 'Data/' + raw_cands_table, \
			format = 'fits')
		# convert metadata to astropy description
		if 'TCOMM1' in raw_cands.keys():
			for j,i in enumerate(raw_cands.keys().copy()):
				raw_cands[i].description = raw_cands.meta.pop('TCOMM%i'%(j+1))

	# mask None values	
	if 'roi' in raw_cands.colnames: 
		raw_cands.remove_column("roi")
	raw_cands = Table(raw_cands, masked = True, copy = False)
	for i in raw_cands.colnames:
		if any(np.isnan(raw_cands[i].astype(float))): 
			raw_cands[i].mask = np.isnan(raw_cands[i].astype(float))
		else: 
			raw_cands[i].mask = raw_cands[i].astype(float) == 1e20


def check_hpms():
	# load whitelist of HPMS source_id's
	global good_raw_cands_HPMS, hpms_file
	if 'raw_cands' not in globals():
		load_raw_pairs()
	print('-------------------------')
	if hpms_file == 'good_HPMS.py':
		# Find good HPMS using the good_HPMS.py script and the full gaia table 
		print('Find good HPMS') 
		good_HPMS_source_ID = GH.main()
		#check if source id in whitlist
		good_raw_cands_HPMS = raw_cands[np.isin(
			raw_cands['source_id'],	good_HPMS_source_ID)]
	elif os.path.isfile(hpms_file) or os.path.isfile(Folder + 'Data/'+hpms_file):
		# load pre determined whitelist
		print('Load good_hpms')
		if os.path.isfile(hpms_file) == False:
			hpms_file = Folder + 'Data/'+hpms_file
		if hpms_file.split(".")[-1] == "vot" :
			hpms = Table.read(hpms_file, format = 'votable')
		else:
			hpms = Table.read(hpms_file, format = hpms_file.split(".")[-1])
		# does the file contains a boolean to indicate good values
		if 'hpms_good' in hpms.colnames:
			good_hpms = hpms[hpms['hpms_good']]
		elif 'fidelity_v1' in hpms.colnames:
			good_sus = hpms['fidelity_v1']> 0.7
			good_phot = (GH.HPMS_check_phot(cat = hpms))
			good=good_sus &good_phot
		else: 
			good = np.ones(len(hpms), dtype = bool)
		#check if source id in whitlist
		good_hpms = hpms[good]



		good_raw_cands_HPMS = raw_cands[np.isin(
			raw_cands['source_id'], good_hpms['source_id'])]
		print('HPMS good:',len(good_hpms),"/",len(hpms))
		print('HPMS index:',len(good_raw_cands_HPMS),"/",len(raw_cands))
		if save_table_process:
			print('save good HPMS')
			print(Folder + 'Results/HPMS.good' + prefix+ form_out[0])
			print(Folder + 'Results/HPMS.bad' + prefix+ form_out[0])
			bad_hpms = hpms[good == False]
			if 'good_sus' in locals():
				bad_hpms["out_sus"] = good_sus[good == False] ==False
				bad_hpms["out_phot"] = good_phot[good == False] ==False

			good_hpms.write(Folder + 'Results/HPMS.good' + prefix\
					+ form_out[0], format = form_out[1], overwrite = True)
		
			bad_hpms.write(Folder + 'Results/HPMS.bad' + prefix\
					+ form_out[0], format = form_out[1], overwrite = True)
	else: 
		#do not performe HPMS filtering 
		print('No hpms file found')
		good_raw_cands_HPMS = raw_cands

def check_bgs():
	# load whitelist of HPMS source_id's
	global	good_raw_cands_BGS
	if 'good_raw_cands_HPMS' not in globals():
		check_hpms()
	print('-------------------------')
	if bgs_file == 'good_BGS.py':
		# Find good BGS using the good_HPMS.py script and the full gaia table 
		print('Find good BGS')
		good_BGS_source_ID, _ , pmDR2 = GB.main()
		#check if source id in whitlist
		good_raw_cands_BGS = good_raw_cands_HPMS[np.isin(
			good_raw_cands_HPMS['ob_source_id'], good_BGS_source_ID)]

		# include DR2- DR3 propermotion
		DR2_pm = np.array([pmDR2[i] if i in pmDR2.keys() else [0,0] \
			for i in good_raw_cands_BGS['ob_source_id']])
		good_raw_cands_BGS['ob_displacement_ra_doubled'] = \
			MaskedColumn(DR2_pm[:,0],unit = 'mas', description = 'Doubled ' \
			+ 'Displacemnet in RA between DR2 and DR3. (cos(DEC) applied)')
		good_raw_cands_BGS['ob_displacement_dec_doubled'] = \
			MaskedColumn(DR2_pm[:,1], unit = 'mas', description = 'Doubled ' \
			+ 'Displacemnet in DEC between DR2 and DR3.')

	elif os.path.isfile(bgs_file): 
		# load pre determined whitelist
		print('Load good_bgs')
		bgs = Table.read(bgs_file, format = 'votable')
		# does the file contains a boolean to indicate good values
		if 'bgs_good' in bgs.colnames:
			good_bgs = bgs[bgs['bgs_good']]
		else: 
			good_bgs = bgs
		#check if source id in whitlist
		good_raw_cands_BGS = good_raw_cands_HPMS[np.isin(
			good_raw_cands_HPMS['obj_source_id'], good_bgs['source_id'])]
	else: 
		#do not performe BGS filtering 
		good_raw_cands_BGS = good_raw_cands_HPMS
		good_BGS_source_ID, _ , pmDR2 = GB.main()
		DR2_pm = np.array([pmDR2[i] if i in pmDR2.keys() else [0,0] \
			for i in good_raw_cands_BGS['ob_source_id']])

		good_raw_cands_BGS['ob_displacement_ra_doubled'] = \
			MaskedColumn(DR2_pm[:,0],unit = 'mas', description = 'Doubled ' \
			+ 'Displacemnet in RA between DR2 and DR3. (cos(DEC) applied)')

		good_raw_cands_BGS['ob_displacement_dec_doubled'] = \
			MaskedColumn(DR2_pm[:,1], unit = 'mas', description = 'Doubled ' \
			+ 'Displacemnet in DEC between DR2 and DR3.')
		print('No bgs file found')
def check_pairs():
	# filter pairs based on the comparison between lens and source data
	# i.e exclude binary stars 
	if 'good_raw_cands_BGS' not in globals():
		check_bgs()
	print('-------------------------')
	print('Find good Pairs')
	global good_raw_cands, do_filter
	# apply filters from good_Pairs.py

	good = np.ones(len(good_raw_cands_BGS), dtype = bool)
	if do_filter == 1: do_filter = 1111111
	i =0
	goodfilter = []
	for ff in GP.All_filter: 
		i+=1
		if (do_filter//10**i)%10:
			good_ff = ff(good_raw_cands_BGS)
			print(ff.__name__, np.sum(good_ff), np.sum(good_ff==False), 
				len(good_ff))
			goodfilter.append([ff.__name__, good_ff])
			good= good & good_ff
		else:
			print("skip", ff.__name__)
	good_raw_cands= good_raw_cands_BGS[good]
	# save  good raw_cands table to disk
	if save_table_process:
		print('save good rawcands')
		print(Folder + 'Results/rawcands.good' + prefix+ form_out[0])
		good_raw_cands.write(Folder + 'Results/rawcands.good' + prefix+ form_out[0], \
			format = form_out[1], overwrite = True)
		print(Folder + 'Results/rawcands.bad' + prefix+ form_out[0])
		bad_raw_cands =good_raw_cands_BGS[good==False]
		for i in goodfilter:
			bad_raw_cands[i[0][7:]] = Column(i[1][good==False]==False) 

		bad_raw_cands.write(Folder + 'Results/rawcands.bad' + prefix+ form_out[0], \
			format = form_out[1], overwrite = True)

def main(redo = 0):
	# initialate loading the raw_cands table and going throug the different
	# filter types
	# returns good raw candidates

	broadcast()

	tt = []
	tt.append(time.time())
	if 'raw_cands' not in globals() or redo < 2:
		load_raw_pairs()
	tt.append(time.time())
	if 'good_raw_cands_HPMS' not in globals()  or redo < 3:
		check_hpms()

	tt.append(time.time())
	if 'good_raw_cands_BGS' not in globals()  or redo < 4:
		check_bgs()

	tt.append(time.time())	
	if 'good_raw_cands' not in globals()  or redo < 5 :
		check_pairs()

	tt.append(time.time())
	tt = np.array(tt)
	cpt = tt[1:] - tt[:-1]
	#print('Load RawData:', *(cpt))

	if make_plots:
		bgs_good = GB.BGS[np.isin(
			GB.BGS['source_id'], good_raw_cands["ob_source_id"])]
		pf.plot_psi_part2(bgs_good)

		pf.plot_pos_err(data=bgs_good)

		pf.plot_sim_px(data=good_raw_cands)

		pf.plot_HPMS_part2(good_raw_cands)

	print('-------------------------')
	print('Good raw cands', len(good_raw_cands) )
	print('-------------------------')

	return good_raw_cands

if __name__ == '__main__': 
	print(__file__)
