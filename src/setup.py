# Amlensing folder
Folder = '/'.join(__file__.split('/')[:-2]) + '/'

# Ignore Warnings
import warnings
warnings.simplefilter("ignore")



#------------------------------------------------------------------------------
#Default values and files 
n_core = 4
prefix = ''
do_filter = 1
# Files has to be located in Data or given as relative path from current working directory.
# Raw Pairs File 
raw_cands_table = 'rawcands.fits' 

# Method to determine good High Proper Motion Stars and BackGround Stars
# If empty string: all source are belived to be good
# If 'good_HPMS.py' & 'good_BGS.py' use this methods to determine good HPMS & good BGS
# If an file is given: read the corresponding file as astropy Table 
# Has to contain Gaia eDR3 source ID's in the 'source_id' colum
# good sources are either flaged by HPMS_good and BGS_good, respectively or all sources are asumed to be good
hpms_file = 'good_HPMS.py'
bgs_file = 'good_BGS.py'

# Files for determine Good HPMS and BGS
HPMS_eDR3_file = 'gaia_edr3_hpms.fits' #all potential HPMS
GCNS_cat_file = 'GCNS_cat.fits' #Gaia catalogue of nearby stars
GCNS_reject_file= 'GCNS_reject.fits' #Gaia catalogue of nearby stars, rejected sources
HPMS_spur_file = 'HPMS_spur.fits' #fidility values from J. Rybizki et al. 2021 (optional)

BGS_eDR3_file = 'gaia_edr3_bgs.fits' #all potential BGS
#gaiaedr3.dr2_neighbourhood for all BGS
DR2_BGS_file = 'gaia_bgs_dr2_neighbours.fits' 
#gaiaedr3 random sample, used for plots only 
random_sample_file = 'random_sample.fits'
dr2_random_file = 'dr2_random.fits'


#Blacklist
Blacklist_file = 'Blacklist.csv'

# create plots on the fly
make_plots = False

# save intermediate tables 
save_table_process = False

# output format 
form_out = ['.fits','fits'] # save data as fits table
#form_out = ['.vot','votable'] # save data as vo table

#------------------------------------------------------------------------------

# use Percentile to determine Errors 
error_percentile = True

# Gaia paralax zero point
zeropoint = 0	 # Used in Filters, 

# Default values for NONE in the raw pairs
mask_values = {'phot_rp_mean_mag': 0,
		'phot_bp_mean_mag': 0,
		'ob_pmra': 0,
		'ob_pmdec': 0,
		'ob_pmra_error': 5,
		'ob_pmdec_error': 5,
		'ob_parallax': zeropoint,
		'ob_parallax_error': 3.,
		'ob_phot_rp_mean_mag': 0,
		'ob_phot_bp_mean_mag': 0}

# Gaia reference epoch
Gaia_epoch = 2016.

#------------------------------------------------------------------------------
# limits for criteria on High Proper Motion Stars 
# used in good_HPMS.py
HPMS_limit = {
	'ruwe': 2.0,
	'px':5,
	'mag': 21,
	'n_obs_sig_g_flux': 3e5,
	'n_obs_sig_g_flux_power': 1.5}
# limits for criteria on BackGround Stars 
# used in good_BGS.py
BGS_limit = {
	'ruwe':2.0,
	'px': -3,
	'pm': 80,
	'mag': 21.5,
	'pos_err':10}
# limits for DR2 Match
# used in good_BGS.py
DR2_limit = {
	'dist': 400, 
	'pm_bad': 1000,
	'mag': 0.3}

# limits for criteria on HPMS-BGS Pairs 
# used in good_Pairs.py
Pair_limit = {
	'px': 0.9,
	'pm_sim_1': 0.8,
	'pm_sim_2': 0.7,
	'pm_tot': 80}
#------------------------------------------------------------------------------




