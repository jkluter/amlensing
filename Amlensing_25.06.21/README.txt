Description 

Uses Rawcands table to determine microlensing events.
Uses several filters based on Gaia_source column's to remove spurious
High Propermotion Stars and BGS
Using parallax, propermotion, and positional shift between GaiaDR2 und Gaia DR3
to detect comoving stars. 
Determine an asumed mass of the HPMS based on its asolute G magnitude. 
Forecast the motion until 2066 (GaiaeDR3 epoch + 50yrs) and determine the minimum angular separation using a nested intervals algorithm


Using 

python Amlensing [- comands]

comands: (defaultvalues are stored in setup.py)
	-r --raw_cands_table file:	use the given file as raw-candidates

	-H --HPMS ["py" | file]:   		use the given file in order to 
									determine good HPMS's. If set to 'py' is 
									given uses the good_HPMS.py script, If 
									file does not exist, no filter will be 
									aplied.
	-B --BGS ["py" |file]: 			use the given file in order to 
									determine good BGS's 
									If set to 'py' is given uses the 
									good_BGS.py script. If file does not 
									exist, no filter will be aplied.
	-n --n_core [int]:				number of cores used for parallel computing
	-m --make_plots [None,bool]:	creat plots on the fly (None = True)
	-s --save_table [None,bool]:	stors different steps as table 
									(None = True)
	-f --filter [None,bool]:		use filters for pairs (None = True)
	-p --prefix [str]:				prefix for stored tables 
	-b --blacklist [file]: 			list of ob_source_id's to remove from 
									results


Requirements
Developed & teste with python 3.7.9 & 3.8.7
numpy >=1.19.5
matplotlib >=3.3.3
joblib >=1.0.0
astropy >=4.2

