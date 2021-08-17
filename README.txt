Description 

Searches for astrometric microlensing events from a list of pairs, each containing a high proper motion stars and a background star.
Uses several filters based on Gaia_source column's to remove spurious
high-proper-motion stars and background sources  
Using parallax, proper motion, and positional displacement between GaiaDR2 und Gaia DR3
to detect co-moving stars. 
Estimates an assumed mass of the HPMS based on its asolute G magnitude. 
Forecast the motion until 2066 (GaiaeDR3 epoch + 50yrs) and determine the minimum angular separation using a nested intervals algorithm


To download the standard files from Klüter et al 2021 use:

python amlensing download_data [-random_sample] [-n] [- comands]

	-random_sample 		download also a random gaia sample. These are only used for plots but not necessary)
	-n			size of the random sample default  = 100,000) 

Data are downloaded to the ./Data directory, following the file names defined in setup.py.



Using 

python amlensing [- comands]

comands: (defaultvalues are stored in setup.py)
	-r --raw_cands_table file:		use the given file as raw-candidates
	-H --HPMS ["py" | file]:   		use the given file in order to determine good HPMS's. 
				 		If set to 'py' is given uses the good_HPMS.py script.
						If file does not exist, no filter will be aplied.
	-B --BGS ["py" |file]: 			use the given file in order to determine good BGS's 
						If set to 'py' is given uses the good_BGS.py script. 
						If file does not exist, no filter will be aplied.
	-n --n_core [int]:			number of cores used for parallel computing
	-m --make_plots [bool]:			creat plots on the fly (None = True)
	-s --save_table [bool]:			saves different steps as table (None = True) 
	-f --filter [bool]:			use filters for pairs (None = True)
	-p --prefix [str]:			prefix for stored tables 
	-b --blacklist [file]: 			list of ob_source_id's to remove from results
						If file does not exist, no sources will be removed.
Default parameters are defined in the setup.py script.


Requirements
Developed & teste with python 3.7.9 & 3.8.7
numpy >=1.19.5
matplotlib >=3.3.3
joblib >=1.0.0
astropy >=4.2
astroquery >= 0.4.1
pyvo >= 1.1

