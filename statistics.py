from  astropy.table import Table, unique, vstack,hstack, join
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from scipy import stats

sus = 0

#sys.path.append(os.getcwd())
sys.path.append(__file__[:-13]+'src')
print(__file__[:-13]+'src')
if len(sys.argv)>1:
	if sys.argv[1] == '-s':
		sus = 1

from setup import Folder
print(Folder)

if sus:
	dr3_file = Folder+ "Results/amlensing.sus.fits"
	file=Folder + 'Results/amlensing.statistics.sus.txt'
else:
	dr3_file = Folder + "Results/amlensing.fits"
	file= Folder+'Results/amlensing.statistics.txt'

results = Table.read(dr3_file, format = "fits")


# cathegory 
all_results = np.ones(len(results), dtype = bool)
hpms = results['pmra']**2+results['pmdec']**2>150**2
g18 = results["ob_phot_g_mean_mag"]<18
WD = abs(results["mass"]- 0.65)< 1e-6
lowmass = results["mass"] < 0.65 - 1e-6
solarMass = results["mass"] > 0.65 + 1e-6
cathegories = {'all':all_results,
	'pm > 150 mas/yr':hpms,
	'G_source < 18 mag':g18,
	'White Dwarfs':WD,
	'Mass < 0.65Msun':lowmass,
	'Mass > 0.65Msun':solarMass}
deltaG6 = results["ob_phot_g_mean_mag"]-  results["phot_g_mean_mag"] < 6
deltaG3 = results["ob_phot_g_mean_mag"]-  results["phot_g_mean_mag"] < 3
# criterium 
shift1 = results['shift_plus']>1
shift05 = (results['shift_plus']>0.5)& (results['shift_plus'] < 1 )
shift1lum = results['shift_lum']>1 
shift05lum = (results['shift_lum']>0.5) & (results['shift_lum'] < 1 )
shift01lum = results['shift_lum']>0.1  
nextdecade = (results['TCA']<2031) & (results['TCA']>2021)
mag01 = results['magnification']>0.001
mag10 = results['magnification']>0.01
mag100 = results['magnification']>0.1

open(file, 'w').close()
file_output= open(file,'a')

for cat in cathegories.keys():
	current = cathegories[cat]
	if 'White' in cat:  	file_output.write("\n")
	if 'pm' in cat:  	file_output.write("\n")
	file_output.write("----------------------------------------------------\n")
	file_output.write(cat+":\n")
	diffstars = len(np.unique(results[current]['source_id']))
	n_events = [np.sum(current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_nextdecade = [np.sum(nextdecade&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	dg0 = np.sum(current & (results["ob_phot_g_mean_mag"]-results["phot_g_mean_mag"] < 0))
	evp = []
	evp1m = []
	for dg in [all_results,deltaG6,deltaG3]:
		evp_list = [np.sum(abs(results['TCA'][current&dg]-y)<0.25)*2 for y in np.arange(2011,2065.,.25)]
		evp.append(np.median(evp_list))
		evp1m_list = [np.sum(abs(results['TCA'][current&shift1&dg]-y)<0.5) for y in np.arange(2025,2065.,.5)]
		evp1m.append(np.median(evp1m_list))
	n_s1 = [np.sum(shift1&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_s1_nextdecade = [np.sum(nextdecade&shift1&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_s05 = [np.sum(shift05&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_s05_nextdecade = [np.sum(nextdecade&shift05&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	
	n_l1 = [np.sum(shift1lum&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_l1_nextdecade = [np.sum(nextdecade&shift1lum&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_l05 = [np.sum(shift05lum&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_l05_nextdecade = [np.sum(nextdecade&shift05lum&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_l01 = [np.sum(shift01lum&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_l01_nextdecade = [np.sum(nextdecade&shift01lum&current&dg) for dg in [all_results,deltaG6,deltaG3]]


	n_mag = [np.sum(mag01&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_mag_nextdecade = [np.sum(nextdecade&mag01&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_mag10 = [np.sum(mag10&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_mag10_nextdecade = [np.sum(nextdecade&mag10&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_mag100 = [np.sum(mag100&current&dg) for dg in [all_results,deltaG6,deltaG3]]
	n_mag100_nextdecade = [np.sum(nextdecade&mag100&current&dg) for dg in [all_results,deltaG6,deltaG3]]

	file_output.write("events:           %i (%i, %i)\n"%tuple(n_events))
	file_output.write("dG<0:             %i \n"%(dg0))
	file_output.write("different lenses: %i \n"%(diffstars))
	file_output.write("-events per year: %i (%i, %i)\n"%tuple(evp))
	file_output.write("-next decade:     %i (%i, %i)\n"%tuple(n_nextdecade))
	file_output.write("-shift > 1 mas:   %i (%i, %i)\n"%tuple(n_s1))
	file_output.write("--events per year:%i (%i, %i)\n"%tuple(evp1m))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_s1_nextdecade))
	file_output.write("-shift > 0.5 mas: %i (%i, %i)\n"%tuple(n_s05))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_s05_nextdecade))
	
	file_output.write("Luminous include: \n")
	file_output.write("-shift > 0.1 mas: %i (%i, %i)\n"%tuple(n_l01))
	file_output.write("-next decade:     %i (%i, %i)\n"%tuple(n_l05_nextdecade))
	file_output.write("-shift > 1 mas:   %i (%i, %i)\n"%tuple(n_l1))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_l1_nextdecade))
	file_output.write("-shift > 0.5 mas: %i (%i, %i)\n"%tuple(n_l05))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_l05_nextdecade))

	file_output.write("Magnification: \n")
	file_output.write("-mag > 1 mmag:    %i (%i, %i)\n"%tuple(n_mag))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_mag_nextdecade))
	file_output.write("-mag > 10 mmag:   %i (%i, %i)\n"%tuple(n_mag10))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_mag10_nextdecade))
	file_output.write("-mag > 100 mmag:  %i (%i, %i)\n"%tuple(n_mag100))
	file_output.write("--next decade:    %i (%i, %i)\n"%tuple(n_mag100_nextdecade))

	file_output.write("----------------------------------------------------\n")




