import numpy as np
from numpy.polynomial.polynomial import polyval
from astropy.table import MaskedColumn
import time
from setup import make_plots
import plot_functions as pf
'''
determine a approximated mass for the lens based on Gaia data
'''
def approx_Mass(tab):
	print('Determine approximate Mass')

	"""
	estimate the mass of an lense from photometric properties
	classify star in White Dwarfs, Red Giant and Main Sequence
	For MS determine approximatly Mass using an rough M-G_abs relation.
	"""
	tt = []
	tt.append(time.time())
	# Absolut Bp & G magnitude
	no_bprp = \
		(tab['phot_bp_mean_mag'] == 1e20) | (tab['phot_rp_mean_mag'] == 1e20)\
		| np.isnan(tab['phot_bp_mean_mag']) | np.isnan(tab['phot_rp_mean_mag']) 
	B_abs = tab['phot_bp_mean_mag'] + 5 * np.log10(tab['parallax'] / 100)
	G_abs = tab['phot_g_mean_mag'] + 5 * np.log10(tab['parallax'] / 100)
	# G - Rp color
	g_rp = tab['phot_g_mean_mag'] - tab['phot_rp_mean_mag']
	tt.append(time.time())

	# Mass - G_abs function
	popt = [7.8623e-03, -2.9089e-01, 1.1825e+00, -3.0118e-01, 1.8895e+00]
	mass = np.exp(popt[0] * G_abs * G_abs + popt[1] * G_abs + popt[2])
	mass[G_abs > 8.85] = np.exp(popt[3] * G_abs[G_abs > 8.85] + popt[4])
	tt.append(time.time())

	# set BD to 0.07 Msun
	mass = np.maximum(mass, 0.07)
	mass_err = 0.1 * mass
	mass_err[mass < 0.1] = 0.03
	tt.append(time.time())
	
	star_type = np.repeat('MS',len(tab)) 
	star_type[mass == 0.07] = 'BD'
	tt.append(time.time())

	# Define WD
	WD_Terms = [7.4,4.5,4]
	WD = (polyval(g_rp,WD_Terms)  < B_abs) & (no_bprp == False)
	mass[WD] = 0.65
	mass_err[WD] = 0.15
	star_type[WD] = 'WD'
	# Define RG
	RG_Terms = [-20,70,-50.]
	RG = (polyval(g_rp,RG_Terms) > B_abs) & (no_bprp == False)

	mass[RG] = 1
	mass_err[RG] = 0.5
	star_type[RG] = 'RG'

	tt.append(time.time())
	# transform to astropy Columns
	star_type = MaskedColumn(star_type, dtype = 'str', description = \
		'Type of the lensing star: WD = White Dwarf, MS = Main Sequence, '
		+ 'RG = Red Giant, BD = Brown Dwarf')
	mass = MaskedColumn(mass,dtype = 'float64', unit = 'M_sun', description = \
		'Mass of the lensing star')
	mass_error = MaskedColumn(mass_err,dtype = 'float64', unit = 'M_sun',\
		description = 'Error of the Mass of the lensing star')
	tt.append(time.time())
	tt = np.array(tt)
	cpt = tt[1:] - tt[:-1]
	if make_plots:
		pf.plot_CMD(tab, g_rp,B_abs,WD_Terms,RG_Terms)

	return mass, mass_error, star_type, cpt
		