import numpy as np
from astropy.table import MaskedColumn
import time
import plt_setup as ps
from setup import make_plots
import matplotlib.pyplot as plt 

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
	WD = 4. * g_rp * g_rp + 4.5 * g_rp + 7.4 < B_abs
	mass[WD] = 0.65
	mass_err[WD] = 0.15
	star_type[WD] = 'WD'
	# Define RG
	RG = -50. * g_rp * g_rp + 70. * g_rp - 20. > B_abs

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
		plt.figure('CMD')
		gg = np.where(tab['phot_bp_mean_mag'] != 0)
		ps.plot(g_rp[gg], B_abs[gg],\
			'o', color = ps.grey, ms = 0.5, label =  'Candidates', zorder = 0)

		xlim = [-0.5,2.5]
		ylim = [24,-2]
		x = np.linspace(*xlim, 1000)
		plt.plot(x, 4. * x * x + 4.5 * x + 7.4, ps.green, label = 'WD limit',zorder = 4)
		plt.plot(x, -50. * x * x + 70. * x - 20., ps.red, label = 'RG limit', zorder = 3)
		plt.xlim(xlim)
		plt.ylim(ylim)
		plt.xlabel(r'$G_{RP} - G$ [mag]')
		plt.ylabel(r'$G_{BP,abs}$ [mag]')
		plt.pause(0.01)

	return mass, mass_error, star_type, cpt
		