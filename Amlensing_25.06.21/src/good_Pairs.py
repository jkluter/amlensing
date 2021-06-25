import numpy as np
import matplotlib.pyplot as plt
import time

#sub packages
from setup import make_plots, Pair_limit
import plt_setup as ps

''' 
list of Filter functions for the combined data of lens and source 
Used in Raw Data
might be improoved 
'''
plt.ion()
ff = []

import importlib
def re():
	import good_Pairs as GP
	importlib.reload(GP)
def broadcast():
	global make_plots, Pair_limit
	from setup import make_plots, Pair_limit
	


def Filter_sim_pm(rp, limit = Pair_limit['pm_sim_1'],\
	limit_sim = Pair_limit['pm_sim_2']):
	# non similar proper motion
	F_pm_1= (rp['pmra']-rp['ob_pmra'])**2 + (rp['pmdec']-rp['ob_pmdec'])**2 \
			> limit_sim**2 * (rp['pmdec']**2 + rp['pmra']**2)
	F_pm_2 = rp['ob_pmra']**2 + rp['ob_pmdec']**2 \
		< limit**2 * (rp['pmdec']**2 + rp['pmra']**2)

	print(np.sum(F_pm_1 & (F_pm_2 == False)))
	print(np.sum(F_pm_2 & (F_pm_1 == False)))

	# Filter is not used for bgs without 5-parameter solution
	F_pm_3 = rp['ob_pmra'].mask
	if make_plots:
		fig = plt.figure("delta_pm", figsize = [6.4,3])
		fig.clf()
		fig.subplots_adjust(right = 0.87, bottom = 0.17)
		ax1 = plt.subplot(121)
		FF = [(F_pm_1 & F_pm_2)]
		
		ps.plot(rp['pmra'] - rp['ob_pmra'],\
			rp['pmdec'] - rp['ob_pmdec'], color = ps.red,  ms = 0.5,\
			label = 'excluded Pairs')
		ps.plot(rp['pmra'][FF] - rp['ob_pmra'][FF], \
			rp['pmdec'][FF] - rp['ob_pmdec'][FF], color = ps.blue, ms = 0.5,\
			label = 'Pairs')
		plt.pause(0.01)
		plt.xlabel(r"$\mu_{\alpha\star} - Sou\_\mu_{\alpha\star}$ [mas/yr]")
		plt.ylabel(r"$\mu_{\delta} - Sou\_\mu_{\delta}$ [mas/yr]")
		plt.xlim([-500,500])
		plt.ylim([-500,500])
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)	
	F_pm_1.fill_value = 0
	F_pm_1 = F_pm_1.filled()	
	F_pm_2.fill_value = 0
	F_pm_2 = F_pm_2.filled()
	return (F_pm_1 & F_pm_2) | F_pm_3

def Filter_sim_pm_DR2(rp, limit = Pair_limit['pm_sim_1'],\
	limit_sim = Pair_limit['pm_sim_2']):

	if "ob_displacement_dec_doubled" not in rp.columns:
		return np.ones(len(rp), dtype = bool)
	# non similar proper motion

	F_pm_1= (rp['pmra']-rp['ob_displacement_ra_doubled'])**2 \
			+ (rp['pmdec']-rp['ob_displacement_dec_doubled'])**2 \
			> limit_sim**2 * (rp['pmdec']**2 + rp['pmra']**2)
	F_pm_2 = rp['ob_displacement_ra_doubled']**2 \
		+ rp['ob_displacement_dec_doubled']**2 \
		< limit**2 * (rp['pmdec']**2 + rp['pmra']**2)
	# Filter is not used for bgs with an 5-parameter solution 
	F_pm_3 = rp['ob_pmra'].mask == False
	if make_plots:
		fig = plt.figure("delta_pm", figsize = [6.4,3])
		fig.subplots_adjust(right = 0.87, bottom = 0.17)
		ax2 = plt.subplot(122)
		FF = [(F_pm_1 & F_pm_2) & (F_pm_3 == False)]
		ps.plot(rp[F_pm_3 == False]['pmra'] \
			- rp[F_pm_3 == False]['ob_displacement_ra_doubled'],\
			rp[F_pm_3 == False]['pmdec'] \
			- rp[F_pm_3 == False]['ob_displacement_dec_doubled'],\
			color = ps.red,  ms = 0.5)
		ps.plot(rp['pmra'][FF] - rp['ob_displacement_ra_doubled'][FF], \
			rp['pmdec'][FF] - rp['ob_displacement_dec_doubled'][FF], \
			color = ps.blue, ms = 0.5)
		plt.xlabel(\
			r"$\mu_{\alpha^{\star}} - \Delta \phi_{\rm \times2,\,\alpha^{\star}}/1{\rm yr}$ [mas/yr]")
		plt.ylabel(r"$\mu_{\delta} - \Delta \phi_{\rm \times2,\,\delta}/1{\rm yr}$ [mas/yr]")
		plt.xlim([-500,500])
		plt.ylim([-500,500])
		ax2.yaxis.tick_right()

		ax2.yaxis.set_label_position("right")

		plt.pause(0.01)
	return (F_pm_1 & F_pm_2) | F_pm_3

def Filter_pm_tot_DR2(rp, limit = Pair_limit['pm_tot']):
	# Filter on the absolut dr2 propermotion
	if "ob_displacement_dec_doubled" not in rp.columns:
		return np.ones(len(rp), dtype = bool)
	F_pmdr2_1 = (rp['ob_displacement_ra_doubled']**2 \
		+ rp['ob_displacement_dec_doubled']**2) < limit*limit
	# Filter is not used for bgs with an 5-parameter solution 
	F_pmdr2_2 = rp['ob_pmra'].mask == False

	return F_pmdr2_1 | F_pmdr2_2


def Filter_sim_px(rp, limit = Pair_limit['px']):
	# non similar parallax
	F_px_1 = rp['ob_parallax']  < limit * rp['parallax']
	F_px_2 = rp['ob_pmra'].mask
	F_px_1.fill_value = 0
	F_px_1=F_px_1.filled()

	if make_plots:


		fig = plt.figure("sim_px")
		ps.semilogx(rp['parallax'], \
			rp['ob_parallax'],  color = ps.red, label='excluded Pairs')
		ps.semilogx(rp[F_px_1]['parallax'], \
			rp[F_px_1]['ob_parallax'], color = ps.blue, label='Pairs')
		plt.xlabel(r'$\varpi$ [mas]')
		plt.ylabel(r'$Sou\_\varpi$ [mas]')
		plt.pause(0.01)
	

	return F_px_1 | F_px_2

All_filter = [Filter_sim_pm, Filter_sim_px, Filter_sim_pm_DR2]

