
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LogLocator, AutoMinorLocator, FuncFormatter
from numpy.polynomial.polynomial import polyval
import setup

formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
#ax.xaxis.set_major_formatter(formatter)
#ax.yaxis.set_major_formatter(psformatter)
ion = False

mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.dpi'] = 200
mpl.rcParams['figure.figsize'] = [6.4,4.8]
mpl.rcParams["xtick.direction"] = "in"
mpl.rcParams["ytick.direction"] = "in"
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['errorbar.capsize'] =2

grey = [0.5,0.5,0.5,1]
green =  np.array([1,143,53,255])/255
red = np.array([225,0,35,255])/255
blue = np.array([0,66,148,255])/255
cyan = np.array([0,176,230,255])/255
darkgrey = [0.2,0.2,0.2,1]
yellow = "gold"
limitcolor = 'orange'
colorchange =  np.array([0.6,0.6,0.6,0.05])

def loglog(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.loglog(x,y,fmt, ms = ms, fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.loglog(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)
	return out


def semilogy(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.semilogy(x,y,fmt, ms = ms,fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.semilogy(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)

	return out


def semilogx(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.semilogx(x,y,fmt, ms = ms, fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.semilogx(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)

	return out

def plot(x,y,fmt="o", ms = 2, color = [0., 0.25882353, 0.58039216, 1.],\
	label= None, **kwargs):
	from  matplotlib.colors import cnames, to_rgba_array
	if isinstance(color, str):
		color = to_rgba_array(cnames[color])[0]
	color2 = color * colorchange
	out = plt.plot(x,y,fmt, ms = ms, fillstyle = 'full',\
		color = color,label = label, **kwargs)
	_ = plt.plot(x,y,fmt, ms = ms/2, fillstyle = 'full',\
		color = color2, **kwargs)
	return out





# Plot functions used in  good_BGS.py 

def plot_gof(BGS,out, limit):
	f = lambda x: -4.80459 + 0.000520143 * np.sqrt(4.964e7 * x + 3.57727e7)
	fig = plt.figure("gof")
	fig.clf()
	loglog(BGS['astrometric_n_good_obs_al'],\
		BGS['astrometric_gof_al'],color = red, \
		ms = 0.5 , label = "excluded sources" )
	loglog(BGS['astrometric_n_good_obs_al'][out],\
		BGS['astrometric_gof_al'][out],color = blue, ms = 0.5, \
		label = "Sou_GoF/sqrt(Sou_N) < %.2f"%f(limit))
	good_ruwe = BGS['ruwe'] < limit
	loglog(BGS['astrometric_n_good_obs_al'][good_ruwe],\
		BGS['astrometric_gof_al'][good_ruwe],color = green, \
		ms = 0.5,  label = "RUWE < %.1f"%limit)		
	xlim = np.array(plt.xlim())
	plt.loglog(xlim, np.sqrt(xlim) * f(limit),color = limitcolor,\
		linewidth = 1)
	# plt.xlim(xlim)
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	plt.ylabel('Sou_GoF')
	plt.xlabel('Sou_N')
	if ion: plt.pause(0.1)

def plot_px(BGS,out):
	fig = plt.figure("px")
	fig.clf()
	_,b,_ = plt.hist(BGS['parallax'], color= red, bins = 1000,rwidth = 0.8)
	plt.hist(BGS[out]['parallax'],  color= blue, bins = b, rwidth = 0.5)
	#plt.yscale("log")
	plt.xlim([-5,5])
	plt.ylabel('#')
	plt.xlabel(r'$SOU\_\varpi$')
	if ion: plt.pause(0.1)

def plot_psi(BGS,psi,out,random_sample = None):
	fig = plt.figure("psi")
	fig.clf()

	semilogy(BGS['phot_g_mean_mag'][out], psi[out],\
		color = blue,ms = 0.5,zorder = 1, label = "All BGS")
	semilogy(BGS['phot_g_mean_mag'][out == False], psi[out == False], \
		color = red, ms = 0.5, label = r"$\Psi > 1$ & G < 18 mag", \
		zorder = 0)
	if random_sample is not None:
		rgamma = np.maximum(pow(10, \
			0.2 * (random_sample['phot_g_mean_mag'] - 18)), 1)
		rpsi = random_sample['astrometric_sigma5d_max'] / (1.2 * rgamma)
		semilogy(random_sample['phot_g_mean_mag'], rpsi,\
			color = grey, ms = 0.05,  zorder = 3, \
			label = "random sample")
		c = np.arange(5.5,22,0.1)
		d = np.array([np.percentile(
			rpsi[np.abs(random_sample['phot_g_mean_mag']-i) < 1],90) 
			for i in c])
		plt.plot(c,d, color = limitcolor, label ="90th percentile",zorder = 100)
	
	plt.ylabel(r"$\Psi$")
	plt.xlabel(r"G [mag]")
	plt.ylim([10**-3,10**8])

	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.1)

def plot_psi_result_part_1(random_sample):
	fig = plt.figure("psi_result")
	fig.clf()
	rgamma = np.maximum(
		pow(10, 0.2 * (random_sample['phot_g_mean_mag'] - 18)), 1)
	rpsi = random_sample['astrometric_sigma5d_max'] / (1.2 * rgamma)
	semilogy(random_sample['phot_g_mean_mag'], rpsi,\
		color = grey, ms = 0.05,  zorder = 0, label = "random sample")
	c = np.arange(5.5,22,0.1)
	d = np.array([np.percentile(
		rpsi[np.abs(random_sample['phot_g_mean_mag']-i) < 1],90) 
		for i in c])
	plt.plot(c,d, color = limitcolor, label = "90th percentile", zorder = 100)
	plt.ylim([10**-3,10**8])
	plt.ylabel(r"$\Psi$")
	plt.xlabel(r"G [mag]")

def plot_DR2_match(DR2_BGS,good,excluded,dist_limit,dr2_random = None):
	fig = plt.figure("DR2")
	fig.clf()
	if dr2_random is not None:
		semilogx(dr2_random["angular_distance"], \
			np.abs(dr2_random["magnitude_difference"]), ms = 0.1,\
			color = grey,label = "random sample")
	loglog(DR2_BGS[good]["angular_distance"], \
		np.abs(DR2_BGS[good]["magnitude_difference"]), ms = 0.5, \
		color = blue,alpha = 0.5, label = "good matches")
	loglog(DR2_BGS[good == False]["angular_distance"], \
		np.abs(DR2_BGS[good == False]["magnitude_difference"]), \
		ms = 0.5, color = green, label = "mismatches")
	loglog(DR2_BGS[excluded]["angular_distance"], \
		np.abs(DR2_BGS[excluded]["magnitude_difference"]), ms = 0.5, \
		color = red, label = "excluded sources")
	xlim = np.array([plt.xlim()[0]*5,dist_limit])
	f =lambda x: 0.3 * pow(x,0.2)
	plt.loglog(xlim, f(xlim),limitcolor, linewidth = 1)
	plt.loglog([400,400],[2e-4,10],limitcolor, linewidth = 1)

	plt.ylim([1e-4,30])
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	plt.xlabel("$\Delta\phi$ [mas]")
	plt.ylabel("$\Delta G$ [mag]")
	if ion: plt.pause(0.1)

def plot_DR2_pm(d3,dr2,pp):
	fig = plt.figure("DR2_pm", figsize = [6.4,3])
	fig.clf()
	fig.subplots_adjust(right = 0.87, bottom = 0.17)
	ax1 = plt.subplot(121)
	try:five = d3['parallax'].mask == False
	except:five= np.isnan(d3['parallax']) == False
	try:twodr2 = dr2['parallax'].mask
	except: twodr2 = np.isnan(dr2['parallax'])
	if all(five):
		twodr2 = dr2['parallax'] == 1e20
	if any(twodr2) == False:
		twodr2 = dr2['parallax'] == 1e20
	fivedr2	= twodr2==False


	plot(d3[five & twodr2]["pmra"], pp[five & twodr2][:,0], \
		color = blue, label = "2par in dr2")
	plot(d3[fivedr2]["pmra"], pp[fivedr2][:,0], color = green,\
		label = "5par in dr2")

	ax2 = plt.subplot(122)
	plot(d3[five & twodr2]["pmdec"], pp[five & twodr2][:,1], \
		color = blue ,label = "2par in dr2")
	plot(d3[fivedr2]["pmdec"], pp[fivedr2][:,1], color = green,\
		label = "5par in dr2")
	
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	ax1.set_xlabel(r'$\mu_{\alpha^{\star},DR3}$ [mas/yr]')
	ax1.set_ylabel(r'$\Delta \phi_{\rm \times2,\,\alpha^{\star}}$ [mas]')
	ax2.set_xlabel(r'$\mu_{\delta,DR3}$ [mas/yr]')
	ax2.set_ylabel(r'$\Delta \phi_{\rm \times2,\,\delta}$ [mas]')
	ax1.set_xlim([-1700,1700])
	ax1.set_ylim([-1700,1700])
	ax2.set_xlim([-1700,1700])
	ax2.set_ylim([-1700,1700])

	leg = ax1.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.1)

def plot_pos_err(BGS=None,F_pos=None,limit=None,data=None): 
	fig = plt.figure("pos_err")
	if BGS is not None:
		fig.clf()
		ax = plt.axes()
		semilogy(BGS['phot_g_mean_mag'][F_pos==False], \
			np.sqrt(BGS['ra_error'][F_pos==False] \
			* BGS['ra_error'][F_pos==False] \
			+ BGS['dec_error'][F_pos==False] \
			* BGS['dec_error'][F_pos==False]), ms = 0.5, color = red, \
			label = "excluded BGS" )
		semilogy(BGS['phot_g_mean_mag'][F_pos], \
			np.sqrt(BGS['ra_error'][F_pos] * BGS['ra_error'][F_pos] + \
			BGS['dec_error'][F_pos] * BGS['dec_error'][F_pos]), ms = 0.5, \
			color = blue, label = \
			r"$\sqrt{\sigma^{2}_{\alpha} + \sigma^{2}_{\delta}} < $" \
			+ "%i mas"%limit)
		xlim=plt.xlim()
		plt.plot(xlim,[10,10], color = limitcolor, 
			label = r'$\sigma_{\rm pos} = 10$',zorder = 5)
		# plt.xlim(xlim)
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		plt.xlabel("G [mag]")
		plt.ylabel(r"$\sqrt{\sigma_{\alpha}^{2} + \sigma_{\delta}^{2}}$ [mas]")
		ax.yaxis.set_major_formatter(formatter)

	if data is not None:
		semilogy(data['phot_g_mean_mag'], \
			np.sqrt(data['ra_error'] * data['ra_error'] + \
			data['dec_error'] * data['dec_error']),  \
			color = green, ms = 0.5, zorder = 2, \
			label = 'good BGS')
		ax = plt.gca()
		minor_locator = AutoMinorLocator(5)
		ax.xaxis.set_minor_locator(minor_locator)
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.1)

# Plot functions used in  good_HPMS.py 

def plot_HPMS(HPMS,good):
	fig = plt.figure("HPMS")
	phot = HPMS['phot_g_mean_mag']<21
	fig.clf()
	pmtot = np.sqrt(HPMS['pmra'] * HPMS['pmra'] + HPMS['pmdec']*HPMS['pmdec'])
	plt.loglog(pmtot[phot&(good==False)],\
	 	HPMS[phot&(good==False)]['parallax'],'o', color =red, ms = 0.5, 
	 	label = 'excluded high-proper-motion-stars')
	loglog(pmtot[phot&good], \
		HPMS[phot&good]['parallax'], color =blue,	ms = 0.5, 
		label = "high-proper-motion-stars")
	plt.xlabel(r"$\mu_{tot}$ [mas/y]")
	plt.ylabel(r"$\varpi$ [mas]")
	leg = plt.legend()
	for lh in leg.legendHandles: 
	 	lh._legmarker.set_alpha(1)
	 	lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.1)

def plot_ruwe(HPMS,out):
	phot  = HPMS['phot_g_mean_mag']<21
	fig = plt.figure("ruwe")
	fig.clf()
	ax = plt.axes()
	semilogy(HPMS['phot_g_mean_mag'][phot],\
		HPMS['ruwe'][phot], color =red, \
		ms = 0.5, label = "excluded HPMS")
	semilogy(HPMS[out& phot]['phot_g_mean_mag'], \
		HPMS[out& phot]['ruwe'], color =blue,\
		ms = 0.5, label = "HPMS")
	xlim=plt.xlim()
	plt.plot(xlim,[2,2], color = 'orange', label = 'ruwe = 2')
	# plt.xlim(xlim)
	plt.ylabel("Ruwe")
	plt.xlabel("G [mag]")
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	ax.yaxis.set_major_formatter(formatter)
	if ion: plt.pause(0.1)

def plot_sig_flux(HPMS,HPMS_in_GCNS_reject_bad,out,power,limit):	
	fig = plt.figure("n_vs_sig_g_flux")
	fig.clf()
	ax = plt.axes()
	f2 = out == False
	phot2 = HPMS_in_GCNS_reject_bad['phot_g_mean_mag']<21
	phot = HPMS['phot_g_mean_mag']<21
	loglog(HPMS[f2&phot]['phot_g_n_obs'],\
		HPMS[f2&phot]['phot_g_mean_flux_over_error'], color =red, \
		ms = 0.5, label = 'excluded high-proper-motion-stars')
	loglog(HPMS[out&phot]['phot_g_n_obs'], \
		HPMS[out&phot]['phot_g_mean_flux_over_error'], color =blue,\
		ms = 0.5, label = "high-proper-motion-stars")
	plt.loglog(HPMS_in_GCNS_reject_bad[phot2]['phot_g_n_obs'], \
		HPMS_in_GCNS_reject_bad[phot2]['phot_g_mean_flux_over_error'], "s", \
		color ="k" , ms = 2, label = "GCNS probability < 0.38")
	plt.ylabel(r"$G_{flux}/\sigma_{G_{flux}}$")
	plt.xlabel(r"$n_{obs}$")
	xlim = plt.xlim()
	ylim = plt.ylim()

	s1 = "%i"%power if power%1==0 else "%.1f"%power
	if limit > 1000: 
		kk = "%.2e"%(limit)
		if kk[2:4]=="00":
			s2 = '%s \cdot 10^{%i}'%(kk[0],int(kk.split('e')[-1]))
		else: 
			s2 = '%s \cdot 10^{%i}'%(kk[0:4],int(kk.split('e')[-1]))
	else:
		s2 = "%i"%limit 
	plt.plot([np.power(limit/ylim[0],1/power),
		np.power(limit/ylim[1],1/power)],ylim, color = limitcolor,
		label = r"$G_{flux}/\sigma_{G_{flux}}\cdot n_{obs}^{%s} < %s$"%(s1,s2)) 
	leg = plt.legend()
	for lh in leg.legendHandles:  
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	ax.xaxis.set_major_formatter(formatter)
	ax.yaxis.set_major_formatter(formatter)
	if ion: plt.pause(0.1)

# Plot functions used in good_Pairs.py

def plot_delta_pm_ax1(rp,FF):
	fig = plt.figure("delta_pm", figsize = [6.4,3])
	fig.clf()
	fig.subplots_adjust(right = 0.87, bottom = 0.17)
	ax1 = plt.subplot(121)
	plot(rp['pmra'] - rp['ob_pmra'],\
		rp['pmdec'] - rp['ob_pmdec'], color = red,  ms = 0.5,\
		label = 'excluded Pairs')
	plot(rp['pmra'][FF] - rp['ob_pmra'][FF], \
		rp['pmdec'][FF] - rp['ob_pmdec'][FF], color = blue, ms = 0.5,\
		label = 'Pairs')
	plt.xlabel(r"$\mu_{\alpha\star} - Sou\_\mu_{\alpha\star}$ [mas/yr]")
	plt.ylabel(r"$\mu_{\delta} - Sou\_\mu_{\delta}$ [mas/yr]")
	plt.xlim([-500,500])
	plt.ylim([-500,500])
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)	
	if ion: plt.pause(0.1)

def  plot_delta_pm_ax2(rp,F_pm_3,FF):
	fig = plt.figure("delta_pm", figsize = [6.4,3])
	fig.subplots_adjust(right = 0.87, bottom = 0.17)
	ax2 = plt.subplot(122)
	plot(rp[F_pm_3]['pmra'] \
		- rp[F_pm_3]['ob_displacement_ra_doubled'],\
		rp[F_pm_3]['pmdec'] \
		- rp[F_pm_3]['ob_displacement_dec_doubled'],\
		color = red,  ms = 0.5)
	plot(rp['pmra'][FF] - rp['ob_displacement_ra_doubled'][FF], \
		rp['pmdec'][FF] - rp['ob_displacement_dec_doubled'][FF], \
		color = blue, ms = 0.5)
	plt.xlabel(r'$\mu_{\alpha^{\star}} - \Delta \phi_{\rm \times2,\,\alpha^{\star}}/1{\rm yr}$ [mas/yr]')
	plt.ylabel(r'$\mu_{\delta} - \Delta \phi_{\rm \times2,\,\delta}/1{\rm yr}$ [mas/yr]')
	plt.xlim([-500,500])
	plt.ylim([-500,500])
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	if ion: plt.pause(0.01)

def plot_sim_px(rp= None,F_px_1=None, data = None):
	fig = plt.figure("sim_px")
	if rp is not None:
		semilogx(rp['parallax'], \
			rp['ob_parallax'],  color = red, label='excluded Pairs')
		semilogx(rp[F_px_1]['parallax'], \
			rp[F_px_1]['ob_parallax'], color = blue, label='Pairs')
		plt.xlabel(r'$\varpi$ [mas]')
		plt.ylabel(r'$Sou\_\varpi$ [mas]')
	if data is not None:
		semilogx(data['parallax'], \
			data['ob_parallax'],  color = yellow, label='Good Candidates')
		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
		plt.ylim([-10,20])
		ax = plt.gca()
		minor_locator = AutoMinorLocator(5)
		ax.yaxis.set_minor_locator(minor_locator)
	if ion: plt.pause(0.01)

# Plot functions used in raw_data
def plot_psi_part2(data):
	plt.figure("psi")
	semilogy(data['phot_g_mean_mag'], data["psi"], \
		color = yellow, ms = 1,  zorder = 2, \
		label = 'good BGS')
		#color = green, ms = 0.5,  zorder = 2
	ax = plt.gca()
	plt.yticks(10**np.arange(*np.log10(plt.yticks()[0])[[0,-1]]+[0,1]))
	y_minor = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0), 
		numticks = 10)
	ax.yaxis.set_minor_locator(y_minor)
	minor_locator = AutoMinorLocator(5)
	ax.xaxis.set_minor_locator(minor_locator)
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.01)

def plot_sim_px_part2(data):
	fig = plt.figure("sim_px")
	semilogx(data['parallax'], \
		data['ob_parallax'],  color = yellow)
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	plt.ylim([-10,20])
	ax = plt.gca()
	minor_locator = AutoMinorLocator(5)
	ax.yaxis.set_minor_locator(minor_locator)
	if ion: plt.pause(0.01)

def plot_HPMS_part2(data):
	plt.figure("HPMS")
	cc = data['ob_parallax']>0
	xlim = plt.xlim()
	ylim = plt.ylim()
	loglog(np.sqrt(data[cc]['ob_pmra']**2+data[cc]['ob_pmdec']**2),
		data[cc]['ob_parallax'],color = green, ms = 0.5, 
		label = 'background stars')
	plt.xlim([1,xlim[1]])
	plt.ylim([0.1,ylim[1]])
	leg = plt.legend()
	for lh in leg.legendHandles: 
		lh._legmarker.set_alpha(1)
		lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.01)


def plot_CMD(tab = None, g_rp= None,B_abs=None,WD_Terms=None,RG_Terms=None):
	if g_rp is not None:
		plt.figure('CMD')
		gg = np.where(tab['phot_bp_mean_mag'] != 0)
		plot(g_rp[gg], B_abs[gg],\
				'o', color = grey, ms = 0.5, label =  'Candidates', zorder = 0)
		xlim = [-0.5,2.5]
		ylim = [24,-2]
		x = np.linspace(*xlim, 1000)
		plt.plot(x,polyval(x,WD_Terms), color=green, label='WD limit',zorder=4)
		plt.plot(x,polyval(x,RG_Terms), color=red, label = 'RG limit', zorder=3)
		plt.xlim(xlim)
		plt.ylim(ylim)
		plt.xlabel(r'$G_{RP} - G$ [mag]')
		plt.ylabel(r'$G_{BP,abs}$ [mag]')
	else:
		plt.figure('CMD')
		#select only events with full photometry
		gg = np.where(tab['phot_bp_mean_mag'] != 0) 

		plt.plot(tab['phot_g_mean_mag'] - tab['phot_rp_mean_mag'],\
			tab['phot_bp_mean_mag'] + 5*np.log10(tab['parallax']/100),\
			'o', color = blue, ms = 0.5, label = 'predicted events', \
			zorder = 2)

		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
	if ion: plt.pause(0.01)


# used for Results

def plot_psi_result_part_2(BGS,result,blacklist, filtered):
	if 'psi_result' in plt.get_figlabels():
		plt.figure('psi_result')

		bl = np.isin(BGS['source_id'],\
			result[blacklist]['ob_source_id'])
		bl2= bl & (BGS['dec'] < -30)

		res = np.isin(BGS['source_id'],result[filtered]['ob_source_id'])
		g = BGS['phot_g_mean_mag'][res]
		g_BL = BGS['phot_g_mean_mag'][bl]
		g_BL2 = BGS['phot_g_mean_mag'][bl2]

		if 'psi' in BGS.colnames:				
			psi = BGS['psi'][res]
			psi_BL = BGS['psi'][bl]
			psi_BL2 = BGS['psi'][bl2]
		else:
			gamma = np.maximum(pow(10, 0.2 * (BGS['phot_g_mean_mag'] - 18)), 1)
			psi_bgs = BGS['astrometric_sigma5d_max'] / (1.2 * gamma)
			psi = psi_bgs[res]
			psi_BL = psi_bgs[bl]
			psi_BL2 = psi_bgs[bl2]
		plt.plot(g,psi,'o',ms = 1, color = blue,\
			label = r'predicted events', zorder = 1)	
		plt.plot(g_BL,psi_BL,'s',ms = 1, color = red, zorder = 2,\
			label = r'removed candidates  $(\delta >-30^{\circ})$')		
		plt.plot(g_BL2,psi_BL2,'s',ms = 1, color = 'm', zorder = 3,\
			label = r'removed candidates  $(\delta <-30^{\circ})$')
		ax = plt.gca()
		plt.yticks(10**np.arange(*np.log10(plt.yticks()[0])[[0,-1]]+[0,1]))
		y_minor = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0),
			numticks = 10)
		ax.yaxis.set_minor_locator(y_minor)
		minor_locator = AutoMinorLocator(5)
		ax.xaxis.set_minor_locator(minor_locator)

		leg = plt.legend()
		for lh in leg.legendHandles: 
			lh._legmarker.set_alpha(1)
			lh._legmarker.set_markersize(5)
	else: print('No plot found for "psi_result"')



def plot_results(result):
	for typ in ['All','WD','RG','MS','BD']:
		for timerange in [[2010,2070],[2014,2032]]:
			for mag_dif in [20,6,3]:
				if mag_dif == 20 : add =''
				else: add = '_DG_'+str(mag_dif)
				if typ == 'All':
					if timerange[1] == 2070:
						fig = plt.figure('Results_'+typ+add,figsize= [12.8,4.8])
						plt.subplots_adjust(\
							left = plt.rcParams['figure.subplot.left']/2,
							right= plt.rcParams['figure.subplot.right']/2+0.5)
					else: fig = plt.figure('Results_2030_'+typ+add)
					which = np.ones(len(result),bool)
				elif typ == 'MS':
					if timerange[1] == 2070:
						fig = plt.figure('Results_'+typ+add, figsize=[12.8,4.8])
						plt.subplots_adjust(\
						left = plt.rcParams['figure.subplot.left']/2,
						right= plt.rcParams['figure.subplot.right']/2+0.5)
					else:fig = plt.figure('Results_2030_'+typ+add)
					which = result['star_type']==typ
				else:
					if timerange[1] == 2070: 
						fig = plt.figure('Results_'+typ+add)
					else: fig = plt.figure('Results_2030_'+typ+add)
					which = result['star_type']==typ
				which = which & ((result['ob_phot_g_mean_mag']
					-result['phot_g_mean_mag']) < mag_dif)
				ax = plt.gca()
				twopar = which & (result['ob_parallax'] == 0)
				fivepar = which & (result['ob_parallax'] != 0)
				plt.semilogy(result[twopar]['TCA'],result[twopar]['shift_plus'],
					'.', color = grey, zorder = -5, 
					label = '2 parameter solution')
				plt.errorbar(result[fivepar]['TCA'],
					result[fivepar]['shift_plus'],fmt = 'o',\
					xerr = result[fivepar]['TCA_error'],
					yerr = [-result[fivepar]['shift_plus_error_m'],\
					result[fivepar]['shift_plus_error_p']],\
					color = 'blue', ms = 2,lw = 0.5, markeredgewidth = 0.2,\
					markeredgecolor = 'k', label = '5 parameter solution')
				plt.ylim([0.05,50])
				plt.yticks([0.1,0.2,0.5,1,2,5,10,20,50],\
					['0.10','0.20','0.50','1.0','2.0','5.0','10','20','50'])
				plt.xlim(timerange)

				leg = plt.legend()
				for lh in leg.legendHandles: 
					try: 
						lh._legmarker.set_alpha(1)
						lh._legmarker.set_markersize(5)
					except: pass
				plt.ylabel(r'$\delta\theta_{+}$ [mas]')
				plt.xlabel(r'$T_{CA}$ [yr]')
				if ion: plt.pause(0.01)
				if mag_dif == 3: 
					c = red 
					z = 2
					label = r'$\Delta G < 3$ mag'
				if mag_dif == 6: 
					z = 1
					c = blue
					label = r'$\Delta G < 6$ mag'
				if mag_dif == 20: 
					z = 0
					c = grey 
					label = 'All events'
				if timerange[1] == 2070:
					if typ == 'All' or typ == 'MS':
						fig = plt.figure('Results_Hist_'+ typ, \
							figsize = [12.8,2.4])
						plt.subplots_adjust(\
							left = plt.rcParams['figure.subplot.left']/2,
							right= plt.rcParams['figure.subplot.right']/2+0.5,
							bottom = plt.rcParams['figure.subplot.bottom']*2)
					else: 
						fig = plt.figure('Results_Hist_' + typ,\
						figsize = [6.4,2.4])
						plt.subplots_adjust(\
							bottom = plt.rcParams['figure.subplot.bottom']*2)
					plt.hist(result['TCA'][which],color = c, bins = 60, \
						range= timerange, rwidth = 0.9, zorder = z, 
						label = label)
				else:
					fig = plt.figure('Results_Hist_2030_' + typ, 
						figsize = [6.4,2.4])
					plt.subplots_adjust(\
							bottom = plt.rcParams['figure.subplot.bottom']*2)
					plt.hist(result['TCA'][which], color = c, range= timerange,
						bins = (timerange[1]-timerange[0])*2,rwidth = 0.9,
						zorder = z, weights = np.ones(len(result[which]))*2,
						label = label)
				leg = plt.legend()

				plt.xlim(timerange)
				plt.ylabel(r'Events per year')
				plt.xlabel(r'$T_{CA}$ [yr]')

	# plot magnification 
	for typ in ['All','WD','RG','MS','BD']:	
		if typ == 'All':
			which = np.ones(len(result),bool)
		else:
			which = result['star_type']==typ

		fig = plt.figure('Magnification_'+typ)
		which = which & (result['magnification']> 0.001)
		ax = plt.gca()
		twopar = which & (result['ob_parallax'] == 0)
		fivepar = which & (result['ob_parallax'] != 0)
		plt.semilogy(result[twopar]['TCA'],result[twopar]['magnification'],'.', 
			color = 'grey',zorder = -5,	label = '2 parameter solution')
		plt.errorbar(result[fivepar]['TCA'],result[fivepar]['magnification'],
			fmt = 'o',xerr = result[fivepar]['TCA_error'],
			yerr = [np.minimum(-result[fivepar]['magnification_error_m'],
			result[fivepar]['magnification']-0.0001),
			result[fivepar]['magnification_error_p']],
			color = 'blue', ms = 2,lw = 0.5, markeredgewidth = 0.2,
			markeredgecolor = 'k', label = '5 parameter solution')
		plt.ylim([0.001,5])
		plt.yticks([0.001,0.005,0.01,0.05,0.1,0.5,1,5],\
			['0.001', '0.005','0.01','0.05','0.1','0.5','1','5'])
		plt.xlim(2010,2070)
		leg = plt.legend()
		for lh in leg.legendHandles: 
			try: 
				lh._legmarker.set_alpha(1)
				lh._legmarker.set_markersize(5)
			except: pass
		plt.ylabel(r'$\Delta m$ [mag]')
		plt.xlabel(r'$T_{CA}$ [yr]')
		if ion: plt.pause(0.01)

def save_figures(Folder):
	for im_name in plt.get_figlabels(): # loop over every image
		f = plt.figure(im_name)
		prefix = setup.prefix
		print('Save Figure: '+Folder+'Images/'+im_name+prefix+'.png')
		f.savefig(Folder+'Images/'+im_name+prefix+'.png', transparent =False)
		#print('Save Figure: '+Folder+'Images/'+im_name+prefix+'.pdf')
		#f.savefig(Folder+'Images/'+im_name+prefix+'.pdf', transparent =False)
