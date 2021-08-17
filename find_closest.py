import numpy as np
import time
from astropy.table import MaskedColumn 
#sub packages
from setup import Gaia_epoch
from astrometry import movePm_parallax, getGCDist_3Vec,epoch_predef,\
	pos_sun_predef,pos_sun, Vector3,dirVecToCelCoos,parallax_correction
from utils import cosdeg, sindeg

def find_Minimum(f, left, right, minInterval = 3e-10, nmax = 1e6):
	'''
	Nested intervals algorithm to find minima of a function in a 
	certain interval
	'''
	r = 0.61803399
	c = 1-r
	x0 = left
	x3 = right
	x1 = c*x3 + r*x0
	x2 = r*x3 + c*x0
	f1 = f(x1)
	f2 = f(x2)
	n = 0
	while np.fabs(x3-x0) > minInterval:#*(np.fabs(x1)+np.fabs(x2)):
		n += 1
		if f2 < f1:
			x0 = x1
			x1 = x2
			x2 = c * x0 + r * x3
			f1 = f2
			f2 = f(x2)
		if f1 <= f2:
			x3 = x2
			x2 = x1
			x1 = c * x3 + r * x0
			f2 = f1
			f1 = f(x1)
		if n> nmax: print(n_max)
		break
	if f1 < f2:
		return x1
	else:
		return x2


def calc_approx_date(tab):
	cd = cosdeg(tab['dec'])
	delta_ra = (tab['ob_ra'] - tab['ra']) * cd * 3.6e6
	delta_dec = (tab['ob_dec'] - tab['dec']) * 3.6e6

	obj_pmra = tab['ob_pmra'].copy()
	obj_pmdec = tab['ob_pmdec'].copy()

	# set DR2 displacement as guess if ob_pm not known 
	if 'ob_displacement_ra_doubled' in tab.colnames:
		obj_pmra[obj_pmra==0] = tab[obj_pmra==0]['ob_displacement_ra_doubled']
		obj_pmdec[obj_pmdec==0] = \
			tab[obj_pmdec==0]['ob_displacement_dec_doubled']

	delta_pmra = (tab["pmra"]-obj_pmra)
	delta_pmdec = (tab["pmdec"]-obj_pmdec)
	delta_pmtot = np.sqrt(delta_pmdec * delta_pmdec + delta_pmra * delta_pmra)
	minDate = (delta_ra * delta_pmra + delta_dec * delta_pmdec) \
		/ delta_pmtot ** 2 + Gaia_epoch
	minDist = np.abs(delta_ra * delta_pmdec - delta_dec * delta_pmra) \
		/ delta_pmtot
	minDate =MaskedColumn(minDate,dtype = 'float64', unit = 'jyear',\
		description = 'approximated date of the closest appoach ')	
	minDist =MaskedColumn(minDist,dtype = 'float64', unit = 'mas', \
		description = 'approximated separation at the closest appoach ')	
	return minDate,minDist


def estimate_Closest_parallax(row, gaia = True):
	"""
	units: 
	Ra,DEC in Deg
	ra_err, Dec_err in mas
	pmRA, pmDec in mas/year, 
	pmRa_err,pmDec_err in mas/year 
	parrallax in mas
	parrallax_err in mas

	distance in mas
	date in year

	mass in solar_mass
	einsteinradii in mas

	shift in mas
	magnifikation in delta mag

	"""

	# Define distance function with and with out parallax


	lens_ra = row['ra']
	lens_dec = row['dec']
	lens_pmra = row['pmra']
	lens_pmdec = row['pmdec']
	lens_parallax = row['parallax']

	obj_ra = row['ob_ra']
	obj_dec = row['ob_dec']

	# set DR2 displacement as guess if ob_pm not known 
	if row['ob_pmra'] == 0 and 'ob_displacement_ra_doubled' in row.colnames:
		obj_pmra = row['ob_displacement_ra_doubled']
		obj_pmdec = row['ob_displacement_dec_doubled' ]
		mu = np.sqrt(obj_pmra*obj_pmra +obj_pmdec*obj_pmdec)
		obj_parallax = 4.740*mu/75 
	else:
		obj_pmra = row['ob_pmra']
		obj_pmdec = row['ob_pmdec']
		obj_parallax = row['ob_parallax']


	dist_pre = lambda i: getGCDist_3Vec(
		movePm_parallax(lens_ra, lens_dec, lens_pmra, lens_pmdec, \
			lens_parallax, epoch_predef[i], pos_sun_predef[i]), \
		movePm_parallax(obj_ra, obj_dec, obj_pmra, obj_pmdec, obj_parallax, \
			epoch_predef[i],pos_sun_predef[i])) * 3.6e6
	dist2 = lambda t: getGCDist_3Vec(
		movePm_parallax(lens_ra, lens_dec, lens_pmra, lens_pmdec, \
			lens_parallax, t), \
		movePm_parallax(obj_ra, obj_dec, obj_pmra, obj_pmdec,\
			obj_parallax,t)) * 3.6e6


	tt = []
	tt.append(time.time())
	Mu_tot = np.sqrt((lens_pmra - obj_pmra)**2 + (lens_pmdec - obj_pmdec)**2)
	#finde	minima limit to approx 1 week 
	t_0 = row['approx_tca']
	tlim = max(2.0, (3 * lens_parallax) / Mu_tot)
	t_min = t_0 - tlim - Gaia_epoch
	t_max = t_0 + tlim - Gaia_epoch
	
	tt.append(time.time())
	nn = np.where((epoch_predef > t_min) & (epoch_predef < t_max))[0]
	dist = list(map(dist_pre, nn))
	number_list = range(len(nn) - 2)
	a = list(filter(lambda x: dist[x+1] - dist[x] < 0 \
		and dist[x+1] - dist[x+2] <0 , number_list))
	tt.append(time.time())
	#evaluate all minima
	closest_time_vec = [999999.,999999.,999999.]
	closest_dist_vec = [999999.,999999.,999999.]
	tt44 = 0
	tt55 = 0
	for j in nn[a]:
		t1 = epoch_predef[j-1] 
		t2 = epoch_predef[j+3]
		tt44 += time.time()
		closest_time = find_Minimum(dist2, t1, t2)
		tt55 += time.time()
		closest_dist = dist2(closest_time)
		if closest_dist <= closest_dist_vec[0]:
			closest_time_vec[2] = closest_time_vec[1]
			closest_time_vec[1] = closest_time_vec[0]
			closest_time_vec[0] = closest_time
			closest_dist_vec[2] = closest_dist_vec[1]
			closest_dist_vec[1] = closest_dist_vec[0]
			closest_dist_vec[0] = closest_dist		

		elif closest_dist <= closest_dist_vec[1]:
			closest_time_vec[2] = closest_time_vec[1]					
			closest_time_vec[1] = closest_time
			closest_dist_vec[2] = closest_dist_vec[1]
			closest_dist_vec[1] = closest_dist		

		elif closest_dist <= closest_dist_vec[2]:
			closest_time_vec[2] = closest_time
			closest_dist_vec[2] = closest_dist	
	tt.append(time.time())

	closest_dist_vec = [x if x < 999999 else -999 for x in closest_dist_vec]
	closest_date_vec = [x + Gaia_epoch if x < 999999 else -999 \
			for x in closest_time_vec]


	if not gaia:
		tt.append(time.time())
		tt = np.array(tt)
		cpt = tt[1:] - tt[:-1]
		cpt = np.append(cpt, (tt55-tt44) / len(a))
		return np.array([closest_date_vec, closest_dist_vec]), cpt

	else:
		"""
		calculate closest approxhes seen from L2 therefore the parallax is
		 multiplied with 1.01
		"""		

		dist3 = lambda t: getGCDist_3Vec(
			movePm_parallax(lens_ra, lens_dec, lens_pmra, lens_pmdec, \
			lens_parallax, t, gaia = True), \
			movePm_parallax(obj_ra, obj_dec, obj_pmra, obj_pmdec, \
			obj_parallax, t, gaia = True)) * 3.6e6	
		dist_pre = lambda i: getGCDist_3Vec(
			movePm_parallax(lens_ra, lens_dec, lens_pmra, lens_pmdec, \
			lens_parallax, epoch_predef[i], pos_sun_predef[i], gaia = True),\
			movePm_parallax(obj_ra, obj_dec, obj_pmra, obj_pmdec, \
			obj_parallax, epoch_predef[i], pos_sun_predef[i], gaia = True)) \
			* 3.6e6

		dist = list(map(dist_pre,nn))
		number_list = range(len(nn) - 2)
		a2 = list(filter(lambda x: dist[x+1] - dist[x] < 0 \
			and dist[x+1] - dist[x+2] < 0 , number_list))

		tt.append(time.time())
		closest_time_gaia_vec = [999999.,999999.,999999.]
		closest_dist_gaia_vec = [999999.,999999.,999999.]
		for j in nn[a2]:
			t1 = epoch_predef[j+0] 
			t2 = epoch_predef[j+2]
			tt44 += time.time()
			closest_time = find_Minimum(dist3,t1,t2)
			tt55 += time.time()
			closest_dist = dist3(closest_time)
			if closest_dist <= closest_dist_gaia_vec[0]:
				closest_time_gaia_vec[2] = closest_time_gaia_vec[1]
				closest_time_gaia_vec[1] = closest_time_gaia_vec[0]
				closest_time_gaia_vec[0] = closest_time
				closest_dist_gaia_vec[2] = closest_dist_gaia_vec[1]
				closest_dist_gaia_vec[1] = closest_dist_gaia_vec[0]
				closest_dist_gaia_vec[0] = closest_dist	

			elif closest_dist <= closest_dist_gaia_vec[1]:
				closest_time_gaia_vec[2] = closest_time_gaia_vec[1]
				closest_time_gaia_vec[1] = closest_time
				closest_dist_gaia_vec[2] = closest_dist_gaia_vec[1]
				closest_dist_gaia_vec[1] = closest_dist	

			elif closest_dist <= closest_dist_gaia_vec[2]:
				closest_time_gaia_vec[2] = closest_time
				closest_dist_gaia_vec[2] = closest_dist
		
		closest_dist_gaia_vec = [x if x < 999999 else -999 \
			for x in closest_dist_gaia_vec]
		closest_date_gaia_vec = [x + Gaia_epoch if x < 999999 else -999 \
			for x in closest_time_gaia_vec]
		tt.append(time.time())
		tt = np.array(tt)
		cpt = tt[1:] - tt[:-1]
		cpt = np.append(cpt, (tt55-tt44) / (len(a)+len(a2)))
		return np.array([closest_date_vec, closest_dist_vec, \
			closest_date_gaia_vec, closest_dist_gaia_vec]), cpt


def estimate_errors_parallax(tab, delta_t_approx = 1/26.,gaia = False):
	'''
	_error estimation considering parallax also
	'''

	cd, sd = cosdeg(tab["dec"]), sindeg(tab["dec"])
	ca, sa = cosdeg(tab["ra"]), sindeg(tab["ra"])
	ocd, osd = cosdeg(tab["ob_dec"]), sindeg(tab["ob_dec"])
	oca, osa = cosdeg(tab["ob_ra"]), sindeg(tab["ob_ra"])



	obj_pmra = tab['ob_pmra'].copy()
	obj_pmdec = tab['ob_pmdec'].copy()
	obj_parallax = tab["ob_parallax"].copy()
	# set DR2 displacement as guess if ob_pm not known 
	if 'ob_displacement_ra_doubled' in tab.colnames:
		obj_pmra[obj_pmra==0] = tab[obj_pmra==0]['ob_displacement_ra_doubled']
		obj_pmdec[obj_pmdec==0] = \
			tab[obj_pmdec==0]['ob_displacement_dec_doubled']
		mu = obj_pmra[obj_parallax==0]**2+obj_pmra[obj_parallax==0]**2
		obj_parallax[obj_parallax==0] = 4.740*mu/75 

	lens_parallax = tab["parallax"]
	obj_parallax_err = tab["ob_parallax_error"]
	lens_parallax_err = tab["parallax_error"]
	delta_pmra = (tab["pmra"]-obj_pmra)
	delta_pmdec = (tab["pmdec"]-obj_pmdec)
	lens_ra_err = (tab["ra_error"])
	lens_dec_err = (tab["dec_error"])
	obj_ra_err = (tab["ob_ra_error"])
	obj_dec_err = (tab["ob_dec_error"])
	lens_pmra_err = (tab["pmra_error"])
	lens_pmdec_err = (tab["pmdec_error"])
	obj_pmra_err = (tab["ob_pmra_error"])
	obj_pmdec_err = (tab["ob_pmdec_error"])


	#!!!! change to exact value
	if gaia:	minDate = tab['L2_TCA']
	else:		minDate = tab['TCA']
	obj_baryVec = Vector3(ocd*oca, ocd*osa, osd)
	lens_baryVec = Vector3(cd*ca, cd*sa, sd)
	earthVec = pos_sun(minDate.data-Gaia_epoch)

	lens_ra,lens_dec = dirVecToCelCoos(parallax_correction(
		lens_baryVec,earthVec,lens_parallax * (gaia * .01 + 1)))
	obj_ra,obj_dec = dirVecToCelCoos(parallax_correction(
		obj_baryVec,earthVec,obj_parallax * (gaia * .01 + 1)))

	delta_ra = (lens_ra-obj_ra) * cd * 3.6e6
	delta_dec = (lens_dec-obj_dec) * 3.6e6



	delta_pmtot = np.sqrt(np.square(delta_pmra) + np.square(delta_pmdec))
	delta_ra_err = np.sqrt(np.square(lens_ra_err) + np.square(obj_ra_err))
	delta_dec_err = np.sqrt(np.square(lens_dec_err) + np.square(obj_dec_err))
	delta_pmra_err = np.sqrt(np.square(lens_pmra_err) \
		+ np.square(obj_pmra_err))
	delta_pmdec_err = np.sqrt(np.square(lens_pmdec_err) \
		+ np.square(obj_pmdec_err)) 



	minDate_err = np.sqrt(\
		pow(pow(delta_pmdec,2) * delta_ra \
		- 2 * delta_pmra * delta_pmdec * delta_dec \
		- pow(delta_pmra,2) * delta_ra,2) \
		/ pow(delta_pmtot,8) * pow(delta_pmra_err,2) \
		+ pow(pow(delta_pmra,2) * delta_dec \
		- 2 * delta_pmra * delta_pmdec * delta_ra \
		- pow(delta_pmdec,2) * delta_dec,2) \
		/ pow(delta_pmtot,8) * pow(delta_pmdec_err,2) \
		+ pow(delta_ra_err * delta_pmra,2) / pow(delta_pmtot,4) \
		+ pow(delta_dec_err * delta_pmdec,2) / pow(delta_pmtot,4))


	minDate_err2 = np.minimum(minDate_err,0.25)

	earthVec2 = pos_sun(minDate.data - Gaia_epoch + minDate_err2.data)

	#relative parallax_error
	parallax_err = np.sqrt(np.square(lens_parallax_err) \
		+ np.square(obj_parallax_err))

	# calculate position at silightly different parallax and epoch, unit = deg
	lens_ra2, lens_dec2 = dirVecToCelCoos(parallax_correction(lens_baryVec, \
		earthVec2, (lens_parallax - obj_parallax) * (gaia * .01 + 1)))
	lens_ra3, lens_dec3 = dirVecToCelCoos(parallax_correction(lens_baryVec, \
		earthVec, (lens_parallax - obj_parallax + parallax_err) * \
		(gaia * .01 + 1)))


	dvra = (lens_ra - lens_ra2) * 3.6e6 * cd / minDate_err2 * minDate_err 
	dvdec = (lens_dec - lens_dec2)*3.6e6 / minDate_err2 * minDate_err 
	dpx_ra = (lens_ra - lens_ra3) * 3.6e6 * cd 
	dpx_dec = (lens_dec -lens_dec3) * 3.6e6


	delta_ra_err = np.sqrt(delta_ra_err**2 + dpx_ra**2 + dvra ** 2)
	delta_dec_err = np.sqrt(delta_dec_err**2 + dpx_dec**2 + dvdec ** 2)
	'''
	minDate = (delta_ra * delta_pmra + delta_dec * delta_pmdec) 
		/ delta_pmtot ** 2
	minDist = (delta_ra * delta_pmdec - delta_dec * delta_pmra ) 
		/ delta_pmtot
	'''
	minDate_err = np.sqrt(\
		pow(pow(delta_pmdec, 2) * delta_ra \
		- 2 * delta_pmra * delta_pmdec * delta_dec \
		- pow(delta_pmra, 2) * delta_ra, 2) \
		/ pow(delta_pmtot, 8) * pow(delta_pmra_err, 2) \
		+ pow(pow(delta_pmra, 2) * delta_dec \
		- 2 * delta_pmra * delta_pmdec * delta_ra \
		- pow(delta_pmdec, 2) * delta_dec, 2) \
		/ pow(delta_pmtot, 8) * pow(delta_pmdec_err, 2) \
		+ pow(delta_ra_err * delta_pmra, 2) / pow(delta_pmtot, 4) \
		+ pow(delta_dec_err * delta_pmdec, 2) / pow(delta_pmtot, 4))

	minDist_err = np.sqrt(
		pow((delta_dec * pow(delta_pmdec, 2) 
		+ delta_ra * delta_pmdec * delta_pmra) * delta_pmra_err, 2) 
		/ pow(delta_pmtot, 6) 
		+ pow((delta_ra * pow(delta_pmra, 2) 
		+ delta_dec * delta_pmra * delta_pmdec) * delta_pmdec_err, 2) 
		/ pow(delta_pmtot, 6)
		+ pow(delta_ra_err * delta_pmdec, 2) / pow(delta_pmtot, 2) 
		+ pow(delta_dec_err * delta_pmra, 2) / pow(delta_pmtot, 2) )



	if gaia: 
		minDate_err = MaskedColumn(minDate_err, dtype = 'float64', \
			unit = 'year', description = 'Error of the Time of the closest ' \
			+ 'approch for Lagrange Point 2 parallax')
		minDist_err = MaskedColumn(minDist_err, dtype = 'float64', \
			unit = 'mas', description = 'Error of the closest distance for ' \
			 + 'Lagrange Point 2 parallax')
	else: 
		minDate_err = MaskedColumn(minDate_err, dtype = 'float64', \
			unit = 'year', description = 'Error of the Time of the closest ' \
			+ 'approch')
		minDist_err = MaskedColumn(minDist_err, dtype = 'float64',\
			unit = 'mas', description = 'Error of the closest distance')
	return minDate_err, minDist_err


