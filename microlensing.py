import numpy as np
from astropy.table import MaskedColumn 
import time

const_Einsteinradius = 2.854062733172987 #(mas/Msun)**0.5

def calculate_Einsteinradius(tab): 
	"""
	calculate the Einstein radius in unit mas from Mass in SolarMasses,
	and parallax in mas 
	ThetaE[rad] = sqrt(4GM/c^2)*sqrt(d_LS/(d_L*d_S))
	ThetaE[mas] = const *sqrt(M/M_sun *((Lensparallax[mas]-Objparallax[mas])) 
	const = 180/pi * 3.6e6 *sqrt(4GM_sun / (c^2* 1pc[m]*1000))
	"""
	tt = []
	tt.append(time.time())
	ThetaE = const_Einsteinradius * np.sqrt(tab['mass']\
		* (tab['parallax'] - tab['ob_parallax']))
	ThetaE_error = ThetaE * 0.5 * np.sqrt(\
		np.square(tab['mass_error'] / tab['mass']) \
		+ (np.square(tab['parallax_error']) \
		+ np.square(tab['ob_parallax_error'])) \
		/ np.square(tab['parallax'] - tab['ob_parallax']))
	tt.append(time.time())
	ThetaE = MaskedColumn(ThetaE,dtype = 'float64', unit = 'mas', \
		description = 'Einstein radius of the event')
	ThetaE_error =MaskedColumn(ThetaE_error,dtype = 'float64', unit = 'mas',\
		description = 'Error of the Einstein radius')	
	tt.append(time.time())
	tt = np.array(tt)
	cpt = tt[1:] - tt[:-1]
	return ThetaE, ThetaE_error, cpt

def calc_shift(dist, ThetaE, dist_error=None, ThetaE_error=None, \
		FL_FS=None):
	'''
	calculate the expected shift of the center if light 
	expect dist and ThetaE in the same unit
	return shift in the same unit as dist and ThetaE
	if dist_error and ThetaE_error is given also the error is propagated 
	FL_FS is a placeholder to match with the other calc_... functions 
	'''

	u = dist / ThetaE

	# maximum shift at u = sqrt(2)
	u_2 = np.maximum(u, np.sqrt(2))

	shift = u_2 / (u_2 ** 2 + 2) * ThetaE
	if dist_error is None: return shift
	elif len(dist_error.shape) == 2:
		u_e = dist_error / ThetaE_error
		u_2_e = np.maximum(u_e, np.sqrt(2))
		shift_e = u_2_e / (u_2_e ** 2 + 2) * ThetaE_error

		shift = np.percentile(shift_e, 50, axis = 1) 
		shift_error_m = np.percentile(shift_e, 15.866, axis = 1) \
			- shift
		shift_error_p = np.percentile(shift_e, 84.866, axis = 1) \
			- shift
		return shift,shift_error_p,shift_error_m
	else:
		# calculate error
		shift_error = shift**2 \
			* np.sqrt(np.square((2 / dist**2 - 1/ThetaE**2) * dist_error)
			+ np.square(2 * dist / ThetaE**3 * ThetaE_error))
		return shift,shift_error

def calc_shift_plus(dist, ThetaE, dist_error=None, ThetaE_error=None, \
	FL_FS=None):
	'''
	calculate the expected shift of the bright image 
	expect dist and ThetaE in the same unit
	return shift in the same unit as dist and ThetaE
	if dist_error and ThetaE_error is given also the error is propagated 
	FL_FS is a placeholder to match with the other calc_... functions 
	'''
	shift_plus = (np.sqrt(dist**2 + 4 * ThetaE**2) - dist)/2

	if dist_error is None: return shift_plus
	elif len(dist_error.shape) == 2:
		shift_e = (np.sqrt(dist_error**2 + 4 * ThetaE_error**2) - dist_error)\
			/ 2
		shift_plus = 	np.percentile(shift_e, 50, axis = 1)
		shift_error_m = np.percentile(shift_e, 15.866, axis = 1) \
			- shift_plus
		shift_error_p = np.percentile(shift_e, 84.866, axis = 1) \
			- shift_plus
		return shift_plus, shift_error_p, shift_error_m
	else:
		# calculate error
		shift_plus_error = np.sqrt(np.square((dist / np.sqrt(dist**2 \
			+ 4 * ThetaE**2) - 1) / 2 * dist_error) \
			+ np.square(4 * ThetaE / 2 / np.sqrt(dist**2 + 4 * ThetaE**2) \
			* ThetaE_error))
		return shift_plus, shift_plus_error


def calc_shift_lum(dist, ThetaE, dist_error=None, ThetaE_error=None, FL_FS=0):
	'''
	calculate the expected shift of combined center of light
	including luminous lens effects 
	expect dist and ThetaE in the same unit
	return shift in the same unit as dist and ThetaE
	if dist_error != 0 also the error is propagated 
	'''
	u = dist / ThetaE 
	# maximum shift at u = sqrt(2) / (1+FL_FS)
	mm = np.where(u < np.sqrt(2) / (1 + FL_FS))
	if hasattr(FL_FS, "__len__"):
		u[mm] = np.sqrt(2) / (1 + FL_FS[mm])
	else: 
		u[mm] = np.sqrt(2) / (1 + FL_FS)

	#shift in mas 
	shift_lum = u * ThetaE / (1 + FL_FS) * (1 + FL_FS * (u**2 + 3 \
		- u * np.sqrt(u**2 + 4)))/(u**2 + 2 + FL_FS * u * np.sqrt(u**2 + 4))


	if dist_error is None:
		return shift_lum
	elif len(dist_error.shape) == 2:
		u_e = dist_error / ThetaE_error 
		u_e = np.abs(u_e)
		FL_FS = FL_FS.reshape([-1,1]) * np.ones(dist_error.shape)
		mm = np.where(u_e < np.sqrt(2) / (1 + FL_FS))
		u_e[mm] = np.sqrt(2) / (1 + FL_FS[mm])

		shift_e = u_e * ThetaE_error / (1 + FL_FS) * (1 + FL_FS \
			* (u_e**2 + 3 - u_e * np.sqrt(u_e**2 + 4)))\
			/(u_e**2 + 2 + FL_FS * u_e * np.sqrt(u_e**2 + 4))
		shift_lum = np.percentile(shift_e, 50, axis = 1)
		shift_error_m = np.percentile(shift_e, 15.866, axis = 1) \
			- shift_lum
		shift_error_p = np.percentile(shift_e, 84.866, axis = 1) \
			- shift_lum
		return shift_lum, shift_error_p, shift_error_m
	else: 
		# calculate error
		u_error = u * np.sqrt(np.square(dist_error / dist) 
		+ np.square(ThetaE_error / ThetaE))
		# d shift_lum / du
		nn = (ThetaE * FL_FS * u * (-u**2 / np.sqrt(u**2 + 4) \
			- np.sqrt(u**2 + 4) + 2 * u))\
			/ ((FL_FS + 1) * (FL_FS * np.sqrt(u**2 + 4) * u + u**2 + 2)) \
			+ (ThetaE * (FL_FS * (u**2 - np.sqrt(u**2 + 4) * u + 3) + 1)) \
			/ ((FL_FS + 1) * (FL_FS * np.sqrt(u**2 + 4) * u + u**2 + 2)) \
			- (ThetaE * u * ((FL_FS * u**2) / np.sqrt(u**2 + 4) + FL_FS \
			* np.sqrt(u**2 + 4) + 2 * u) \
			* (FL_FS * (u**2 - np.sqrt(u**2 + 4) * u + 3) + 1))\
			/ (FL_FS + 1) / np.square(FL_FS * np.sqrt(u**2 + 4) * u + u**2 + 2)

		shift_lum_error = np.sqrt(nn**2 * u_error**2\
			+ shift_lum**2 * ThetaE_error**2 / ThetaE**2)
		return shift_lum, shift_lum_error

def calc_magnification(dist, ThetaE, dist_error, ThetaE_error,FL_FS):
	'''
	calculate the expected magnification including luminous lens effects 
	expect dist and ThetaE in the same unit
	return shift in the same unit as dist and ThetaE
	if dist_error != 0 also the error is propagated 
	'''
	u = dist / ThetaE

	A = (u **2 + 2) / (u * np.sqrt(u**2 + 4))
	magnification = 5 * np.log((FL_FS + A) / (FL_FS + 1)) / np.log(100)
	if dist_error is None: 
		return magnification
	elif len(dist_error.shape) == 2:
		FL_FS = FL_FS.reshape([-1,1]) * np.ones(dist_error.shape)
		u_e = dist_error / ThetaE_error 
		u_e = np.abs(u_e)
		A_e = (u_e**2 + 2) / (u_e * np.sqrt(u_e**2 + 4))
		mag_e = 5 * np.log((FL_FS + A_e)/(FL_FS + 1)) / np.log(100)
		magnification = np.percentile(mag_e, 50, axis = 1)
		mag_error_m = np.percentile(mag_e, 15.866, axis = 1) \
			- magnification
		mag_error_p = np.percentile(mag_e, 84.866, axis = 1) \
			- magnification
		return magnification, mag_error_p, mag_error_m
	else:
		# calculate error
		u_error = u * np.sqrt(np.square(dist_error / dist) 
		+ np.square(ThetaE_error / ThetaE))
		'''
		(u^2 + 2)/(u*sqrt(u^2+4))' 
		= ((u^2 + 2)' (u*sqrt(u^2+4)) - (u^2 + 2) (u*sqrt(u^2+4))') 
			/(u^2*(u^2+4))
		= (2*u^2 - (u^2+2)* (1+ u^2/u^2+4 )/(u^2*sqrt(...))
		= ((u^2-2)* (u^2+4) - (u^2+2)*u^2)/(u^2*sqrt...^3)
		= -8/(u^2*sqrt(u^2+4)^3)
		'''
		A_error = 8 / (u**2 * pow(u**2 + 4,3/2)) * u_error
		magnification_error = 5 / np.log(100) * A_error / (FL_FS + A)
		return magnification, magnification_error



def calculate_Effect(tab, approx = False, gaia = False, \
			error_percentile = False, nMC = 50000):
	"""
	calculate the microlensing effect for an given distance and Einsteinradius
	 both in mas
	It return the distance in values of the Einstein radius, the 
	shift in mas, the magnification in delta mag as well as their errors in 
	the same units.
	Formated as astropy.table.MaskedColumn object
	if approx: use approx dist and only determine shift_plus 
	if gaia: use L2 distance and Epoch, and include L2 in table discription 
	"""
	if approx or 'dist' not in tab.colnames:
		# only determine shift_plus 
		approx_shift_plus = calc_shift_plus(tab['approx_dist'], tab['ThetaE'])
		return MaskedColumn(approx_shift_plus, dtype = 'float64', \
			unit = 'mas', \
			description ='Approximated maximal astrometric shift of image (+)')
	else:
		ti = time.time()
		if gaia: 
			suffix = 'L2_'
			suffix_descr = ' for Lagrange Point L2' 
		else: 
			suffix = ''
			suffix_descr = '' 

		# get distance
		dist = tab[suffix + 'dist']
		dist_error = tab[suffix + 'dist_error']
		# if error_percentile: 
		# 	np.random.seed(550095)
		# 	dist_error2 = np.random.normal(dist, dist_error,\
		# 		(nMC,len(dist))).T
		# 
		ThetaE = tab['ThetaE']
		ThetaE_error = tab['ThetaE_error']
		# if error_percentile: 
		# 	np.random.seed(630201)
		# 	ThetaE_error2 = np.random.normal(ThetaE, ThetaE_error, \
		# 		(nMC,len(dist))).T
		# determine Flux Ratio
		FL_FS = pow(100,(tab['ob_phot_g_mean_mag'] - tab['phot_g_mean_mag'])/5)


		# determine normed impact parameter u
		u = dist / ThetaE
		if error_percentile: 

			u = np.zeros(len(tab))
			u_error = np.zeros(len(tab))
			u_error_p = np.zeros(len(tab))
			u_error_m = np.zeros(len(tab))
			shift = np.zeros(len(tab))
			shift_error = np.zeros(len(tab))
			shift_error_p = np.zeros(len(tab))
			shift_error_m = np.zeros(len(tab))
			shift_plus = np.zeros(len(tab))
			shift_plus_error = np.zeros(len(tab))
			shift_plus_error_p = np.zeros(len(tab))
			shift_plus_error_m = np.zeros(len(tab))
			shift_lum = np.zeros(len(tab))
			shift_lum_error = np.zeros(len(tab))
			shift_lum_error_p = np.zeros(len(tab))
			shift_lum_error_m = np.zeros(len(tab))
			magnification = np.zeros(len(tab))
			magnification_error = np.zeros(len(tab))
			magnification_error_p = np.zeros(len(tab))
			magnification_error_m = np.zeros(len(tab))
			np.random.seed(550095)
			v=np.arange(len(tab))
			for i in range(int(nMC/10000)):
				print("calculate_Effect %s %im:%is"%(i, (time.time()-ti) // 60,\
					(time.time()-ti) % 60))
				vv = (v<len(tab)/int(nMC/10000)*(i+1)) \
					& (v>=len(tab)/int(nMC/10000)*i) 


				dist_error2 = np.random.normal(dist[vv], dist_error[vv],\
				(nMC,int(np.sum(vv)))).T
				ThetaE_error2 = np.random.normal(ThetaE[vv], ThetaE_error[vv], \
					(nMC,int(np.sum(vv)))).T
				u_e = dist_error2 / ThetaE_error2
				u[vv] = np.percentile(u_e, 50, axis = 1)

				u_error_m[vv] = np.percentile(u_e, 15.866, axis = 1) \
					- u[vv]
				u_error_p[vv] = np.percentile(u_e, 84.866, axis = 1) \
					- u[vv]
				del u_e
				shift[vv], shift_error_p[vv],shift_error_m[vv] = calc_shift(
					dist[vv], ThetaE[vv], dist_error2, ThetaE_error2,FL_FS[vv])
				shift_plus[vv], shift_plus_error_p[vv],shift_plus_error_m[vv]= \
					calc_shift_plus(dist[vv], ThetaE[vv], dist_error2, \
					ThetaE_error2, FL_FS[vv])
				shift_lum[vv], shift_lum_error_p[vv], shift_lum_error_m[vv] = \
					calc_shift_lum(dist[vv], ThetaE[vv], dist_error2, \
					ThetaE_error2, FL_FS[vv])
				magnification[vv], magnification_error_p[vv],\
					magnification_error_m[vv] =  calc_magnification(
					dist[vv], ThetaE[vv], dist_error2, ThetaE_error2,FL_FS[vv])
			
			u_error = np.maximum(u_error_p, -u_error_m)
			shift_error = np.maximum(shift_error_p, -shift_error_m)
			shift_plus_error = np.maximum(shift_plus_error_p, \
					-shift_plus_error_m)
			shift_lum_error = np.maximum(shift_lum_error_p, -shift_lum_error_m)
			magnification_error = np.maximum(magnification_error_p, \
					-magnification_error_m)


			# transfer to Column objects in order to include in 
			# astropy.table.Table object
			u = MaskedColumn(u, dtype = 'float64', unit = '', description = \
				'Closest distance in einstein radii' + suffix_descr)
			u_error = MaskedColumn(u_error, dtype = 'float64', unit = '', \
				description ='Error of the closest distance in einstein radii'\
				+ suffix_descr)
			u_error_p = MaskedColumn(u_error_p, dtype = 'float64', unit = '', \
				description = \
				'Positiv error of the closest distance in einstein radii'\
				+ suffix_descr)
			u_error_m = MaskedColumn(u_error_m, dtype = 'float64', unit = '', \
				description = \
				'Negativ error of the closest distance in einstein radii'\
				+ suffix_descr)
			shift = MaskedColumn(shift, dtype = 'float64', unit = 'mas',\
				description =\
				'Expected shift of the center of light' \
				+ suffix_descr)
			shift_error = MaskedColumn(shift_error, dtype = 'float64', \
				unit = 'mas', description = \
				'Error of the shift of the center of light' + suffix_descr)
			shift_error_p = MaskedColumn(shift_error_p, dtype = 'float64', \
				unit = 'mas', description = \
				'Positiv error of the shift of the center of light' \
				+ suffix_descr)
			shift_error_m = MaskedColumn(shift_error_m, dtype = 'float64', \
				unit = 'mas', description = \
				'Negativ error of the shift of the center of light' \
				+ suffix_descr)
			
			shift_plus = MaskedColumn(shift_plus, dtype = 'float64',\
				description = 'Expected shift of the major image (+)'\
				+ suffix_descr,unit = 'mas')
			shift_plus_error = MaskedColumn(shift_plus_error, dtype='float64',\
				unit = 'mas', description = \
				'Error of the shift of the major image (+)' + suffix_descr)
			shift_plus_error_p = MaskedColumn(shift_plus_error_p, 
				dtype='float64', unit = 'mas', description = \
				'Positiv error of the shift of image (+)' + suffix_descr)
			shift_plus_error_m = MaskedColumn(shift_plus_error_m,\
				dtype='float64', unit = 'mas', description = \
				'Negativ error of the shift of the major image (+)' \
				+ suffix_descr)
			
			shift_lum = MaskedColumn(shift_lum, dtype='float64', unit = 'mas',\
				description = 
				'Expected maximal shift includig luminous-lens effect' \
				+ suffix_descr)
			shift_lum_error = MaskedColumn(shift_lum_error, dtype = 'float64',\
				unit = 'mas', description = \
				'Error of the shift includig lens-luminosity effect' \
				+ suffix_descr)
			shift_lum_error_p = MaskedColumn(shift_lum_error_p, \
				dtype = 'float64',unit = 'mas', description = \
				'Positiv error of the shift includig luminous-lens effect'\
				+ suffix_descr)
			shift_lum_error_m = MaskedColumn(shift_lum_error_m,\
				dtype = 'float64', unit = 'mas', description = \
				'Negativ error of the shift includig luminous-lens effect' \
				+ suffix_descr)

			magnification = MaskedColumn(magnification, dtype = 'float64',\
				unit = 'mag', description = \
				'Magnification in magnitudes' + suffix_descr)
			magnification_error = MaskedColumn(magnification_error, 
				dtype = 'float64', unit = 'mag', description = \
				'Error of the magnification in magnitudes' + suffix_descr)
			magnification_error_p = MaskedColumn(magnification_error_p, 
				dtype = 'float64', unit = 'mag', description = \
				'Positiv error of the magnification in magnitudes' + \
				suffix_descr)
			magnification_error_m = MaskedColumn(magnification_error_m, 
				dtype = 'float64', unit = 'mag', description = \
				'Negativ error of the magnification in magnitudes' + \
				suffix_descr)
			return u, u_error,u_error_p,u_error_m, shift, shift_error, \
				shift_error_p, shift_error_m, shift_plus, shift_plus_error, \
				shift_plus_error_p, shift_plus_error_m, shift_lum, \
				shift_lum_error, shift_lum_error_p, shift_lum_error_m,\
				magnification, magnification_error, magnification_error_p,\
				magnification_error_m

		else: 

			u_error = u * np.sqrt((dist_error/dist)**2 \
				+ (ThetaE_error/ThetaE)**2)

			# calculate effects
			shift, shift_error = calc_shift(dist, ThetaE, dist_error,\
				ThetaE_error,FL_FS)
			shift_plus, shift_plus_error = calc_shift_plus(dist, ThetaE,\
				dist_error, ThetaE_error,FL_FS)
			shift_lum, shift_lum_error = calc_shift_lum(dist, ThetaE, \
				dist_error, ThetaE_error,FL_FS)
			magnification, magnification_error = calc_magnification(dist, \
				ThetaE, dist_error, ThetaE_error,FL_FS)

			# transfer to Column objects in order to include in 
			# astropy.table.Table object
			u = MaskedColumn(u, dtype = 'float64', unit = '', description = \
				'Closest distance in einstein radii' + suffix_descr)
			u_error = MaskedColumn(u_error, dtype = 'float64', unit = '', \
				description ='Error of the closest distance in einstein radii'\
				+ suffix_descr)
			shift = MaskedColumn(shift, dtype = 'float64', unit = 'mas',\
				description =\
				'Expected shift of the center of light' \
				+ suffix_descr)
			shift_error = MaskedColumn(shift_error, dtype = 'float64', \
				unit = 'mas', description = \
				'Error of the shift of the center of light' + suffix_descr)
			shift_plus = MaskedColumn(shift_plus, dtype = 'float64',\
				description = 'Expected shift of image (+)' \
				+ suffix_descr,unit = 'mas')
			shift_plus_error = MaskedColumn(shift_plus_error, dtype='float64',\
				unit = 'mas', description = 'Error of the shift of image (+)' \
				+ suffix_descr)
			shift_lum = MaskedColumn(shift_lum, dtype='float64', unit = 'mas',\
				description = 
				'Expected shift includig luminous-lens effect' \
				+ suffix_descr)
			shift_lum_error = MaskedColumn(shift_lum_error, dtype = 'float64',\
				unit = 'mas', description = \
				'Error of the shift includig luminous-lens effect' \
				+ suffix_descr)
			magnification = MaskedColumn(magnification, dtype = 'float64',\
				unit = 'mag', description = \
				'Magnification in magnitudes' + suffix_descr)
			magnification_error = MaskedColumn(magnification_error, 
				dtype = 'float64', unit = 'mag', description = \
				'Error of the magnification in magnitudes' + suffix_descr)

			return u, u_error, shift, shift_error, shift_plus, \
				shift_plus_error, shift_lum, shift_lum_error, magnification, \
				magnification_error



