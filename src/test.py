import numpy as np

from astropy.coordinates import get_sun
from astropy.time import Time
Gaia_epoch = 2016
timeGaia = 0.
T = Gaia_epoch+timeGaia

ra = 250
for dec in [23,-12, 60, -45.2342]:
	parallax = 1500
	rad2pi = lambda x: x+(x-abs(x))/x*np.pi
	spherToCart = lambda x,y: np.array([np.cos(x/180.*np.pi)*np.cos(y/180.*np.pi),np.sin(x/180.*np.pi)*np.cos(y/180.*np.pi),np.sin(y/180.*np.pi)])
	Carttospher = lambda x,y,z: np.array([rad2pi(np.arctan2(y,x)), np.arcsin(z/np.sqrt(x**2+y**2+z**2))])*180./np.pi

	sd, cd = np.sin(dec/180.*np.pi), np.cos(dec/180.*np.pi)
	sa, ca = np.sin(ra/180.*np.pi), np.cos(ra/180.*np.pi)
	baryVec =spherToCart(ra,dec)
	earthcoord = get_sun(Time(timeGaia + Gaia_epoch, format = 'jyear'))

	earthVec = spherToCart(earthcoord.ra.deg, earthcoord.dec.deg)\
				* earthcoord.distance.pc
	geoVec = baryVec + earthVec * parallax/1e3

	print(Carttospher(*baryVec))
	print(Carttospher(*geoVec))
	print((Carttospher(*geoVec)-Carttospher(*baryVec))*3.6e6)
	# print((baryVec))
	# print((geoVec))
	# print(((geoVec)-(baryVec))*1e8)



	scsc = np.array([np.sin(ra/180.*np.pi), np.cos(ra/180.*np.pi), \
					np.sin(dec/180.*np.pi), np.cos(dec/180.*np.pi)])
	gs = -1.0 * get_sun(Time(T, format = 'jyear')).cartesian.xyz.to('AU').value
	earth = gs.T
	loc_vec = np.array([(scsc[0] * earth[0] - scsc[1] *  earth[1])					/ np.maximum(scsc[3],1e-6),
		scsc[1] * scsc[2] *  earth[0] + scsc[0] * scsc[2]*  earth[1] - scsc[3] *  earth[2]])

	radec_vec = np.array([ra,dec])
	px_vec = parallax
	radec_cur  = radec_vec + px_vec/3.6e6 * loc_vec
	print(radec_vec)
	print(radec_cur)
	print((radec_cur-radec_vec)*3.6e6)
	# print(spherToCart(*radec_vec))
	# print(spherToCart(*radec_cur))
	# print((spherToCart(*radec_cur)-spherToCart(*radec_vec))*1e8)
	print('---------------')
