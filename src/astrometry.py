import numpy as np
from astropy.coordinates import get_sun
from astropy.time import Time

#sub packages
from setup import Gaia_epoch,zeropoint
from utils import DEG		

def spherToCart(theta, phi):
	"""returns a 3-cartesian unit vector pointing to longitude theta,
	latitude phi.

	The angles are in rad.
	"""
	cp = np.cos(phi)
	return np.cos(theta) * cp, np.sin(theta) * cp, np.sin(phi)

def dirVecToCelCoos(dirVec):
	"""returns alpha, delta in degrees for the direction vector dirVec.

	>>> dirVecToCelCoos(computeUnitSphereCoords(25.25, 12.125))
	(25.25, 12.125)
	>>> dirVecToCelCoos(computeUnitSphereCoords(25.25, 12.125)*16)
	(25.25, 12.125)
	>>> "%g,%g"%dirVecToCelCoos(computeUnitSphereCoords(25.25, 12.125)+
	... computeUnitSphereCoords(30.75, 20.0))
	'27.9455,16.0801'
	"""

	dirVec = dirVec.normalized()
	alpha = np.arctan2(dirVec.y, dirVec.x)
	alpha[alpha < 0] = alpha[alpha < 0] + 2 * np.pi
	return alpha/DEG, np.arcsin(dirVec.z)/DEG



class Vector3(object):
	"""is a 3d vector that responds to both .x... and [0]...
	>>> x, y = Vector3(1,2,3), Vector3(2,3,4)
	>>> x+y
	Vector3(3.000000,5.000000,7.000000)
	>>> 4*x
	Vector3(4.000000,8.000000,12.000000)
	>>> x*4
	Vector3(4.000000,8.000000,12.000000)
	>>> x*y
	20
	>>> "%.6f"%abs(x)
	'3.741657'
	>>> print abs((x+y).normalized())
	1.0
	"""
	def __init__(self, x, y=None, z=None):
		if isinstance(x, tuple):
			self.coos = np.array([*x]).reshape(3,-1)
		else:
			self.coos = np.array([x, y, z]).reshape(3,-1)
	'''
	def __repr__(self):
		return "Vector3(%f,%f,%f)"%tuple(self.coos)
	
	def __str__(self):
		def cutoff(c):
			if abs(c)<1e-10:
				return 0
			else:
				return c
		rounded = [cutoff(c) for c in self.coos]
		return "[%.2g,%.2g,%.2g]"%tuple(rounded)
	'''
	def __getitem__(self, index):
		return Vector3(*self.coos[:,index])
	def __mul__(self, other):
		"""does either scalar multiplication if other is not a Vector3, or
		a scalar product.
		"""
		if isinstance(other, Vector3):
			return self.x*other.x + self.y*other.y + self.z*other.z
		else:
			return Vector3(self.x*other, self.y*other, self.z*other)
	__rmul__ = __mul__
	def __truediv__(self, scalar):
		return Vector3(self.x/scalar, self.y/scalar, self.z/scalar)
	def __add__(self, other):
		return Vector3(self.x+other.x, self.y+other.y, self.z+other.z)
	def __sub__(self, other):
		return Vector3(self.x-other.x, self.y-other.y, self.z-other.z)
	def __abs__(self):
		return np.sqrt(self.x**2 + self.y**2 + self.z**2)
	def cross(self, other):
		return Vector3(self.y*other.z - self.z*other.y,
			self.z*other.x - self.x*other.z,
			self.x*other.y - self.y*other.x)
	def normalized(self):
		return self/abs(self)
	def getx(self): return self.coos[0]
	def setx(self, x): self.coos[0] = x
	x = property(getx, setx)
	def gety(self): return self.coos[1]
	def sety(self, y): self.coos[1] = y
	y = property(gety, sety)
	def getz(self): return self.coos[2]
	def setz(self, z): self.coos[2] = z
	z = property(getz, setz)

def getGCDist_3Vec(pos1, pos2):
	"""returns the distance along a great circle between two points. 
	The distance is in degrees, the input positions are Vektor3.
	"""
	scalarprod = (pos1.normalized()) * (pos2.normalized())
	if scalarprod > 0.9:
		dif = pos1.normalized() - pos2.normalized()
		return 2e0 * np.arcsin(abs(dif) / 2e0) / DEG
	return np.arccos(scalarprod)/DEG

def movePm_3Vec(raDeg, decDeg, pmra, pmdec, timeGaia, foreshort = 0):
	"""
	returns cartesian coordinates for an object with pos ra, dec 
	and pm pmra after Gaia epoch.
	pmra has to have cos(dec) applied, position in deg prop.motion in mas 
	the time unit is yours to choose.
	"""
	ra, dec = raDeg*DEG, decDeg*DEG
	sd, cd = np.sin(dec), np.cos(dec)
	sa, ca = np.sin(ra), np.cos(ra)
	pmra, pmdec = pmra * DEG / 3.6e6 , pmdec * DEG / 3.6e6
	muAbs = np.sqrt(pmra**2 + pmdec**2);
	muTot = muAbs + 0.5e0 * foreshort * timeGaia;
	if muAbs < 1e-20: # no proper motion 
		dirVec = Vector3(cd*ca, cd*sa, sd)
		return dirVec
	# this is according to  
	dirA = pmra / muAbs;
	dirD = pmdec / muAbs;
	sinMot = np.sin(muTot * timeGaia);
	cosMot = np.cos(muTot * timeGaia);
	dirVec = Vector3(
		-sd * ca * dirD * sinMot - sa * dirA * sinMot + cd * ca * cosMot,
		-sd * sa * dirD * sinMot + ca * dirA * sinMot + cd * sa * cosMot,
		+ cd * dirD * sinMot + sd * cosMot)
	return dirVec


def parallax_correction(baryVec, earthVec, parallax):
	""" returns	parrallax corrected cartesian Coordinates as Vector3 
		baryVec and earthVec as Vector3, parallax in mas	
	"""
	geoVec = baryVec + earthVec * parallax/1e3
	return geoVec


def movePm_parallax(raDeg, decDeg, pmra, pmdec, parallax, timeGaia, \
		earthVec = None, foreshort = 0, gaia = False):
	"""
	returns cartesian coordinates for an object with pos ra, dec 
	and pm pmra after Gaia epoch.
	pmra has to have cos(dec) applied, position in deg, prop.motion and 
	parallax in mas , the time in units of years after 2015.5
	"""
	baryVec = movePm_3Vec(raDeg, decDeg, pmra, pmdec, timeGaia, \
		foreshort = foreshort)
	if parallax >= zeropoint:
		if earthVec is None:
			earthVec = pos_sun(timeGaia)
		if gaia == True:
			geoVec = parallax_correction(baryVec, earthVec, parallax * 1.01)
		else:
			geoVec = parallax_correction(baryVec, earthVec, parallax)
		return geoVec
	else: 
		return baryVec



def pos_sun(timeGaia):
	"""
	retruns earth positions for t in julian years (365.25 days) after
	Gaia_reference epoch
	"""
	earthcoord = get_sun(Time(timeGaia + Gaia_epoch, format = 'jyear'))
	return Vector3(spherToCart(earthcoord.ra.rad, earthcoord.dec.rad)) \
			* earthcoord.distance.pc

epoch_predef = np.arange(-15, 70, 1/52)
pos_sun_predef = pos_sun(epoch_predef)

