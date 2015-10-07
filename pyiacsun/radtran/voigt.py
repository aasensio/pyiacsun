__all__ = ['voigt']
from scipy.special import wofz
###########################################
# Returns the Voigt function for an axis of wavelengths l and damping parameter a
###########################################
def voigt(l, a):
	z = l + 1j*a
	I = wofz(z).real
	return I
