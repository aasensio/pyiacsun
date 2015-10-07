# Author: cdiazbas@iac.es
# Code: Translation of IDL file
#       @vivivum: https://github.com/vivivum/MilosIDL/
#		Author: D. Orozco Suarez & J.C. del Toro Iniesta


def fvoigt(damp,vv):

	"""
	Extract spectral data from the Kitt Peak FTS-Spectral-Atlas
	as provided by H. Neckel, Hamburg.

	INPUTS: 
			DAMP: A scalar with the damping parameter
			VV: Wavelength axis usually in Doppler units.

	OUTPUTS: 
			H: Voigt function
			F: Faraday-Voigt function

	NOTES: 
			A rational approximation to the complex error function is used
			after Hui, Armstrong, and Wray(1978, JQSRT 19, 509). H and F are 
			the real and imaginary parts of such function, respectively.
			The procedure is inspired on that in SIR (Ruiz Cobo & del Toro 
			Iniesta 1992, ApJ 398, 385). On its turn, that routine was taken
			from modifications by A. Wittmann (1986) to modifications by S.K.
			Solanki (1985) to an original FORTRAN routine written by J.W. Harvey
			and A. Nordlund.
	"""

	import numpy as np

	A = [122.607931777104326, 214.382388694706425, 181.928533092181549,\
		93.155580458138441, 30.180142196210589, 5.912626209773153,\
		0.564189583562615]

	B = [122.60793177387535, 352.730625110963558, 457.334478783897737,\
		348.703917719495792, 170.354001821091472, 53.992906912940207,\
		10.479857114260399,1.]

	z = np.array(damp*np.ones(len(vv)) + -abs(vv)*1j)

	Z = ((((((A[6]*z+A[5])*z+A[4])*z+A[3])*z+A[2])*z+A[1])*z+A[0])/\
	(((((((z+B[6])*z+B[5])*z+B[4])*z+B[3])*z+B[2])*z+B[1])*z+B[0])

	h = np.real(Z)
	f = np.sign(vv)*np.imag(Z)*0.5

	return [h,f]



if __name__ == "__main__":

	import numpy as np
	import matplotlib.pyplot as plt

	uvals = np.linspace(-20., 20., 200)
	a = 2.E-1

	[h,f] = fvoigt(a,uvals)

	plt.plot(uvals,f,'k-')
	plt.plot(uvals,h,'r-')
	plt.show()



