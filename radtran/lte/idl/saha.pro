;-----------------

	function SAHA,t,pe,ul,uu,chi	

; PURPOSE:
;	 it computes the ionization equilibrium between two ionizations states
;	 according to saha's equation.
;
; INPUT:
;	t 	temperature (k)
;	pe	electron pressure (dny/cm2)
;	ul	partition function for the ionization l
;	uu	partition function for the ionization u=l+1
;	chi	ionization energy  from l to u (ev)
;	
; OUTPUT:
;	saha	number of atoms with ionization u=l+1/number with l
; 
; MODIFICATION HISTORY:
;		Comming from sos (fortran)
;	August 31, 95	IDL (jorge)

;
	const=.6665*uu/ul
	c1=-5040.*chi
; let's go
	saha=const*t^2.5*10.^(c1/t)/pe
;
	return,saha
	end
