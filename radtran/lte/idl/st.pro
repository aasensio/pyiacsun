;--------------------------

	function st,pe,ph

; PURPOSE:
; 	it computes the absorption coefficient per h atom (cm**2) 
;	produced by free electrons (thompson scattering).
;
; INPUT:
;	pe	electron pressure (dny/cm**2)
;	ph	partial pressure of neutral h atom (dy/cm**2)
;
; OUTPUT:
;	st	thompson absorption coefficient
;
; MODIFICATION HISTORY:
; started:	23/2/90
; finished:	23/2/90		fortran (for sos)
;		Sep. 1st, 95	fortran > idl (jorge)

; definitions
	const=.6652e-24
;
	st=const*(pe/ph)
	return,st
	end
