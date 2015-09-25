;-----------------------


	function srh,rlambda


; PURPOSE:
;	Calculation of the  rayleigh scattering on hydrogen atoms 
;	per h atom (cm**2) as given by
; 
;		landi degl'innocenti: 1976, astron. astro. supl. ser. 25, 379
;		( he quote dalgarno: 1962, geophys. corp. ann. techn. rep. no.
;		  62-28-a; the polynomial has an error)
; INPUT:
;	rlambda	wavelenght (a)
;
; output:
;	srh	continuum apsorption coefficient
;
; MODIFICATION HISTORY:
; (started:	20/2/90
; finished:	21/2/90)
;	Sept 1st, 95	fortran (sos) > idl. Jorge


; definitions
;
	c=[5.799d-13,1.422d-6,2.784d0]
;
	srh=(c(0)+(c(1)+c(2)/rlambda^2.)/rlambda^2.)/rlambda^4.
	return,float(srh)
	end
