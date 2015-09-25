;--------------------------


	function srh2,rlambda,ph,ph2


;PURPOSE:
;	
; 	To compute the rayleigh scattering on h2 molecules 
;	per h atom (cm**2) as given by
;
; 	landi degl'innocenti: 1976, astron. astro. supl. ser. 25, 379
;		( he quote dalgarno & williams: 1962, astrophys. j. 136, 960)
;
; INPUT:
;	rlambda	wavelenght (a)
;	ph	hydrogen partial pressure (dy/cm**2)
;	ph2	h2         "        "       "
;	
;OUTPUT
;	srh2	continuum apsorption coefficient
;
; MODIFICATION HISTORY:
; (started:	23/2/90		fortran (sos)
; finished:	23/2/90)
;	Sept 1st, 95		fortran > idl (jorge)


; definitions & parameters
	cc=[8.14d-13,1.28d-6,1.61]
; let's go
	srh2=(cc(0)+(cc(1)+cc(2)/rlambda^2.)/rlambda^2.)/rlambda^4.
	srh2=srh2*(ph2/ph)
;
	return,float(srh2)
	end
