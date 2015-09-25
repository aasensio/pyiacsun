;-----------------------
	
	function kh_bf,t,pe,rlambda

; PURPOSE:
; 	It computes negative hydrogen bound-free continuum absorption coeff.
; 	per h atom (cm**2), according to the following reference;
;
; 	john: 1989, a&a 193, 189 (he has a couple of mistakes)
;		
;
; INPUT:
;	t	temperature (k)
;	pe	electron pressure
;	rlambda	wavelenght (a) (.gt. 1250a)
;
; output:
;	kh_bf	continuum apsorption coefficient
;
; MODIFICATION HISTORY:
;	(started:	?
; 	finised:	20/2/90)	for SOS
;			Sept 1st, 95	fortran > idl (jorge)

; definitions & constants
	lambda_0=1.6419d0
	alfa=1.439d4
	const=.75d-18
	cc=[152.519d0,49.534d0,-118.858d0,92.536d0,-34.194d0,4.982d0]
; h_bf
	lambda=rlambda/1.e4
	if(lambda lt lambda_0)then begin
	   com=1./lambda-1./lambda_0
	   e=findgen(6)/2.	;idl stile
	   sigma=total(cc*com^e);idl stile
	   sigma=lambda^3.*com^1.5*sigma
	   kh_bf=t^(-2.5)*exp((alfa/lambda_0)/t)*(1.-exp(-(alfa/lambda)/t))
	   kh_bf=(const*sigma)*kh_bf
	endif else begin
	   kh_bf=0.d0*t
	endelse
; return
	kh_bf=kh_bf*pe
	return,float(kh_bf)
	end
