; This code gets the intensity through a plane-parallel atmosphere integrating the
; transfer equation using a short-characteristics scheme
;--------------------------------------------------------------
;--------------------------------------------------------------
; MONOCHROMATIC FORMAL SOLUTION OF THE TRANSFER EQUATION USING SHORT-CHARACTERISTICS
; GIVEN THE OPACITY AND EMMISIVITY 
;--------------------------------------------------------------
;--------------------------------------------------------------
; Intensity units depend on the parameter's units
;--------------------------------------------------------------
; All variables have index 0 in the bottom of the atmosphere
; and index N-1 in the top
;--------------------------------------------------------------
; x : spatial points in the atmosphere
; chi : opacity
; eta : emmisivity
; mu : angle at which one wants to observe
; I0 : boundary condition at the bottom
; Inten(0..size(x)) : vector of intensity through the atmosphere
;--------------------------------------------------------------

function shortcar, x, chi, eta, mu, I0
	
	tamano = n_elements(x)
	Inten = dblarr(tamano)
	Inten(0) = I0

; Points where parabolic SC can be used (interior points)
	for k=1, tamano-2 do begin

		etatot = dblarr(3)
		chitot = dblarr(3)
		Stot = dblarr(3)
		
		km = k-1
		kp = k+1
		
		dm = abs((abs(x(km))-abs(x(k))) / mu)
		dp = abs((abs(x(k))-abs(x(kp))) / mu)
		
		chitot(0) = chi(km)
		chitot(1) = chi(k)
		chitot(2) = chi(kp)
		
		etatot(0) = eta(km)
		etatot(1) = eta(k)
		etatot(2) = eta(kp)
		
		if (chitot(0) ne 0.d0) then Stot(0) = etatot(0) / chitot(0)
		if (chitot(1) ne 0.d0) then Stot(1) = etatot(1) / chitot(1)
		if (chitot(2) ne 0.d0) then Stot(2) = etatot(2) / chitot(2)
				
		dtm = 0.5d0 * (chitot(0) + chitot(1)) * dm
		dtp = 0.5d0 * (chitot(1) + chitot(2)) * dp
		
		if (dtm > 0.001) then exu = 1.d0*exp(-dtm) else exu = 1.d0 - dtm + dtm^2/2.d0
		
		Inten(k) = Inten(k-1) * exu + total(parabcoef(dtm,dtp) * Stot)
		
	endfor
	
; Points where linear SC have to be used (last point)

	Stot = dblarr(2)
	chitot = dblarr(2)
	etatot = dblarr(2)
	
	dm = abs((abs(x(tamano-2))-abs(x(tamano-1))) / mu)
	
	chitot(0) = chi(tamano-2)
	chitot(1) = chi(tamano-1)
	
	etatot(0) = eta(tamano-2)
	etatot(1) = eta(tamano-1)
	
	if (chitot(0) ne 0) then Stot(0) = etatot(0) / chitot(0)
	if (chitot(1) ne 0) then Stot(1) = etatot(1) / chitot(1)
	
	dtm = 0.5d0 * (chitot(0) + chitot(1)) * dm
	
	if (dtm > 0.001) then exu = exp(-dtm) else exu = 1.d0 - dtm + dtm^2/2.d0
	
	Inten(tamano-1) = Inten(tamano-2) * exu + total(linearcoef(dtm) * Stot(0:1))
	
	return, Inten[tamano-1]
end