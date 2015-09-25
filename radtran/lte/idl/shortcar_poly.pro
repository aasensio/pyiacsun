;--------------------------------------------------------------
;--------------------------------------------------------------
; POLYCHROMATIC FORMAL SOLUTION OF THE TRANSFER EQUATION USING SHORT-CHARACTERISTICS
; GIVEN THE OPACITY AND EMMISIVITY 
;--------------------------------------------------------------
;--------------------------------------------------------------
; Intensity units depend on the parameter's units
;--------------------------------------------------------------
; All variables have index 0 in the bottom of the atmosphere
; and index N-1 in the top
;--------------------------------------------------------------
; x : spatial points in the atmosphere
; chi : opacity in format chi(freq,x)
; eta : emmisivity in format eta(freq,x)
; mu : angle at which one wants to observe
; I0 : boundary condition at the bottom in format I0(freq)
; Inten(*,0..size(x)) : vector of intensity through the atmosphere for all the frequencies
;--------------------------------------------------------------

function shortcar_poly, x, chi, eta, mu, I0
	
	tamano = n_elements(x)
	n_freq = n_elements(I0)
	Inten = dblarr(n_freq,tamano)
	Inten(*,0) = I0
	dtm = dblarr(n_freq)
	dtp = dblarr(n_freq)

; Points where parabolic SC can be used (interior points)
	for k=1, tamano-2 do begin

		etatot = dblarr(n_freq,3)
		chitot = dblarr(n_freq,3)
		Stot = dblarr(n_freq,3)
		exu = dblarr(n_freq)
		
		km = k-1
		kp = k+1
		
		dm = abs((abs(x(km))-abs(x(k))) / mu)
		dp = abs((abs(x(k))-abs(x(kp))) / mu)

		
		chitot(*,0) = chi(*,km)
		chitot(*,1) = chi(*,k)
		chitot(*,2) = chi(*,kp)
		
		etatot(*,0) = eta(*,km)
		etatot(*,1) = eta(*,k)
		etatot(*,2) = eta(*,kp)
		
		inds = where(chitot(*,0) ne 0.d0)
		if (inds(0) ne -1) then Stot(inds,0) = etatot(inds,0) / chitot(inds,0)
		inds = where(chitot(*,1) ne 0.d0)
		if (inds(0) ne -1) then Stot(inds,1) = etatot(inds,1) / chitot(inds,1)
		inds = where(chitot(*,2) ne 0.d0)
		if (inds(0) ne -1) then Stot(inds,2) = etatot(inds,2) / chitot(inds,2)
						
		dtm = 0.5d0 * (chitot(*,0) + chitot(*,1)) * dm
		dtp = 0.5d0 * (chitot(*,1) + chitot(*,2)) * dp
		
		inds = where(dtm gt 0.001)
		if (inds(0) ne -1) then exu(inds) = 1.d0*exp(-dtm(inds))
		inds = where(dtm lt 0.001)
		if (inds(0) ne -1) then exu(inds) = 1.d0 - dtm(inds) + dtm(inds)^2/2.d0

		Inten(*,k) = Inten(*,k-1) * exu + total(parab_polycoef(dtm,dtp,n_freq) * Stot,2)
	endfor
	
; Points where linear SC must be used (last point)

	Stot = dblarr(n_freq,2)
	chitot = dblarr(n_freq,2)
	etatot = dblarr(n_freq,2)
	exu = dblarr(n_freq)
	
	dm = abs((abs(x(tamano-2))-abs(x(tamano-1))) / mu)
	
	chitot(*,0) = chi(*,tamano-2)
	chitot(*,1) = chi(*,tamano-1)
	
	etatot(*,0) = eta(*,tamano-2)
	etatot(*,1) = eta(*,tamano-1)
	
	inds = where(chitot(*,0) ne 0.d0)
	if (inds(0) ne -1) then Stot(inds,0) = etatot(inds,0) / chitot(inds,0)
	inds = where(chitot(*,1) ne 0.d0)
	if (inds(0) ne -1) then Stot(inds,1) = etatot(inds,1) / chitot(inds,1)
	
	dtm = 0.5d0 * (chitot(*,0) + chitot(*,1)) * dm
	
	inds = where(dtm gt 0.001)
	if (inds(0) ne -1) then exu(inds) = 1.d0*exp(-dtm(inds))
	inds = where(dtm lt 0.001)
	if (inds(0) ne -1) then exu(inds) = 1.d0 - dtm(inds) + dtm(inds)^2/2.d0
	
	Inten(*,tamano-1) = Inten(*,tamano-2) * exu + total(linear_polycoef(dtm,n_freq) * Stot,2)
	
	return, Inten[*,tamano-1]
end