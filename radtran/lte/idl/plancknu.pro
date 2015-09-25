function plancknu, T, nu, mks=mks, rayleigh=rayleigh
; T [K]
; nu [Hz]
; Returns the Planck's function in erg cm^-2 s^-1 Hz^-1 sr^-1
	if (not keyword_set(mks)) then begin
		hplanck = 6.626d-27   ;(erg*s)
		clight = 2.99793d10   ;(cm/s)
		kboltz = 1.381d-16  ;(erg/K)		
	endif else begin
		hplanck = 6.626d-34   ;(J*s)
		clight = 2.99793d8   ;(m/s)
		kboltz = 1.381d-23  ;(J/K)
	endelse

	if (keyword_set(rayleigh)) then begin
		return, 2.d0*kboltz*nu/clight^2*T
	endif else begin	
		factor1 = (2.d0 * hplanck / clight^2) * nu^3

		h_kt = hplanck / (kboltz * T)

		factor2 = 1.d0 / (exp(h_kt * nu) - 1.d0)
	
		return, factor1 * factor2	
	endelse
	
end
