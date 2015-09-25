function calcdamping, T, Pe, nH, lambda, doppler_width, alpha=alpha, sigma0=sigma0, mass=mass

	PME = 9.10938188d-28
	PH = 6.62606876d-27
	PC = 2.99792458d10
	PUMA = 1.66053873d-24
	PK = 1.3806503d-16
	PE = 4.8032d-10
		
; Radiative damping
	gamma_nat = 0.22233d0 / (lambda*1.d-8)^2
	
; Stark damping
	gamma_stark = 19.4 + 2.d0/3.d0*(-14.d0) + alog10(Pe) - 5.d0/6.d0*alog10(T)
	gamma_stark = 10.d0^gamma_stark

	gamma6 = nH * 0.d0
	
; Hydrogen damping	
	if (keyword_set(sigma0) and keyword_set(alpha) and keyword_set(mass)) then begin
		a0 = (PH/(2.d0*!DPI))^2 / (PME*PE^2)
		sigma = sigma0 * a0^2
		v0 = 1.d6
		m1 = mass
		m2 = 1.008d0
		reduced_mass = (m1*m2) / (m1+m2) * PUMA
		t1 = (alpha/2.d0)*alog(4.d0/!DPI)
		t2 = lngamma((4.d0-alpha)/2.d0)
		t3 = alpha*alog(v0)
		t4 = alog(sigma)
		
		gamma6 = t1 + t2 + t3 + t4 + 0.5d0*(1.d0-alpha)*alog((8.d0*PK*T)/(!DPI*reduced_mass))
		gamma6 = exp(gamma6)
				
		gamma6 = 2.d0 * gamma6 * nH
		
	endif
	
	cte = 1.d-8 / (4.d0 * !DPI * PC)
	return, cte * (gamma6+gamma_nat+gamma_stark) * lambda^2 / (doppler_width * lambda / PC)
end
