@hydros
@gaspress
@kappa_c
@computeHeight
@shortcar_poly
@shortcar

; Solve the RT equation for the given model and spectral line

function synth_lte, mu, atmosphere, line, spectral, cont
common ATMDAT,w,ab,ei1,ei2,sym

	PK = 1.3806503d-16
	UMA = 1.66053873d-24
	PE0 = 4.8032d-10
	PC = 2.99792458d10
	PME = 9.10938188d-28
	PH = 6.62606876d-27
	OPA = !DPI * PE0^2 / (PME * PC)
	alpha_cte = PE0 / (4.d0*!DPI*PME*PC^2)
		
; Find the order of the logtau(500) axis
	if (atmosphere.ltau500[0] gt atmosphere.ltau500[0]) then begin
		order = 0
	endif else begin
		order = 1
	endelse
	
	atmdat,1
	
	T = atmosphere.T
	vmic = atmosphere.vmic
	ltau500 = atmosphere.ltau500
	vmac = atmosphere.vmac
	
; Put the atmosphere in hydrostatic equilibrium
	hydros, ltau500, T, Pe, pg	
	
; Compute the partial pressures of H, H-, H+, H2, H2+
	gaspress, Pe, T, p_h,phmi,phpl,ph2,ph2pl,pg
	nH = p_H / (PK*T)
	
	n = n_elements(T)
	
	opacity500 = kappa_c(Pe,T,5000.d0,htoverv)
	opacity500 = opacity500 * htoverv
	tau500 = 10.d0^ltau500
		
; Compute height from tau500
	z = computeHeight(tau500, opacity500)	
		
; Spectral line of interest
	nlines = n_elements(line.element)
	
	temp = {lambda: ptr_new(), specI: ptr_new(), contI: 0.d0}
	
	spec = replicate(temp, nlines)
	
	for j = 0, nlines-1 do begin
	
		element = line.element[j]
		ioniz = line.ioniz[j]
		lambda0 = line.lambda0[j]
		nu0 = line.nu0[j]
		Elow = line.Elow[j]
		gf = line.gf[j]
		alpha_ABO = line.alpha_ABO[j]
		sigma_ABO = line.sigma_ABO[j]
		lambda_from = line.lambda_from[j]
		lambda_to = line.lambda_to[j]
		lambda_step = line.lambda_step[j]
		
		mass = w[element-1]
				
		doppler_width = sqrt(2.d0*PK*T/(mass*UMA)+vmic^2)		
			
; Carry out the synthesis
		nlambda = (lambda_to-lambda_from)/lambda_step + 1
		lambda = findgen(nlambda) * lambda_step + lambda_from
		nu = PC / (lambda*1.d-8)
		
		eta = fltarr(nlambda,n)
		etaC = fltarr(nlambda,n)
		epsilon = fltarr(nlambda,n)
		epsilonC = fltarr(nlambda,n)
		
		damping = fltarr(n)
		opacityL = fltarr(n)
		
; Continuum opacity
		opacityC = kappa_c(Pe,T,lambda0,htoverv)
		opacityC = opacityC * htoverv		
			
; Calculate line opacity at each depth
		for i = 0, n-1 do begin
			damping[i] = calcdamping(T[i], Pe[i], nH[i], lambda0, $
				doppler_width[i], alpha=alpha_ABO, sigma0=sigma_ABO, mass=mass)
												
			partition, T[i], element, u1, u2, u3
			n1overn0 = saha(T[i], Pe[i], u1, u2, ei1[element-1])
			n2overn1 = saha(T[i], Pe[i], u2, u3, ei2[element-1])
			
			case ioniz of
				1 : begin
						niovern = 1.d0 / (1.d0 + n1overn0 + n2overn1 * n1overn0)
						ui = u1
					end
				2 : begin
						niovern = 1.d0 / (1.d0 + 1.d0/n1overn0 + n2overn1)
						ui = u2
					end
				else : begin
							print, 'Unknown ionization state'
							stop
						end
			endcase
			
			opacityL[i] = OPA * gf / ui * exp(-Elow/(PK * T[i])) * $
				(1.d0 - exp(-PH/PK * nu0 / T[i])) * niovern * nH[i] * ab[element-1]
				
			dnu = nu0 * doppler_width[i] / PC
			
			doppler_shift = nu0 * vmac[i] / PC
			
			v = (nu-nu0) / dnu
			va = doppler_shift / dnu
			prof = voigt(damping[i], v-va)
					
; Stokes I
			eta[*,i] = opacityC[i] + opacityL[i] * prof / (dnu * sqrt(!DPI))
			etaC[*,i] = opacityC[i]
			
			epsilon[*,i] = plancknu(T[i],nu) * eta[*,i]
			epsilonC[*,i] = plancknu(T[i],nu) * etaC[*,i]
					
		endfor
				
		if (order eq 0) then begin
; Boundary condition
			I0 = plancknu(T[0],nu)
			Intensity = shortcar_poly(z, eta, epsilon, mu, I0)	
		endif else begin
; Boundary condition
			I0 = plancknu(T[n-1],nu)
			Intensity = shortcar_poly(reverse(z), reverse(eta,2), reverse(epsilon,2), mu, I0)
			Intensity0 = shortcar(reverse(z), reverse(reform(etaC[0,*])), reverse(reform(epsilonC[0,*])), mu, I0[0])
		endelse
		
		spec[j].lambda = ptr_new(lambda)
		spec[j].specI = ptr_new(Intensity)
		spec[j].contI = Intensity0		
		
	endfor

	return, spec
end