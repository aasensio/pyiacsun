@synth_lte
; Sample file that shows how to synthesize Stokes profiles in an atmosphere in LTE
; and putting the model in hydrostatic equilibrium
; The plot makes use of the Coyote Graphics library. If you don't have it (and you should),
; just change the last line
pro sample

	PK = 1.3806503d-16
	UMA = 1.66053873d-24
	PE = 4.8032d-10
	PC = 2.99792458d10
	PME = 9.10938188d-28
	PH = 6.62606876d-27
	
; Read the model atmosphere
	nlines = file_lines('hsra.model')-2
	openr,2,'hsra.model'
	t = ''
	readf,2,t
	readf,2,t
	model = fltarr(8,nlines)
	readf,2,model
	close,2
	n = n_elements(model[0,*])
			
	atmosphere = create_struct('ltau500',reform(model[0,*]), 'T',reform(model[1,*]), 'vmic',reform(model[3,*] * 1.d5), 'vmac',replicate(0.d0,n))
				
; Read the spectral lines
	nlines = file_lines('lines')-1
	openr,2,'lines'
	readf,2,t
	lines = fltarr(11,nlines)
	readf,2,lines
	close,2
	
	line = create_struct('element', reform(lines[0,*]), 'ioniz', reform(lines[1,*]), $
		'lambda0', reform(lines[2,*]), 'nu0', reform(PC / (lines[2,*]*1.d-8)), $
		'Elow', reform(lines[3,*] * PH * PC), 'gf', reform(10.d0^lines[4,*]),$
		'alpha_ABO', reform(lines[5,*]), 'sigma_ABO', reform(lines[6,*]), $
		'lambda_from', reform(lines[7,*]), 'lambda_to', reform(lines[8,*]), $
		'lambda_step', reform(lines[9,*]), 'gbar', reform(lines[10,*]))
		
	spectral = create_struct('lambda', 0L)

	Stokes = synth_lte(1.d0, atmosphere, line, spectral, Ic)
	

	cgplot, (*Stokes[0].lambda), (*Stokes[0].specI) / Stokes[0].contI
	cgoplot, (*Stokes[1].lambda), (*Stokes[1].specI) / Stokes[1].contI
	
	
	stop
end
