
	pro test

; PURPOSE:
; 	It allows to compare the absorption coefficients
;	computed by KAPPA_C with that listed in the HSRA
;	paper.
; def
	; HSRA data in ~jos/synthesis/model_set.pro
	common HSRA,nsh,tauhs,ths,pehs,pghs,kabchs,zz,rhohs
	tablas_hsra
	rlambda=5000.
	atmdat,1
	time=systime(1)
	kappa=KAPPA_C(pehs,ths,rlambda,htoverv,rmu)
	stop
	print,'time=',systime(1)-time
	;plot_oo,tauhs,kappa/rmu,xtitle='!7s!d!65!n',ytitle='!7j!d!6c5!6!n',$
	;	title='_____ HSRA          ........... KAPPA_C',$
;	;	charsize=1.5
	;oplot,tauhs,kabchs,line=1
	yy=1./(kabchs*rhohs)*(kappa*htoverv)
	;yy=abs((1-yy)/(1+yy))
	;tit='!6(!7j!6-!7j!d!6HSRA!n)/(!7j!6+!7j!d!6HSRA!n)'
	tit='!6!7j!6/!7j!d!6HSRA!n'
	plot_oi,tauhs,yy,xtitle='!7s!d!65!n',ytitle=tit,$
		charsize=1.5
	end

;	
;;;;;;;;;;;;;;;;;;;

	function KAPPA_C,pe,t,rlambda,htoverv,rmu,gp=gp,verbose=verbose,scatt=scatt

; PURPOSE:
; 	Provides the  continuum absorption coefficient per h particle 
;	(including ions and molecules)
;
; 	The two final parameters allow the computation of cont. abs. 
;	coeff. per gram and per unit volume.
;
; DESCRIPTION:
; it takes into account as source of opacity:
;
; 	negative hydrogen ion(bf+ff), neutral h, rayleigh scattering on h,
;	rayleigh scattering on h2 molecules, thompson scattering on electrons
;	total magnesium (bf+ff)
; 	(see references in the subroutines)
;
; INPUT:
;	pe	electronic pressure (dyne/cm**2) (it can be an array!)
;	t	temperature (k)	(it can be an array!)
;	rlambda	wavelegnth (a)	(single wavelength)
;	
; OUTPUT:
;	kappa_c	absorption coefficient per h particle: []=cm**2/nt
;	htoverv number of h particles per unit volume
;		[kappa_c*ntoverv]=cm**2/cm**3=cm**-1
;	rmu	mass per h particle
;		[kappa_c/rmu]=cm**2/gr
;	gp	total pressure (cgs), which is required by kappa_l
;		(it is a keyword!)
; REQUIRED SUBROUTINES: (* means enclosed here) 
;	GASPRESS(*),ELECTRON(*), SAHA(*)
;	ATMDAT, PARTITION
;
; MODIFICATION HISTORY:
; 		started:	20/2/90		jorge sanchez almeida
;		finished:	6/3/90		(This two dates correspon to the
;						FORTRAN version used by SOS)
;		Sep 95	translated to IDL (Jorge)
;			the contribution of 
;			has been commented outi; its contribution
;			is not significant and is highly time consuming
;			its computations
;		Oct 95	Electrons are just supplied  by the
;			elementes considered by HSRA


; definition
;	boltzmann constant
	 atmdat,1
	kb=1.38062d-16	
;
;	mass per h particle; form hsra model atmosphere 
;	gingerich et al.: 1971,
;	solar phys. 18,347
	;mu=2.368d-24
	mu=2.3848d-24;, according to ATMDAT
;
; let's start
	GASPRESS,pe,t,ph,phmi,phpl,ph2,ph2pl,gp
	
;	total h over neutral h
	factor=1+(phmi+phpl+2*(ph2+ph2pl))/ph	
; number of total h/volume
	htoverv=factor*ph/t/kb
; gr per total h (per h particle)
	rmu=mu
;
	d1=kh_ff(t,pe,rlambda)/factor
	d2=kh_bf(t,pe,rlambda)/factor
	;d3=khbf(t,rlambda)/factor
	d3=kh(t,rlambda)/factor
	d4=st(pe,ph)/factor
	d5=srh(rlambda)/factor
	d6=srh2(rlambda,ph,ph2)/factor
	d7=kmg(t,pe,rlambda) 
	
; Return background opacity and scattering separately if indicated	
	if (arg_present(scatt)) then begin
		scatt = d4+d5+d6
		kappa_c = d1+d2+d3+d7
	endif else begin
		kappa_c=d1+d2+d3+d4+d5+d6+d7
	endelse
	
	if (keyword_set(verbose)) then begin
	 	  print, 'PH, PH-, PH+, PH2, PH2+, PG'
		  print, ph, phmi, phpl, ph2, ph2pl, gp
	 	  print, 'H- bound-free                       : ', d2
		  print, 'H- free-free                        : ', d1
		  print, 'H bound-free + free-free            : ', d3
		  print, 'Thompson scattering                 : ', d4
		  print, 'Rayleigh scattering by H atoms      : ', d5
		  print, 'Rayleigh scattering by H2 molecules : ', d6
		  print, 'Bound-free + free-free of Mg        : ', d7
		  print, 'Total                               : ', kappa_c
	endif
;
	;print,'test whether kappa_c is double precision and why?
	;stop
	return,kappa_c
	end

;;-----------------------
;
;	function khbf,t,rlambda
;
;
;; PURPOSE:
;; 	neutral hydrogen atom bound-free continuum absorption coeff.
;; 	per h atom (cm**2) as given by
;;
;; 	landi degl'innocenti: 1976, astron. astro. supl. ser. 25, 379
;;		(except for the sum over states)
;; input:
;;	t	temperature (k) (.lt.12000k)
;;	rlambda	wavelenght (a)
;;
;; output:
;;	khbf	continuum apsorption coefficient
;;
;; MODIFICATION HISTORY:
;; 	(started:	20/2/90
;; 	finished:	21/2/90)	fortran, for sos
;;			Sep 1st, 95	fortran > idl (jorge)
;
;
;; definitions
;	r=1.096776d-3
;	c1=1.5777216d5
;	c2=1.438668d8
;     	const=1.045d-26
;; let's go
;	theta1=c1/t		
;	theta3=theta1*2
;	theta2=c2/t/rlambda
;;	lowest level which can be photoionized
;	n0=1+fix(sqrt(rlambda*r)) 
;;	sum over states which can be photoionized
;	if(n0 le 8)then begin
;	   sum=exp(theta1/(n0*n0))/(n0*n0*n0)
;	   for i=n0+1,8 do begin
;		sum=sum+exp(theta1/(i*i))/(i*i*i)
;	   endfor
;	   ;sum=sum+(exp(theta1/72.25d0)-1.)/theta3
;	   sum=sum+exp(theta1/72.25d0)/theta3
;	endif else begin
;	   ;sum=(exp(theta1/(n0-.5)^2.)-1.)/theta3
;	   sum=exp(theta1/(n0-.5)^2.)/theta3
;	endelse
;;
;	khbf=(1-exp(-theta2))*exp(-theta1)*(rlambda*rlambda*rlambda)
;	khbf=const*khbf*sum
;;
;	return,float(khbf)
;	end
;








