;-----------------------
	
	function kh_ff,t,pe,rlambda

; PURPOSE:
;	It computes the negative hydrogen free-free continuum 
;	absorption coeff. per h atom (cm**2) as given by
;
;		john: 1989, a&a 193, 189
;
; INPUT:
;	t	temperature (k) (.lt.10080k & .gt.1400k)
;	pe	electron pressure
;	rlambda	wavelenght (a) (gt.1880a)
;
; OUTPUT:
;	kh_ff	continuum apsorption coefficient
; MODIFICATION HISTORY:
;
;	(started:	?
;	finished:	20/2/90)	for sos (fortran)
;			Sept. 1st, 95	fortran > idl (jorge)


; parameters for h- ff
	a1=[0.d0,2483.346d0,-3449.889d0,2200.04d0,-696.271d0,88.283d0]
	b1=[0.d0,285.827d0,-1158.382d0,2427.719d0,-1841.4d0,444.517d0]
	c1=[0.d0,-2054.291d0,8746.523d0,-13651.105d0,8624.97d0,$
              -1863.864d0]
	d1=[0.d0,2827.776d0,-11485.632d0,16755.524d0,-10051.53d0,$
              2095.288d0]
	e1=[0.d0,-1341.537d0,5303.609d0,-7510.494d0,4400.067d0,$
              -901.788d0]
	f1=[0.d0,208.952d0,-812.939d0,1132.738d0,-655.02d0,132.985d0]
	a2=[518.1021d0,473.2636d0,-482.2089d0,115.5291d0]
	b2=[-734.8666d0,1443.4137d0,-737.1616d0,169.6374d0]
	c2=[1021.1775d0,-1977.3395d0,1096.8827d0,-245.649d0]
	d2=[-479.0721d0,922.3575d0,-521.1341d0,114.243d0]
	e2=[93.1373d0,-178.9275d0,101.7963d0,-21.9972d0]
	f2=[-6.4285d0,12.36d0,-7.0571d0,1.5097d0]
; kh_ff
;	 wavelength from a to microns
	lambda=rlambda/1d4 
	theta=5040./t
;
	if(lambda lt .3645)then begin
	   com=a2*lambda^2.+b2+c2/lambda+d2/lambda^2.+$
     		e2/lambda^3.+f2/lambda^4.
	   e=(findgen(4)+2.)/2.
	   kh_ff=theta^e(0)*com(0) 
	   for i=1,3 do kh_ff=kh_ff+theta^e(i)*com(i) 
	endif else begin
	   com=a1*lambda^2.+b1+c1/lambda+d1/lambda^2.+$
     		e1/lambda^3.+f1/lambda^4.
	   e=(findgen(6)+2.)/2.
	   kh_ff=theta^e(0)*com(0) 
	   for i=1,5 do kh_ff=kh_ff+theta^e(i)*com(i) 
	endelse
	kh_ff=1e-29*(kh_ff*pe)
;
	return,float(kh_ff)
	end
