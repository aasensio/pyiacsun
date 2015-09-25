;-----------------------

	function kh,t,rlambda


; PURPOSE:
; 	neutral hydrogen atom bound-free (+free-free, I thing
;	by comparison with MIHALAS 1967,  Methods in Computational
; 	Physics, eq 78).	continuum absorption coeff.
; 	per h atom (cm**2) as given by
;
; 	landi degl'innocenti: 1976, astron. astro. supl. ser. 25, 379
; input:
;	t	temperature (k) (.lt.12000k)
;	rlambda	wavelenght (a)
;
; output:
;	kh	continuum apsorption coefficient
;
; MODIFICATION HISTORY:
;		adapted from khbf	Feb, 89 (jorge)


; definitions
	r=1.096776d-3
	c1=1.5777216d5
	c2=1.438668d8
     	const=1.045d-26
; let's go
	theta1=c1/t		
	theta3=theta1*2
	theta2=c2/t/rlambda
;	lowest level which can be photoionized
	n0=1+fix(sqrt(rlambda*r)) 
;	sum over states which can be photoionized
	if(n0 le 8)then begin
	   sum=exp(theta1/(n0*n0))/(n0*n0*n0)
	   for i=n0+1,8 do begin
		sum=sum+exp(theta1/(i*i))/(i*i*i)
	   endfor
	   sum=sum+(exp(theta1/81.)+.117)/theta3
	endif else begin
	   sum=(exp(theta1/n0^2.)+.117)/theta3
	endelse
;
; approximate g_ff taken from mihalas eq (80) @ theta=1,x=.5
;
	kh=(1-exp(-theta2))*exp(-theta1)*(rlambda*rlambda*rlambda)
	kh=const*kh*sum
;
	return,float(kh)
	end
