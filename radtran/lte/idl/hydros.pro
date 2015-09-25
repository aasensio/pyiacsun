;+++++++++++++
;
	pro HYDROS,xx,temp,pe,pg
;
; PURPOSE:
;	To provide the gas pressure and electron
;	pressure which are in hydrostatic equilibrium (HE)
;	with a run temperature(log(tau)).
;	Right now, it is just valid for the sun.
;
; INPUT:
;	xx	integration grid (xx=log(tau5)!!!).
;		Inportant: the grid may not be
;		equispaced.
;	temp	temperatures
;
; OUTPUT:
;	pe	electron pressure in HE
;	pg	gas (total) pressure in HE
;		To check the HE condition,
;		deriv(xx,alog10(pg)) = g*tau/pg/k_c
;
; NOTES:
;	- INPORTANT: There is no single hydrostatic
;	equilibrium pe(tau) and pg(tau). It depends
;	(specially in the upper layers) on the
;	bounday condition used to solve the
;	differential equation!!!!!

;	-The hydrostatic equilibrium equation
;	has been rewritten as an equation for 
;	the electron pressure (i.e., it is not
;	using the iterative techniques described
;	in textbooks like Gray 1992, The observations
;	an analysis of stellar atmospheres).
;
;	-The DE is integrated using a simple
;	point-slope formula (as suggested by Basilio).
;	There is an 'old' version which uses a 
;	predictor corrector method, but it is slower
;	and, if required, one can gain more precision
;	increasing the points of the grid in this
;	routine.
;	
;	-Electrons are supplied by the atoms used
;	by the HSRA model atmosphere (table II, SPh 18, p. 357)
;	
;	-Abundances and the rest of atomoc parameters
;	are taken from atmdat
;	
; MODIFICATION HYSTORY:
;	Oct 31,95	Jorge
;+
;
;
; atomic data
	common ATMDAT,wgt,abu,ei1,ei2,sym
;
;	elements which release electrons (according
;	to the HSRA model atmosphere; SPh 18, p. 357)
	;H, He, C, O, Na, Mg, Al, Si, S, Fe
	nel=[1,2,6,8,11,12,13,14,16,26]
	nnel=n_elements(nel)
	;aahs=[1.,0.1,3.55e-4,5.9e-4,1.51e-6,3.02e-5,2.51e-6,3.55e-5,$
	;	1.62e-5,3.16e-5] & aahs=aahs/float(total(aahs))
;
; sorting the arrays to match that
;	required for the integration
	order=sort(xx)
	x=xx(order)
	t=temp(order)
	step=x(1)-x(0)
;
; defining the various arrays
	nt=n_elements(t)
	pe=fltarr(nt)
	pg=fltarr(nt)
	phi=fltarr(nnel,nt)
	tdphi=fltarr(nnel,nt)
	one=replicate(1.,nt)
;
	for i=0,nnel-1 do begin
	   PARTITION,t,nel(i),u1,u2,u3
	   phi(i,*)=saha(t,one,u1,u2,ei1(nel(i)-1))
	   tt=t*1.001
	   tdphi(i,*)=saha(tt,one,u1,u2,ei1(nel(i)-1))
	endfor
	tdphi=(tdphi-phi)/.001
	aa=abu(nel-1)/1.10165;=total(abu)
	;print,'still HSRA abu'
	;aa=aahs
;
	dtdx=deriv(xx,alog10(t))
;
;
; let's integrate the differential equation
;	for dlog(pe)/dlog(tau)
	y=fltarr(nt)
	dy=fltarr(nt)
	dz=fltarr(nt)
	tau=10.^x 
	g=2.7398e4;cm/s^2
;
;	setting initial values
	pe0=0.001
	k_c=kappa_c(pe0,t(0),5000.,x1,mu)/mu
	pe(0)=sqrt(g*tau(0)/k_c*pe0)
	delta=abs(pe(0)/pe0-1.)
	pe0=pe(0)
	index=0
	while delta ge 0.1 do begin 
	   pe2=pe(0)/2.
	   pg(0)=total(aa*phi(*,0)/(phi(*,0)+pe2))
	   pg(0)=pe2*(1.+1./pg(0))
	   dpgdt=total(aa*tdphi(*,0)/(phi(*,0)+pe2)^2.)
	   dpgdt=(pg(0)-pe2)^2*dpgdt*dtdx(0)
	   der=total(aa*phi(*,0)/(phi(*,0)+pe2)^2.)
	   der=pg(0)+(pg(0)-pe2)^2*der
	   k_c=kappa_c(pe2,t(0),5000.,x1,mu)/mu
	   dy(0)=(g*tau(0)/k_c+dpgdt)/der
	   pe0=pe(0)
	   if(dy(0) le 0.)then pe(0)=pe0/2. else $
	      pe(0)=(dy(0)+1.)*pe2
	   delta=abs(pe(0)/pe0-1.)
	   index=index+1
	   if(index eq 20)then delta=0.
	endwhile
	y(0)=alog10(pe(0))
	if(index eq 20)then print,'HYDROS: pe(0) might be wrong'
;
;	LET'S go
	for i=1,nt-1 do begin
	   y(i)=y(i-1)+(xx(i)-xx(i-1))*dy(i-1)
	   pe(i)=10.^y(i)
	   pg(i)=total(aa*phi(*,i)/(phi(*,i)+pe(i)))
	   pg(i)=pe(i)*(1.+1./pg(i))
	   dpgdt=total(aa*tdphi(*,i)/(phi(*,i)+pe(i))^2.)
	   dpgdt=(pg(i)-pe(i))^2*dpgdt*dtdx(i)
	   der=total(aa*phi(*,i)/(phi(*,i)+pe(i))^2.)
	   der=pg(i)+(pg(i)-pe(i))^2*der
	   k_c=kappa_c(pe(i),t(i),5000.,x1,mu)/mu
	   dy(i)=(g*tau(i)/k_c+dpgdt)/der
	endfor
;
; return with the right order
	pg=pg(order)
	pe=pe(order)
	end
;+++++++++++++
;
;
	pro TEST
;
; PURPOSE: To test HYDROS
; data from HSRA
	common HSRA,nhs,tauhs,ths,pehs,pghs,k_cm,zm
	common CHAPMAN,taum,tm,pem
	common HOLMU,xhm,pehm,thm
	common SOLPLA,ltausp,tsp,pesp,bsp
	common SOLNET,ltausn,tsn,pesn,bsn
	common MALTBY,tauma,tma,nema
;
; selecting model
	print,'Which model? 1=HSRA(def) 2=CHAPMAN  3=HOLMU'
	print,'             4=SOL-PLAGE 5=SOL-NET 6=MAL-SPOT'
	read,index
	case index of
;
1:	begin
	tablas_hsra
	tm=ths
	pem=pehs
	xm=alog10(tauhs)
	pgm=pghs
	title='HSRA'
	end
;
2:	begin
	tablas_cha
	xm=alog10(taum)
	title='CHAPMAN'
	end
;
3:	begin
	tablas_holmu
	xm=xhm
	pem=10.^pehm
	tm=thm
	title='HOLMU'
	end
;
4:	begin
	tablas_solpla
	xm=ltausp
	tm=tsp
	pem=pesp
	title='S-PLAGE'
	end
;
5:	begin
	tablas_solnet
	xm=ltausn
	tm=tsn
	pem=pesn
	title='S-NETWORK'
	end
;
6:	begin
	k_b=1.3806e-16
	tablas_mal
	pem=k_b*nema*tma
	xm=alog10(tauma)
	tm=tma
	title='M-UMBRA'
	end
;
else:	begin	
	print,'No such model'
	return
	end
	endcase
;

;
; range for the integration
	xmin=-4.0
	xmax=1.
	nx=100
	xx=findgen(nx)/(nx-1.)*(xmax-xmin)+xmin
	t=spline(xm,tm,xx)
	ppe=spline(xm,pem,xx)
	pos=where(xx lt xm(0)) 
	if(pos(0) ne -1)then t(pos)=tm(0)
	tau=10^xx
;
	g=2.7398e4;cm/s^2
	atmdat,1

;
	time0=systime(1)
	HYDROS,xx,t,pe,pg
	print,'time(s)=',float(systime(1)-time0)
; plots
;
plot:
	!p.charsize=1.5
	!p.multi=[0,3,1]
	k_c=kappa_c(pe,t,5000.,x1,mu)/mu
	plot_io,xx,g*tau/pg/k_c
	oplot,xx,deriv(xx,alog10(pg)),line=1
;
	plot_io,xx,pg,line=1,ytitle='Pg',xtitle='log(tau5)'
	if(n_elements(pgm) eq 0) then pgm=tm+1e15
	oplot,xm,pgm 
;
	;plot,xx,ppe/pe-1,$
	;	ytitle='Pe/Pe!dHYDROS!n-1',xtitle='log(tau5)',$
	;	xrange=[xx(0),xx(nx-1)],title=title
	plot_io,xx,ppe,ytitle='Pe',xtitle='log(tau5)',$
		xrange=[xx(0),xx(nx-1)],title=title
	oplot,xx,pe,line=1
	if(n_elements(first) eq 0)then begin
		first=1
		set_plot,'ps'
		device,/landscape
		goto,plot
	endif else begin
		device,/close
		set_plot,'x'
	endelse
	!p.multi=0
;
	end

