;--------------------------

	pro GASPRESS,pe,t,ph,phmi,phpl,ph2,ph2pl,pg

; PURPOSE:
; it gives the partial pressure of h, h-, h+, h2, h2+ as well as the
; total gas pressure. it requires information about atmospheric
; abundances which is provided by the subroutine atmdat
;
; references: mihalas: 1967, in "methods in computational physics"
;		      alder (ed.), academic press, p.1
;	      wittmann: 1974, solar phys. 35, 11		
;
; input:
;	pe	electron pressure (dy*cm**(-2))
;	t	temperature	(k)
;
; output:
;	ph	partial pressure of h (dny/cm2)
;	phmi	   "       "     "  h-   "
;	phpl	   "       "     "  h+   "
;	ph2	   "       "     "  h2   "
;	ph2pl	   "       "     "  h2+  "
;	pg	total gas pressure
;
;
; started:	21/2/90
; finised:	22/2/90
;
; it requires the subroutines: electron, saha, atmdat, partition
;
; MODITIOCATION HISTORY:
;	Sept 1st, 95	fortran > idl (jorge)
;	Included updated dissociation constant from
;		
;	

;
; parameters for h2, h2+, h+ and h-
	ch2=[12.533505d0,-4.9251644d0,5.6191273d-2,-3.2687661d-3]
	ch2pl=[1.1206998d1,-2.7942767d0,-7.9196803d-2,2.4790744d-2]
	chpl=[-13.595d0,2.5d0,-0.4772d0]
	chmi=[-0.747d0,2.5d0,.1249d0]
;	total abundance by number (it is not that of atmdat but comes from
;	that adopted in the hsra model atmosphere; gingerich et al.: 1971,
;	solar phys. 18, 347)
	totabu=.101d0
;
; let's go with the equilibrium constants
	theta=5040./t
;
	coh2=ch2(0)+(ch2(1)+(ch2(2)+ch2(3)*theta)*theta)*theta	
	coh2=10.^coh2
;
	coh2pl=ch2pl(0)+(ch2pl(1)+(ch2pl(2)+ch2pl(3)*theta)*theta)*theta
	coh2pl=10.^coh2pl
;	just the saha equation
	cohpl=chpl(0)*theta+chpl(1)*alog10(t)+chpl(2)
	cohpl=10.^cohpl
;	just the saha equation
	cohmi=chmi(0)*theta+chmi(1)*alog10(t)+chmi(2)
	cohmi=10.^cohmi
; now it solves a system of six equations with six unknows.
; 	contribution to pe from elements which are not h (two ionizations)
; 	(the relevant elemets have been taken from mihalas(1967))
;	Warning, I've commented out those elements which 
;	are not used for the HSRA (table II, SPh 18, 357)
	g1=ELECTRON(2,t,pe)	
;					! he
	g1=g1+ELECTRON(6,t,pe)	
;					! c
	g1=g1+ELECTRON(7,t,pe)	
;					! n
	g1=g1+ELECTRON(8,t,pe)	
;					! o
	g1=g1+ELECTRON(11,t,pe)	
;					! na
	g1=g1+ELECTRON(12,t,pe)	
;					! mg
	g1=g1+ELECTRON(13,t,pe)	
;					! al
	g1=g1+ELECTRON(14,t,pe)	
;					! si
	g1=g1+ELECTRON(16,t,pe)	
;					! s
	g1=g1+ELECTRON(19,t,pe)	
;					! k
	g1=g1+ELECTRON(20,t,pe)	
;					! ca
	g1=g1+ELECTRON(24,t,pe)	
;					! cr
	g1=g1+ELECTRON(26,t,pe)	
;					! fe
;	
	g2=cohpl/pe	
	g3=pe/cohmi
	g4=pe/coh2pl
	g5=pe/coh2
;
	a=1+g2+g3
	b=2*(1+g2*g4/g5)
	c=g5
	d=g2-g3
	e=g2*g4/g5
;
	c1=c*b*b+a*d*b-e*a*a
	c2=2*a*e-d*b+a*b*g1
	c3=(-1)*(e+b*g1)
;
	f1=(-1)*c2/2/c1
	f2=sqrt((c2/2/c1)^2.-c3/c1)
	pos=where(c1 ge 0.)
	if(pos(0) ne -1)then f1(pos)=f1(pos)+f2(pos)
	pos=where(c1 lt 0.)
	if(pos(0) ne -1)then f1(pos)=f1(pos)-f2(pos)
	f2=g2*f1
	f3=g3*f1
	f5=(1-a*f1)/b
	f4=e*f5
	f6=pe/(f2-f3+g1+f4)
;
	ph=float(f1*f6)
	phpl=float(f2*f6)
	phmi=float(f3*f6)
	ph2pl=float(f4*f6)
	ph2=float(f5*f6)
	pg=pe+float(f6*(totabu+(f1+f2+f3+f4+f5)))
;
	end	
