;--------------------

      pro PARTITION,tt,nel,u1,u2,u3


; PURPOSE:
;     It supplies the partition function of elements (hydrogen to uranium)
;     It was adapted from the routine ATMDAT by
;     A. D. Wittmann, goettingen (1975).(1974, solar phys. 35, 11)
;
; 	In fact the original routine has been split into two, ATMDAT,
;	which provides real atomica parameters and PARTITION,
;	which computes partition functions (which depends on temperature
;	and pressure).
;
; input:
;	tt	temperature(k) (this parameter can be an array!!!!)
;	nel	atomic number	(just one single element!!!!)
;	
; output:
;	u1	partition function, neutral element
;	u2	    "         "   , 1st ionization
;	u3          "         "   , 2nd     "
;		(source: bolton: 1970, astrophys. j. 161, 1187
;			 aller & everett: 1972, astrophys. j. 172, 447
;			 and....?)
; MODIFICATION  HISTORY:
; started:	1975, wittmann
; finised:	22/2/90, just adding comments (j. sanchez almeida)	
;		6/3/90   change of abundances to that quoted by
;			 grevesse: 1984, phys. scripta t8, 49
;		August 31st, 95		translated to IDL form 
;					fortran (sos). Jorge
;
;		March 2000 Patch H and He at moderate-high temperatures


; errors
	if(n_elements(nel) gt 1) then stop,'PARTITION: nel > 1'
;
	
      t=tt
      nt=n_elements(t)
      x=alog(5040./t)
      y=1.e-3*t
;	CASE goes just to one point
	case nel of
;
;	! h
1:	begin
	; from Irwin, 1981, ApJS, 45, 621
	coeff=[-2.61655891d2,1.63428326d2,-4.06133526d1,5.03282928d0,-3.10998364d-1,$
	   7.66654594d-3]
	zz=alog(tt)
	u1=(((((coeff(5)*zz+coeff(4))*zz)+coeff(3))*zz+coeff(2))*zz+coeff(1))*zz+coeff(0)
	u1=exp(u1)
	pos=where(tt gt 16000.)
	if(pos(0) ne -1)then u1(pos)=exp(7.19420668d-1)
	u2=replicate(1.,nt)
	u3=replicate(0.,nt)
	;
	;	original Wittmans code
	;u1=2.0+0.0*t 
	;pos=where((t gt 1.3e4) and (t lt 1.62e4))
	;if(pos(0) ne -1) then u1(pos)=1.51+3.8e-5*t(pos)  
	;pos=where(t gt 1.62e4)
	;if(pos(0) ne -1)then u1(pos)=11.41d0+t(pos)*(-1.1428d-3+t(pos)*3.52d-8)
        ;u2=1.000+0.*t
      	;u3=0.000+0.*t
	end
;
;	! he
2:	begin
	u1=1.+0.*t
      	;pos=where(t gt 3e4)
      	;if(pos(0) ne -1)then u1(pos)=14.8+t(pos)*(-9.4103e-4+t(pos)*1.6095e-8)
      	u2=2.000+0.*t
      	u3=1.000+0.*t
	end
;
;	! Li
3:     	begin
	u1=2.081-y*(6.8926e-2-y*1.4081e-2) 
      	pos=where(t gt 6.e3) 
      	if(pos(0) ne -1)then u1(pos)=3.4864+t(pos)*(-7.3292e-4+t(pos)*8.5586e-8)
      	u2=1.+0.*t
      	u3=2.+0.*t
	end
;
;	! Be
4:     	begin
	u1=(1.+0.*t) > (.631+7.032e-5*t)  ; IDL choose the largest at each 
      	u2=2.+0.*t
      	u3=1.+0.*t
	end
;
;	! B
5:     	begin
	u1=5.9351+1.0438e-2*y 
      	u2=1.+0.*t
      	u3=2.+0.*t
	end
;
;	! C
6:     	begin
	u1=8.6985+y*(2.0485e-2+y*(1.7629e-2-3.9091e-4*y)) 
	pos=where(t gt 1.2e4)
	if(pos(0) ne -1)then u1(pos)=13.97+t(pos)*(-1.3907e-3+t(pos)*9.0844e-8)
      	u2=5.838+1.6833e-5*t
	pos=where(t gt 2.4e4)
      	if(pos(0) ne -1)then u2(pos)=10.989+t(pos)*(-6.9347e-4+t(pos)*2.0861e-8)
      	u3=1.+0.*t
      	pos=where(t gt 1.95e4)
      	if(pos(0) ne -1)then u3(pos)=-0.555+8e-5*t(pos) else u3=1.00
	end
;
;	! N
7:	begin
	u1=3.9914+y*(1.7491e-2-y*(1.0148e-2-y*1.7138e-3)) 
      	pos=where(t gt 8800. and t le 1.8e4)
        if(pos(0) ne -1)then u1(pos)=2.171+2.54e-4*t(pos)
	pos=where(t gt 1.8e4)
      	if(pos(0) ne -1)then u1(pos)=11.396+t(pos)*(-1.7139e-3+t(pos)*8.633e-8)
      	u2=8.060+1.420e-4*t
      	pos=where(t gt 3.3e4) 
      	if(pos(0) ne -1)then u2(pos)=26.793+t(pos)*(-1.8931e-3+t(pos)*4.4612e-8)
      	u3=5.9835+t*(-2.6651e-5+t*1.8228e-9)
      	pos=where(t lt 7310.5)
	if(pos(0) ne -1)then u3(pos)=5.89
	end
;
;	! O
8:	begin
	u1=8.29+1.10e-4*t 
      	pos=where(t gt 1.9e4)
	if(pos(0) ne -1)then u1(pos)=66.81+t(pos)*(-6.019e-3+t(pos)*1.657e-7)
      	u2=(4.+0.*t) > (3.51+8.e-5*t)
      	pos=where(t gt 3.64e4)
	if(pos(0) ne -1)then u2(pos)=68.7+t(pos)*(-4.216e-3+t(pos)*6.885e-8)
      	u3=7.865+1.1348e-4*t
	end
;
;	! F
9:	begin
	u1=4.5832+y*(.77683+y*(-.20884+y*(2.6771e-2-1.3035e-3*y))) 
      	pos=where(t gt 8750. and t le 2.e4)
	if(pos(0) ne -1)then u1(pos)=5.9+0.*pos
      	pos=where(t gt 2e4)
	if(pos(0) ne -1)then u1(pos)=15.16+t(pos)*(-9.229e-4+t(pos)*2.312e-8)
      	u2=8.15+8.9e-5*t
      	u3=2.315+1.38e-4*t
	end
;
;	! Ne
10:	begin
	u1=1.+0.*t
      	pos=where(t gt 2.69e4)
	if(pos(0) ne -1)then u1(pos)=26.3+t(pos)*(-2.113e-3+t(pos)*4.359e-8)
      	u2=5.4+4.e-5*t
      	u3=7.973+7.956e-5*t
	end
;
;	! Na
11:    	begin
	u1= (2.+0.*t) > (1.72+9.3e-5*t)
      	pos=where(t gt 5400. and t le 8.5e3)
	if(pos(0) ne -1)then u1(pos)=-0.83+5.66e-4*t(pos)
      	pos=where(t gt 8.5e3)
	if(pos(0) ne -1) then u1(pos)=$
		4.5568+t(pos)*(-1.2415e-3+t(pos)*1.3861e-7)
      	u2=1.+0.*t
      	u3=5.69+5.69e-6*t
	end
;
;	! Mg
12:	begin
	u1=1.+exp(-4.027262-x*(6.173172+x*(2.889176+x*(2.393895+.784131*x))))
      	pos=where(t gt 8e3)
	if(pos(0) ne -1)then u1(pos)=2.757+t(pos)*(-7.8909e-4+t(pos)*7.4531e-8)
      	u2=2.+exp(-7.721172-x*(7.600678+x*(1.966097+.212417*x)))
      	pos=where(t gt 2e4)
	if(pos(0) ne -1)then u2(pos)=7.1041+t(pos)*(-1.0817e-3+t(pos)*4.7841e-8)
      	u3=1.+0.*t
      	end
;
;	! Al
13:	begin
	u1=5.2955+y*(.27833-y*(4.7529e-2-y*3.0199e-3)) 
      	u2=(1.+0.*t) > (.725+3.245e-5*t)
      	pos=where(t gt 2.24e4)
	if(pos(0) ne -1)then u2(pos)=61.06+t(pos)*(-5.987e-3+t(pos)*1.485e-7)
      	u3=(2.+0.*t) > (1.976+3.43e-6*t)
      	pos=where(t gt 1.814e4)
	if(pos(0) ne -1)then u3(pos)=3.522+t(pos)*(-1.59e-4+t(pos)*4.382e-9)
	end
;
;	! Si
14:	begin
	u1=6.7868+y*(.86319+y*(-.11622+y*(.013109-6.2013e-4*y))) 
      	pos=where(t gt 1.04e4)
	if(pos(0) ne -1)then u1(pos)=86.01+t(pos)*(-1.465e-2+t(pos)*7.282e-7)
      	u2=5.470+4.e-5*t
      	pos=where(t gt 1.8e4)
	if(pos(0) ne -1)then u2(pos)=26.44+t(pos)*(-2.22e-3+t(pos)*6.188e-8)
      	u3=(1.+0.*t) > (.911+1.1e-5*t)
      	pos=where(t gt 3.33e4)
	if(pos(0) ne -1)then u3(pos)=19.14+t(pos)*(-1.408e-3+t(pos)*2.617e-8)
	end
;
;	! P
15:	begin
	u1=4.2251+y*(-.22476+y*(.057306-y*1.0381e-3)) 
      	pos=where(t gt 6.e3)
	if(pos(0) ne -1)then u1(pos)=1.56+5.2e-4*t(pos)
      	u2=4.4151+y*(2.2494+y*(-.55371+y*(.071913-y*3.5156e-3)))
      	pos=where(t gt 7250.)
	if(pos(0) ne -1)then u2=4.62+5.38e-4*t(pos)
      	u3=5.595+3.4e-5*t
	end
;
;	! S
16:	begin
	u1=7.5+2.15e-4*t 
      	pos=where(t gt 1.16e4)
	if(pos(0) ne -1)then u1(pos)=38.76+t(pos)*(-4.906e-3+t(pos)*2.125e-7)
      	u2=2.845+2.43e-4*t
      	pos=where(t gt 1.05e4)
	if(pos(0) ne -1)then u2(pos)=6.406+t(pos)*(-1.68e-4+t(pos)*1.323e-8)
      	u3=7.38+1.88e-4*t
	end
;
;	! Cl
17:	begin
	u1=5.2+6.e-5*t 
      	pos=where(t gt 1.84e4)
	if(pos(0) ne -1)then u1(pos)=-81.6+4.8e-3*t(pos)
      	u2=7.0+2.43e-4*t
      	u3=2.2+2.62e-4*t
	end
;
;	! Ar
18:    	begin
	u1=1.+0.*t
      	u2=5.20+3.8e-5*t
      	u3=7.474+1.554e-4*t
	end
;
;	! K
19:	begin
	u1=1.9909+y*(.023169-y*(.017432-y*4.0938e-3)) 
      	pos=where(t gt 5800.)
	if(pos(0) ne -1)then u1(pos)=-9.93+2.124e-3*t(pos)
      	u2=1.+0.*t
      	u3=5.304+1.93e-5*t
	end
;
;	! Ca
20:	begin
	u1=1.+exp(-1.731273-x*(5.004556+x*(1.645456+x*(1.326861+.508553*x))))
      	u2=2.+exp(-1.582112-x*(3.996089+x*(1.890737+.539672*x)))
      	u3=1.000+0.*t
	end
;
;	! Sc
21:	begin
	u1=4.+exp(2.071563+x*(-1.2392+x*(1.173504+.517796*x))) 
      	u2=3.+exp(2.988362+x*(-.596238+.054658*x))
      	u3=10.+0.*t 
	end
;
;	! Ti
22:	begin
	u1=5.+exp(3.200453+x*(-1.227798+x*(.799613+.278963*x))) 
      	pos=where(t lt 5.5e3)
	if(pos(0) ne -1)then u1(pos)=16.37+t(pos)*(-2.838e-4+t(pos)*5.819e-7)
      	u2=4.+exp(3.94529+x*(-.551431+.115693*x))
      	u3=16.4+8.5e-4*t
	end
;
;	! V
23:	begin
	u1=4.+exp(3.769611+x*(-.906352+x*(.724694+.1622*x))) 
      	u2=1.+exp(3.755917+x*(-.757371+.21043*x))
      	u3=-18.+1.03e-2*t
      	pos=where(t lt 2.25e3)
	if(pos(0) ne -1)then u3(pos)=2.4e-3*t(pos)
	end
;
;	! Cr
24:	begin
	u1=7.+exp(1.225042+x*(-2.923459+x*(.154709+.09527*x))) 
      	u2=6.+exp(.128752-x*(4.143973+x*(1.096548+.230073*x)))
      	u3=10.4+2.1e-3*t
	end
;
;	!Mn
25:	begin
	u1=6.+exp(-.86963-x*(5.531252+x*(2.13632+x*(1.061055+.265557*x)))) 
      	u2=7.+exp(-.282961-x*(3.77279+x*(.814675+.159822*x)))
      	u3=10.+t*0.
	end
;
;	! Fe
26:	begin
	u1=9.+exp(2.930047+x*(-.979745+x*(.76027+.118218*x))) 
      	pos=where(t lt 4e3)
      	if(pos(0) ne -1)then u1(pos)=15.85+t(pos)*(1.306e-3+t(pos)*2.04e-7)
      	pos=where(t gt 9e3)
      	if(pos(0) ne -1)then u1(pos)=39.149+t(pos)*(-9.5922e-3+t(pos)*1.2477e-6)
      	u2=10.+exp(3.501597+x*(-.612094+.280982*x))
      	pos=where(t gt 1.8e4)
      	if(pos(0) ne -1)then u2(pos)=68.356+t(pos)*(-6.1104e-3+t(pos)*5.1567e-7)
      	u3=17.336+t*(5.5048e-4+t*5.7514e-8)
	end
;
;	! Co
27:	begin
	u1=8.65+4.9e-3*t 
      	u2=11.2+3.58e-3*t
      	u3=15.0+1.42e-3*t
	end
;
;	! Ni
28:	begin
	u1=9.+exp(3.084552+x*(-.401323+x*(.077498-.278468*x))) 
      	u2=6.+exp(1.593047-x*(1.528966+.115654*x))
      	u3=13.3+6.9e-4*t
	end
;
;	!Cu
29:	begin
	u1= (2.+0.*t) > (1.50+1.51e-4*t) 
      	pos=where(t gt 6250.)
	if(pos(0) ne -1)then u1(pos)=-.3+4.58e-4*t(pos)
      	u2= (1.+0.*t) > (.22+1.49e-4*t)
      	u3=8.025+9.4e-5*t
	end
;
;	!Zn
30:	begin
	u1= (1.+0.*t) > (.632+5.11e-5*t) 
      	u2=2.00+0.*t
      	u3=1.00+0.*t
	end
;
;	!Ga
31:	begin
	u1=1.7931+y*(1.9338+y*(-.4643+y*(.054876-y*2.5054e-3))) 
      	pos=where(t gt 6.e3)
	if(pos(0) ne -1) then u1(pos)=4.18+2.03e-4*t(pos)
      	u2=1.0+0.*t
      	u3=2.0+0.*t
	end
;
;	!Ge
32:	begin
	u1=6.12+4.08e-4*t 
      	u2=3.445+1.78e-4*t
      	u3=1.1+0.*t
	end
;
;	! As
33:	begin
	u1=2.65+3.65e-4*t 
      	u2=-.25384+y*(2.284+y*(-.33383+y*(.030408-y*1.1609e-3)))
      	pos=where(t gt 1.2e4)
	if(pos(0) ne -1) then u2(pos)=8.+0.*pos
      	u3=8.+0.*t 
	end
;
;	! Se
34:	begin
	u1=6.34+1.71e-4*t 
      	u2=4.1786+y*(-.15392+3.2053e-2*y)
      	u3=8.+0.*t 
	end
;
;	!Br
35:	begin
	u1=4.12+1.12e-4*t 
      	u2=5.22+3.08e-4*t
      	u3=2.3+2.86e-4*t
	end
;
;	! Kr
36:	begin
	u1=1.00 +0.*t
      	u2=4.11+7.4e-5*t
      	u3=5.35+2.23e-4*t
	end
;
;	! Rb
37:	begin
	u1=(2.+0.*t) > (1.38+1.94e-4*t) 
      	pos=where(t gt 6250.)
	if(pos(0) ne -1)then u1(pos)=-14.9+2.79e-3*t(pos)
      	u2=1.000+0.*t
      	u3=4.207+4.85e-5*t
	end
;
;	! Sr
38:	begin
	u1=.87127+y*(.20148+y*(-.10746+y*(.021424-y*1.0231e-3))) 
      	pos=where(t gt 6500.)
	if(pos(0) ne -1)then u1(pos)=-6.12+1.224e-3*t(pos)
      	u2=(2.+0.*t) > (.84+2.6e-4*t)
      	u3=1.0+0.*t 
	end
;
;	! Y
39:	begin
	u1=.2+2.58e-3*t 
      	u2=7.15+1.855e-3*t
      	u3=9.71+9.9e-5*t
	end
;
;	!Zc
40:	begin
	u1=76.31+t*(-1.866e-2+t*2.199e-6)  
      	pos=where(t lt 6236.)
	if(pos(0) ne -1)then u1(pos)=6.8+t(pos)*(2.806e-3+t(pos)*5.386e-7)
      	u2=4.+exp(3.721329-.906502*x)
      	u3=12.3+1.385e-3*t
	end
;
;	! Nb
41:	begin
	u1=(1.+0.*t) > (-19.+1.43e-2*t)
      	u2=-4.+1.015e-2*t
      	u3=25.+0.*t 
	end
;	Mo
42:	begin
	u1=(7.+0.*t) > (2.1+1.5e-3*t) 
      	pos=where(t gt 7.e3)
	if(pos(0) ne -1)then u1(pos)=-38.1+7.28e-3*t(pos)
      	u2=1.25+1.17e-3*t
      	pos=where(t gt 6900.)
	if(pos(0) ne -1)then u2(pos)=-28.5+5.48e-3*t(pos)
      	u3=24.04+1.464e-4*t
	end
;
;	Tc
43:	begin
	u1=4.439+y*(.30648+y*(1.6525+y*(-.4078+y*(.048401-y*2.1538e-3)))) 
      	pos=where(t gt 6.e3)
	if(pos(0) ne -1)then u1(pos)=24.+0.*pos
      	u2=8.1096+y*(-2.963+y*(2.369+y*(-.502+y*(.049656-y*1.9087e-3))))
      	pos=where(t gt 6.e3)
	if(pos(0) ne -1)then u2(pos)=17.+0.*pos
      	u3=220.+0.*t 
	end
;
;	Ru
44:	begin
	u1=-3.+7.17e-3*t 
      	u2=3.+4.26e-3*t
      	u3=22.+0.*t
	end
;
;
45:	begin
	u1=6.9164+y*(3.8468+y*(.043125-y*(8.7907e-3-y*5.9589e-4))) 
      	u2=7.2902+y*(1.7476+y*(-.038257+y*(2.014e-3+y*2.1218e-4)))
      	u3=30.+t*0. 
	end
;
;	Pd
46:	begin
	u1=(1.+0.*t) > (-1.75+9.86e-4*t) 
      	u2=5.60+3.62e-4*t
      	u3=20.+t*0. 
	end
;
;	Ag
47:	begin
	u1=(2.+t*0.) > (1.537+7.88e-5*t) 
      	u2=(1.+t*0.) > (0.73+3.4e-5*t)
      	u3=6.773+1.248e-4*t
	end
;
;	Cd
48:	begin
	u1=(1.+t*0.) > (.43+7.6e-5*t) 
      	u2=2.00+t*0.
      	u3=1.00+t*0.
	end
;
;	In
49:	begin
	u1=2.16+3.92e-4*t 
      	u2=1.0+t*0.
      	u3=2.00+t*0.
	end
;
;	Sn
50:	begin
	u1=2.14+6.16e-4*t 
      	u2=2.06+2.27e-4*t
      	u3=1.05+0.*t
	end
;
;	Sb
51:	begin
	u1=2.34+4.86e-4*t 
      	u2=.69+5.36e-4*t
      	u3=3.5+0.*t
	end
;
;	Te
52:	begin
	u1=3.948+4.56e-4*t 
      	u2=4.2555+y*(-.25894+y*(.06939-y*2.4271e-3))
      	pos=where(t gt 1.2e4)
	if(pos(0) ne -1)then u2(pos)=7.+0.*pos
      	u3=5.+0.*t 
	end
;
;	I
53:	begin
	u1=(4.+0.*t) > (3.8+9.5e-5*t) 
      	u2=4.12+3.e-4*t
      	u3=7.+0.*t
	end
;
;	Xe
54:	begin
	u1=1.00+0.*t 
      	u2=3.75+6.876e-5*t
      	u3=4.121+2.323e-4*t
	end
;
;	Cs
55:	begin
	u1=(2.+0.*t) > (1.56+1.67e-4*t) 
      	pos=where(t gt 4850.)
	if(pos(0) ne -1)then u1(pos)=-2.680+1.04e-3*t(pos)
      	u2=1.000+t*0.
      	u3=3.769+4.971e-5*t
	end
;
;	Ba
56:	begin
	u1=(1.+0.*t) > (-1.8+9.85e-4*t) 
      	pos=where(t gt 6850.)
	if(pos(0) ne -1)then u1(pos)=-16.2+3.08e-3*t(pos)
      	u2=1.11+5.94e-4*t
      	u3=1.00+0.*t
	end
;
;	La
57:	begin
	u1=15.42+9.5e-4*t 
      	pos=where(t gt 5060.)
	if(pos(0) ne -1)then u1(pos)=1.+3.8e-3*t(pos)
      	u2=13.2+3.56e-3*t
      	u3=12.+0.*t 
	end
;
;	Ce
58:	begin
	u1=9.+exp(5.202903+x*(-1.98399+x*(.119673+.179675*x))) 
      	u2=8.+exp(5.634882-x*(1.459196+x*(.310515+.052221*x)))
      	u3=9.+exp(3.629123-x*(1.340945+x*(.372409+x*(.03186-.014676*x))))
	end
;
;	Pr 
59:	begin
	u2=9.+exp(4.32396-x*(1.191467+x*(.149498+.028999*x))) 
      	u1=u2 
      	u3=10.+exp(3.206855+x*(-1.614554+x*(.489574+.277916*x)))
	end
;
;	nb
60:	begin
	u1=9.+exp(4.456882+x*(-2.779176+x*(.082258+x*(.50666+.127326*x)))) 
      	u2=8.+exp(4.689643+x*(-2.039946+x*(.17193+x*(.26392+.038225*x))))
      	u3=u2 
	end
;
;	Pm
61:	begin
	u1=20.+0.*t 
      	u2=25. +0.*t
      	u3=100. +0.*t
	end
;
;	Sm
62:	begin
	u1=1.+exp(3.549595+x*(-1.851549+x*(.9964+.566263*x))) 
      	u2=2.+exp(4.052404+x*(-1.418222+x*(.358695+.161944*x)))
      	u3=1.+exp(3.222807-x*(.699473+x*(-.056205+x*(.533833+.251011*x))))
	end
;
;	Eu
63:	begin
	u1=8.+exp(1.024374-x*(4.533653+x*(1.540805+x*(.827789+.286737*x))))
      	u2=9.+exp(1.92776+x*(-1.50646+x*(.379584+.05684*x)))
      	u3=8.+0.*t 
	end
;
;
64:	begin
	u1=5.+exp(4.009587+x*(-1.583513+x*(.800411+.388845*x))) 
      	u2=6.+exp(4.362107-x*(1.208124+x*(-.074813+x*(.076453+.055475*x))))
      	u3=5.+exp(3.412951-x*(.50271+x*(.042489-4.017e-3*x)))
	end
;
;	Tb
65:	begin
	u1=16.+exp(4.791661+x*(-1.249355+x*(.570094+.240203*x))) 
      	u2=15.+exp(4.472549-x*(.295965+x*(5.88e-3+.131631*x)))
      	u3=u2 
	end
;
;	dy
66:	begin
	u1=17.+exp(3.029646-x*(3.121036+x*(.086671-.216214*x))) 
      	u2=18.+exp(3.465323-x*(1.27062+x*(-.382265+x*(.431447+.303575*x))))
      	u3=u2 
	end
;
;	Ho
67:	begin
	u3=16.+exp(1.610084-x*(2.373926+x*(.133139-.071196*x))) 
      	u1=u3
      	u2=u3 
	end
;
;	Er
68:	begin
	u1=13.+exp(2.895648-x*(2.968603+x*(.561515+x*(.215267+.095813*x)))) 
      	u2=14.+exp(3.202542-x*(.852209+x*(-.226622+x*(.343738+.186042*x))))
      	u3=u2 
	end
;
;	tm
69:	begin
	u1=8.+exp(1.021172-x*(4.94757+x*(1.081603+.034811*x))) 
      	u2=9.+exp(2.173152+x*(-1.295327+x*(1.940395+.813303*x)))
      	u3=8.+exp(-.567398+x*(-3.383369+x*(.799911+.554397*x)))
	end
;
;	Yb
70:	begin
	u1=1.+exp(-2.350549-x*(6.688837+x*(1.93869+.269237*x))) 
      	u2=2.+exp(-3.047465-x*(7.390444+x*(2.355267+.44757*x)))
      	u3=1.+exp(-6.192056-x*(10.560552+x*(4.579385+.940171*x)))
	end
;
;	Lu
71:	begin
	u1=4.+exp(1.537094+x*(-1.140264+x*(.608536+.193362*x))) 
      	u2=amax1(1.,0.66+1.52e-4*t)
      	pos=where(t gt 5250.)
	if(pos(0) ne -1)then u2(pos)=-1.09+4.86e-4*t(pos)
      	u3=5.+0.*t 
	end
;
;	Hf
72:	begin
	u1=4.1758+y*(.407+y*(.57862-y*(.072887-y*3.6848e-3))) 
      	u2=-2.979+3.095e-3*t
      	u3=30. +0.*t
	end
;
;	Ta
73:	begin
	u1=3.0679+y*(.81776+y*(.34936+y*(7.4861e-3+y*3.0739e-4))) 
      	u2=1.6834+y*(2.0103+y*(.56443-y*(.031036-y*8.9565e-4)))
      	u3=15.+0.*t
	end
;
;	W
74:	begin
	u1=.3951+y*(-.25057+y*(1.4433+y*(-.34373+y*(.041924-y*1.84e-3)))) 
      	pos=where(t gt 1.2e4)
	if(pos(0) ne -1)then u1(pos)=23.+0.*pos
      	u2=1.055+y*(1.0396+y*(.3303-y*(8.4971e-3-y*5.5794e-4)))
      	u3=20.+0.*t
	end
;
;	Re
75:	begin
	u1=5.5671+y*(.72721+y*(-.42096+y*(.09075-y*3.9331e-3))) 
      	pos=where(t gt 1.2e4)
	if(pos(0) ne 0)then u1(pos)=29.+0.*pos
      	u2=6.5699+y*(.59999+y*(-.28532+y*(.050724-y*1.8544e-3)))
      	pos=where(t gt 1.2e4)
	if(pos(0) ne -1)then u2(pos)=22.+0*pos
      	u3=20.+0.*t
	end
;
;	Os
76:	begin
	u1=8.6643+y*(-.32516+y*(.68181-y*(.044252-y*1.9975e-3))) 
      	u2=9.7086+y*(-.3814+y*(.65292-y*(.064984-y*2.8792e-3)))
      	u3=10.+0.*t
	end
;
;	Ir
77:	begin
	u1=11.07+y*(-2.412+y*(1.9388+y*(-.34389+y*(.033511-1.3376e-3*y)))) 
      	pos=where(t gt 1.2e4)
	if(pos(0) ne -1)then u1(pos)=30.+0.*pos
      	u2=15.+0.*t
      	u3=20.+0.*t
	end
;
;	Pt
78:	begin
	u1=16.4+1.27e-3*t 
      	u2=6.5712+y*(-1.0363+y*(.57234-y*(.061219-2.6878e-3*y)))
      	u3=15.
	end
;
;	Au
79:	begin
	u1=1.24+2.79e-4*t 
      	u2=1.0546+y*(-.040809+y*(2.8439e-3+y*1.6586e-3))
      	u3=7.+0.*t
	end
;
;	Hg
80:	begin
	u1=1.0+0.*t 
      	u2=2.+0.*t0
      	u3=(1.+0.*t) > (.669+3.976e-5*t)
	end
;
;	Tl
81:	begin
	u1=(2.+0.*t) > (0.63+3.35e-4*t) 
      	u2=1.0+0.*t
      	u3=2.+0.*t
	end
;
;	Pb
82:	begin
	u1=(1.+0.*t) > (0.42+2.35e-4*t) 
      	pos=where(t gt 6125.)
	if(pos(0) ne -1)then u1(pos)=-1.2+5.e-4*t(pos)
      	u2=(2.+0.*t) > (1.72+7.9e-5*t)
      	u3=1.0+0.*t
	end
;
;	Bi
83:	begin
	u1=2.78+2.87e-4*t 
      	u2=(1.+0.*t) > (.37+1.41e-4*t)
      	u3=2.5 +0.*t
	end
;
;	Po
84:	begin
	u1=5.+0.*t 
      	u2=5.+0.*t
      	u3=4.+0.*t
	end
;
;	At
85:	begin
	u1=4.+0.*t  
      	u2=6.+0.*t
      	u3=6.+0.*t
	end
;
;	Rn
86:	begin
	u1=1. +0.*t
      	u2=4.+0.*t
      	u3=6.+0.*t
	end
;
;	Fr
87:	begin
	u1=2. +0.*t
      	u2=1.+0.*t
      	u3=4.5+0.*t
	end
;
;	Ra
88:	begin
	u1=1. +0.*t
      	u2=2.+0.*t
      	u3=1.+0.*t
	end
;
;	Ac
89:	begin
	u1=6. +0.*t
      	u2=3.+0.*t
      	u3=7.+0.*t
	end
;
;	Th
90:	begin
	u1=8. +0.*t
      	u2=8.+0.*t
      	u3=8.+0.*t
	end
;
;	Pa
91:	begin
	u1=50. +0.*t
      	u2=50.+0.*t
      	u3=50.+0.*t
	end
;
;	U
92:	begin
	u1=25. +0.*t
      	u2=25.+0.*t
      	u3=25.+0.*t
	end
;
else:	print,'PARTITION: no available data'
      endcase
      end
