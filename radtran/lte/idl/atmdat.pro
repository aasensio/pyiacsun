;-----------------

      pro ATMDAT,nel,wgt,abu,ei1,ei2,symbol


; PURPOSE:
;     supplies atomic parameters for 92 elements (hydrogen to uranium)
;     adapted from the routine with the same name of
;     a.d. wittmann, goettingen (1975).(1974, solar phys. 35, 11)
;
; 	In fact the original routine has been split into two, ATMDAT,
;	which provides real atomica parameters and PARTITION,
;	which computes partition functions (dependiing on temperature
;	and pressure).
;
; INPUT:
;	nel	atomic number (this parameter can be an array!)
;			
;	
; OUTPUT:
;	wgt	atomic weight (c12=12.00)
;	abu	abundance, by number, referred to h (not in logaritms)
;		(source: from h to fe, grevesse: 1984 phys. scripta t8, 49)
;	ei1	1st ionization potential (ev)
;	ei2	2nd	"	  "        "
;		(source: allen astrophysical quantities, p.37?)
;	symbol	lowercase symbol corresponding to the elements
;
; COMMON:
;	common ATMDAT,wgt,abu,ei1,ei2,sym
;		therefore, you can also access this data set
;		by including this common block after a firt
;		call to ATMDAT
;
; MODIFICATION  HISTORY:
; started:	1975, wittmann
; finised:	22/2/90, just adding comments (j. sanchez almeida)	
;		6/3/90   change of abundances to that quoted by
;			 grevesse: 1984, phys. scripta t8, 49
;
;		August 31st, 95		translated to IDL form 
;					fortran (sos). This time atmdat
;					has been split into two routines.
;					Jorge
;
;		May 2000, symbols included

; definitions
	common ATMDAT,w,ab,e1,e2,sym
;
; checking out available data
	na=92
	if(nel lt 1 or nel gt na) then stop,'ATMDAT: no data'
;	
;	still from the FORTRAN IV version?
;c      character*2 atm
;c      dimension atm(na)/'h','he','li','be','b','c','n','o','f','ne',
;c     *'na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr',
;c     *'mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr',
;c     *'rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in',
;c     *'sn','sb','te','i','xe','cs','ba','la','ce','pr','nd','pm',
;c     *'sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w',
;c     *'re','os','ir','pt','au','hg','tl','pb','bi','po','at','rn',
;c     *'fr','ra','ac','th','pa','u'/,
;
	sym=['h','he','li','be','b','c','n','o','f','ne',$
     'na','mg','al','si','p','s','cl','ar','k','ca','sc','ti','v','cr',$
     'mn','fe','co','ni','cu','zn','ga','ge','as','se','br','kr',$
     'rb','sr','y','zr','nb','mo','tc','ru','rh','pd','ag','cd','in']
     sym=[sym,'sn','sb','te','i','xe','cs','ba','la','ce','pr','nd','pm',$
     'sm','eu','gd','tb','dy','ho','er','tm','yb','lu','hf','ta','w',$
     're','os','ir','pt','au','hg','tl','pb','bi','po','at','rn',$
     'fr','ra','ac','th','pa','u']
;
      w=[1.008,4.003,6.939,9.012,10.811,12.011,14.007,16.,18.998,$
     20.183,22.99,24.312,26.982,28.086,30.974,32.064,35.453,39.948,39.102,$
     40.08,44.956,47.90,50.942,51.996,54.938,55.847,58.933,58.71,63.54,$
     65.37,69.72,72.59,74.92,78.96,79.91,83.80,85.47,87.62,88.905,91.22,$
     92.906,95.94, 99.00,101.07,102.9,106.4,107.87,112.40,114.82]
     w=[w,118.69,121.75,127.6,126.9,131.3,132.9,137.34,138.91,140.12,140.91,$
     144.24,147.00,150.35,151.96,157.25,158.92,162.50,164.93,$
     167.26,168.93,173.04,174.97,178.49,180.95,183.85,186.2,190.2,192.2,$
     195.09,196.97,200.59,204.37,207.19,208.98,210.,211.,222.,$
     223.,226.1,227.1,232.04,231.,238.03]
;
      e1=[13.595,24.58,5.39,9.32,8.298,11.256,14.529,13.614,17.418,$
     21.559,5.138,7.644,5.984,8.149,10.474,10.357,13.012,15.755,4.339,$
     6.111,6.538,6.825,6.738,6.763,7.432,7.896,7.863,7.633,7.724,9.391,$
     5.997,7.88,9.81,9.75,11.840,13.996,4.176,5.692,6.377,6.838,6.881]
     e1=[e1,7.10,7.28,7.36,7.46,8.33,7.574,8.991,5.785,7.34,8.64,9.01,10.454,$
     12.127,3.893,5.210,5.577,5.466,5.422,5.489,5.554,5.631,5.666,6.141,$
     5.852,5.927,6.018,6.101,6.184,6.254,5.426,6.650,7.879,7.980,7.870,$
     8.70,9.10,9.00,9.22,10.43,6.105,7.415,7.287,8.43,9.30,10.745,$
     4.,5.276,6.9,6.,6.,6.]
     ;
      e2=[0.,54.403,75.62,18.21,25.15,24.376,29.59,35.11,34.98,41.07,$
     47.290,15.03,18.823,16.34,19.72,23.405,23.798,27.62,31.81,11.868,$
     12.891,13.63,14.205,16.493,15.636,16.178,17.052,18.15,20.286,17.96,$
     20.509,15.93,18.63,21.50,21.60,24.565,27.50,11.027,12.233,13.13]
     e2=[e2,14.316,16.15,15.26,16.76,18.07,19.42,21.48,16.904,18.86,14.63,$
     16.50,18.60,19.09,21.20,25.10,10.001,11.060,10.850,10.550,10.730,$
     10.899,11.069,11.241,12.090,11.519,11.670,11.800,11.930,12.050,$
     12.184,13.900,14.900,16.2,17.7,16.60,17.00,20.00,18.56,20.50,18.75,$
     20.42,15.03,16.68,19.,20.,20.,22.,10.144,12.1,12.,12.,12.]
     ;
      eps=[12.00,11.00,1.00,1.15,2.60,8.69,7.99,8.91,4.56,8.00,6.33,7.58,$
     6.47,7.55,5.45,7.21,5.50,6.58,5.12,6.36,3.10,5.02,4.0,5.67,5.45,$
     7.67,4.92,6.25,4.21,4.60,2.88,3.63,2.30,3.20,2.60,3.20,2.60,2.90]
     eps=[eps,2.10,2.75,1.95,2.16,0.00,1.83,1.20,1.50,0.85,1.85,1.65,1.90,1.10,$
     2.00,1.40,2.00,1.10,2.09,1.13,1.55,0.71,1.50,0.00,1.00,0.70,1.12,$
     .40,1.06,.50,.80,.26,1.00,.76,.80,.30,1.10,.60,1.00,1.00,1.75,.75,$
     .90,.90,1.93,.70,-8.,-8.,-8.,-8.,-8.,-8.,0.28,-8.,0.]
;
      wgt=w(nel-1)
      ei1=e1(nel-1)
      ei2=e2(nel-1)
      ab=10^(eps-12.)
      abu=ab(nel-1)
      symbol=sym(nel-1)
;
	end

