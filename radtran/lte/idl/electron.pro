;--------------------------

	function ELECTRON,nel,t,pe

; it computes the number of electrons which a given element release to the
; medium over the number of h atoms. (it is required by gaspress)
;
; input:
;	nel	atomic number
;	t 	tempereture(k)
;	pe	electron pressure
;
; output:
;	electron electron number/total number of h atoms
;
; MODIFICATION HISTORY:
;	Sep 1, 95	Frotran > IDL (jorge)

;
; definitions
	common ATMDAT,wgt,abu,ei1,ei2,sym
;
; let's go	
	PARTITION,t,nel,u1,u2,u3
	n1overn0=saha(t,pe,u1,u2,ei1(nel-1))
	n2overn1=saha(t,pe,u2,u3,ei2(nel-1))
	n1overn=n1overn0/(1.+n1overn0*(1.+n2overn1))
	n2overn=n2overn1*n1overn	
	electron=abu(nel-1)*(n1overn+2*n2overn)
	
;
	return,electron
	end
