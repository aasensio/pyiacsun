module math_vars
implicit none
	integer, parameter :: nfac = 301
	real(kind=8) :: fact(0:nfac)
	real(kind=8) :: identity(4,4)
	
end module math_vars


module math_functions
use math_vars
use vars
use dualClass
implicit none
contains
		
!----------------------------------------------------------------
! Calculates the factorial of the integers up to 301 and save it into the fact array
!----------------------------------------------------------------
   subroutine factrl
	integer :: i
      fact(0) = 1.d0
      do i=1,301
			fact(I) = fact(I-1) * dble(I)
		enddo
   END subroutine factrl

!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
   function w3js(j1,j2,j3,m1,m2,m3)
   integer :: m1, m2, m3, j1, j2, j3
	integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
	real(kind=8) :: w3js, cc, denom, cc1, cc2

!		w3js = w3js_regge(j1/2,j2/2,j3/2,m1/2,m2/2,m3/2)
!		return
      w3js = 0.d0
      if (m1+m2+m3 /= 0) return
      ia = j1 + j2
      if (j3 > ia) return
      ib = j1 - j2
      if (j3 < abs(ib)) return
		if (abs(m1) > j1) return
		if (abs(m2) > j2) return
		if (abs(m3) > j3) return
		
      jsum = j3 + ia
      ic = j1 - m1
      id = j2 - m2
      
		ie = j3 - j2 + m1
		im = j3 - j1 - m2
		zmin = max0(0,-ie,-im)
		ig = ia - j3
		ih = j2 + m2
		zmax = min0(ig,ih,ic)
		cc = 0.d0

		do z = zmin, zmax, 2
			denom = fact(z/2)*fact((ig-z)/2)*fact((ic-z)/2)*fact((ih-z)/2)*&
				fact((ie+z)/2)*fact((im+z)/2)
      	if (mod(z,4) /= 0) denom = -denom
			cc = cc + 1.d0/denom
		enddo

		cc1 = fact(ig/2)*fact((j3+ib)/2)*fact((j3-ib)/2)/fact((jsum+2)/2)
      cc2 = fact((j1+m1)/2)*fact(ic/2)*fact(ih/2)*fact(id/2)*fact((j3-m3)/2)*fact((j3+m3)/2)
      cc = cc * sqrt(1.d0*cc1*cc2)
		if (mod(ib-m3,4) /= 0) cc = -cc
		w3js = cc
		if (abs(w3js) < 1.d-8) w3js = 0.d0		
1000 		return
   end function w3js
			
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function fvoigt_zeeman(da, dv)
	real(kind=8) :: da
	real(kind=8) :: dv(:), fvoigt_zeeman(2,size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n
		
		n = size(dv)
		do i = 1, n
			z = cmplx(dv(i), da)
			t = cmplx(da, -dv(i))
			s = dabs(dv(i)) + da
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else 
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			fvoigt_zeeman(1,i) = dble(w4)
			fvoigt_zeeman(2,i) = aimag(w4)
		enddo                        

	end function fvoigt_zeeman		

!----------------------------------------------------------------
! Voigt function
!----------------------------------------------------------------				
	function voigtDual(a, v)
		type (dual) :: a, v(:)
		type (dual) :: voigtDual(2,size(v))
		integer :: ierr, i
		real(kind=8) :: H(2,size(v)), ga(size(v)), gv(size(v)), PI = 3.14159265359d0
			
			H = fvoigt_zeeman(a%rp, v%rp)

			ga = 2.d0 * (-1.d0 / sqrt(PI) + a%rp * H(1,:) + v%rp * H(2,:))
			gv = 2.d0 * (-v%rp * H(1,:) + a%rp * H(2,:))

			voigtDual(1,:)%rp = H(1,:)
			voigtDual(2,:)%rp = H(2,:)

			do i = 1, size(v)
				voigtDual(1,i)%ip = gv(i) * v(i)%ip + ga(i) * a%ip
				voigtDual(2,i)%ip = -ga(i) * v(i)%ip + gv(i) * a%ip
			enddo
		
	end function voigtDual
	
!----------------------------------------------------------------
! This function returns the strength of the Zeeman components
! Do it using table 3.1 of Landi degl'Innocenti & Landolfi (2004)
! to avoid overflow with the factorials if using the 3-j symbol
!----------------------------------------------------------------
	function strength_zeeman(J_up,J_low,M_up,M_low)
	real(kind=8) :: J_up, J_low, M_up, M_low, strength_zeeman, strength_zeeman2

! 		strength_zeeman2 = 3.d0 * w3js(int(2.d0*J_up),int(2.d0*J_low),2,int(2.d0*M_up),&
! 					-int(2.d0*M_low),int(2.d0*(M_low-M_up)))**2

		if (J_up == J_low+1) then
			if (M_up == M_low+1) then
				strength_zeeman = (3.d0*(J_low+M_low+1)*(J_low+M_low+2)) / (2.d0*(J_low+1.d0)*(2.d0*J_low+1.d0)*(2.d0*J_low+3.d0))
 				return
			endif
			if (M_up == M_low) then
				strength_zeeman = (3.d0*(J_low-M_low+1)*(J_low+M_low+1)) / ((J_low+1.d0)*(2.d0*J_low+1.d0)*(2.d0*J_low+3.d0))
 				return
			endif
			if (M_up == M_low-1) then
				strength_zeeman = (3.d0*(J_low-M_low+1)*(J_low-M_low+2)) / (2.d0*(J_low+1.d0)*(2.d0*J_low+1.d0)*(2.d0*J_low+3.d0))
 				return
			endif
		endif

		if (J_up == J_low) then
			if (M_up == M_low+1) then
				strength_zeeman = (3.d0*(J_low-M_low)*(J_low+M_low+1)) / (2.d0*J_low*(J_low+1.d0)*(2.d0*J_low+1.d0))
 				return
			endif
			if (M_up == M_low) then
				strength_zeeman = (3.d0*M_low**2) / (J_low*(J_low+1.d0)*(2.d0*J_low+1.d0))
 				return
			endif
			if (M_up == M_low-1) then
				strength_zeeman = (3.d0*(J_low+M_low)*(J_low-M_low+1)) / (2.d0*J_low*(J_low+1.d0)*(2.d0*J_low+1.d0))
 				return
			endif
		endif

		if (J_up == J_low-1) then
			if (M_up == M_low+1) then
				strength_zeeman = (3.d0*(J_low-M_low)*(J_low-M_low-2)) / (2.d0*J_low*(2.d0*J_low-1.d0)*(2.d0*J_low+1.d0))
 				return
			endif
			if (M_up == M_low) then
				strength_zeeman = (3.d0*(J_low-M_low)*(J_low+M_low)) / (J_low*(2.d0*J_low-1.d0)*(2.d0*J_low+1.d0))
 				return
			endif
			if (M_up == M_low-1) then
				strength_zeeman = (3.d0*(J_low+M_low)*(J_low+M_low-1)) / (2.d0*J_low*(2.d0*J_low-1.d0)*(2.d0*J_low+1.d0))
 				return
			endif
		endif
					
	end function strength_zeeman

!----------------------------------------------------------------
! This function initializes some things
!----------------------------------------------------------------		
		subroutine init_maths
		integer :: i
			call factrl
      end subroutine init_maths

end module math_functions
