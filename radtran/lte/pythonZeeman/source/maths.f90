module mathsModule
use globalModule, only : PHC, PK, PHK, PI, PME, PH, PC, PH2C3, PHC2, EV_ERG, BOHR_R, UMA, identity
use atomicPartitionModule, only : partitionAtomic
implicit none


interface saha
	module procedure sahaVector, sahaScalar	
end interface

interface electronRelease
	module procedure electronReleaseVector, electronReleaseScalar
end interface

interface gasPressure
	module procedure gasPressureVector, gasPressureScalar
end interface

interface planckFrequency
	module procedure planckFrequencyVector, planckFrequencyScalar
end interface

interface planckWavelength
	module procedure planckWavelengthVector, planckWavelengthScalar
end interface
	
contains

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
	
! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
		subroutine splin1(x,y,yp1,ypn,y2)
		real*8, INTENT(IN) :: x(:), y(:), yp1, ypn
		real*8, INTENT(INOUT) :: y2(size(x))
		integer :: n, i, k
		real*8 :: p, qn, sig, un, u(size(x))

			n = size(x)

			if (yp1 > .99d30) then
				y2(1) = 0.d0
				u(1) = 0.d0
			else
				y2(1) = -0.5d0
				u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			endif

			do i = 2, n-1
				sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
				p = sig * y2(i-1)+2.d0
				y2(i) = (sig-1.d0)/p
				u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
					(x(i+1)-x(i-1))-sig*u(i-1))/p
			enddo
			if (ypn > .99d30) then
				qn = 0.d0
				un = 0.d0
			else
				qn = 0.5d0
				un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
			endif

			y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

			do k = n-1, 1, -1
				y2(k) = y2(k)*y2(k+1)+u(k)
			enddo

		end subroutine splin1

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
		subroutine spline(xa,ya,x,y)
		real*8, INTENT(INOUT) :: y(:)
		real*8, INTENT(IN) :: xa(:), ya(:), x(:)
		real*8 :: y2a(size(xa))
		integer :: n_x, n, i, k, khi, klo
		real*8 :: a, b, h, extrap

			n = size(xa)
			n_x = size(x)
			call splin1(xa,ya,1.d30,1.d30,y2a)

			do i = 1, n_x

! Downward extrapolation
				if (x(i) < xa(1)) then
!					y(i) = ya(1)
					y(i) = ya(1) + (ya(1)-ya(2))/(xa(1)-xa(2)) * (xa(1) - x(i))
				else

! Upward extrapolation
				if (x(i) > xa(n)) then
!					y(i) = ya(n)
					y(i) = ya(n) + (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) * (x(i) - xa(n))
				else
! In range
						klo = 1
						khi = n
1						if(khi-klo > 1) then
							k = (khi+klo)/2
							if (xa(k) > x(i)) then
								khi = k
							else
								klo = k
							endif
							go to 1
						endif

						h = xa(khi)-xa(klo)

						if (h == 0.d0) then
							! print *, 'bad xa input in spline'
							stop
						endif
						a = (xa(khi)-x(i))/h
						b = (x(i)-xa(klo))/h

						y(i) = a*ya(klo)+b*ya(khi)+((a**3.d0-a)*y2a(klo)+(b**3.d0-b)*y2a(khi))*(h**2.d0)/6.d0
					endif
				endif
			enddo

		end subroutine spline
		
!-------------------------------------------------------------------
! Calculate the natural logarithm of the Gamma function
!-------------------------------------------------------------------
	function gammln(xx)
	real(kind=8) :: gammln, xx
	integer :: j 
	real(kind=8) :: ser,tmp,x,y
	real(kind=8) :: cof(6) = (/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0,-1.231739572450155d0,&
		.1208650973866179d-2,-.5395239384953d-5/)
	real(kind=8) :: stp = 2.5066282746310005d0
		
		x = xx 
		y = x 
		tmp = x+5.5d0 
		tmp = (x+0.5d0)*log(tmp)-tmp 
		ser = 1.000000000190015d0 
		do j=1,6 
			y = y+1.d0 
			ser = ser + cof(j)/y 
		enddo 
		gammln = tmp + log(stp*ser/x)
	end function gammln

	
! ---------------------------------------------------------
! Returns Planck's function for a frequency and a set of temperatures in cgs
! INPUT:
!     - nu : frequency
!     - T : temperature
! ---------------------------------------------------------
   function planckFrequencyVector(nu,T)
	real*8, INTENT(IN) :: nu, T(:)
   real*8 :: planckFrequencyVector(size(T))
   real*8 :: temporal(size(T))

		where (T /= 0.d0)
			planckFrequencyVector = PHC * nu**3.d0 / (dexp(PHK * nu / T) - 1.d0)
		elsewhere
			planckFrequencyVector = 0.d0
		endwhere

   end function planckFrequencyVector

! ---------------------------------------------------------
! Returns Planck's function for a frequency and temperature in cgs
! INPUT:
!     - nu : frequency (Hz)
!     - T : temperature (K)
! ---------------------------------------------------------
   function planckFrequencyScalar(nu,T)
   real*8 :: planckFrequencyScalar
   real*8, INTENT(IN) :: nu, T
   real*8 :: temporal
      if (T /= 0.d0) then
         planckFrequencyScalar = PHC * nu**3.d0 / (dexp(PHK * nu / T) - 1.d0)
      else
			planckFrequencyScalar = 0.d0
      endif

   end function planckFrequencyScalar

! ---------------------------------------------------------
! Returns Planck's function for a wavelength and a set of temperatures in cgs
! INPUT:
!     - lambda : wavelength
!     - T : temperature
! ---------------------------------------------------------
   function planckWavelengthVector(lambda,T)
	real*8, INTENT(IN) :: lambda, T(:)
   real*8 :: planckWavelengthVector(size(T))
   real*8 :: temporal(size(T))

		where (T /= 0.d0)
			planckWavelengthVector = PHC2 / lambda**5.d0 / (dexp(PHK * PC / (lambda * T)) - 1.d0)
		elsewhere
			planckWavelengthVector = 0.d0
		endwhere

   end function planckWavelengthVector

! ---------------------------------------------------------
! Returns Planck's function for a wavelength and temperature in cgs
! INPUT:
!     - lambda : wavelength (cm)
!     - T : temperature (K)
! ---------------------------------------------------------
   function planckWavelengthScalar(lambda,T)
   real*8 :: planckWavelengthScalar
   real*8, INTENT(IN) :: lambda, T
   real*8 :: temporal

      if (T /= 0.d0) then
         temporal = dexp(PHK * PC / (lambda * T) - 1.d0)
         planckWavelengthScalar = PHC2 / lambda**5.d0 / temporal
      else
			planckWavelengthScalar = 0.d0
      endif

   end function planckWavelengthScalar

!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function close_to(wave_total, wave)
	real(kind=8) :: wave_total(:), wave
	real(kind=8), allocatable :: diff(:)
	integer :: n, i, location(1)
	integer :: close_to

		do i = 1, size(wave_total)
			if (wave_total(i) > wave) then
				close_to = i
				return
			endif
		enddo
! 		location = minloc(dabs(wave_total-wave))
! 		close_to = location(1)

	end function close_to

!-------------------------------------------------------------------
! Returns which is the index of the array wave_total closer to wave
!-------------------------------------------------------------------
	function extremes(lambda, leftLambda, rightLambda)
	real(kind=8) :: lambda(:), leftLambda, rightLambda
	integer :: left, right, extremes(2)


		left = close_to(lambda, leftLambda)
		right = close_to(lambda, rightLambda)

		extremes(1) = min(left,right)
		extremes(2) = max(left,right)

	end function extremes

!-----------------------------------------------------------------
! Calculates the ratio between two ionization stages using the Saha equation
!-----------------------------------------------------------------
	function sahaVector(t, pe, ul ,uu, chi)
	real(kind=8) :: t(:), pe(:), ul(:), uu(:), chi, sahaVector(size(t))
	real(kind=8) :: constant

		constant = 2.d0 * (2.d0*PI*PME)**(1.5d0) / PH**3 * PK**(2.5d0)
		sahaVector = constant * uu / ul * t**2.5d0 * 10.d0**(-5040.d0*chi/t) / pe

	end function sahaVector
	
!-----------------------------------------------------------------
! Calculates the ratio between two ionization stages using the Saha equation
!-----------------------------------------------------------------
	function sahaScalar(t, pe, ul ,uu, chi)
	real(kind=8) :: t, pe, ul, uu, chi, sahaScalar
	real(kind=8) :: constant

		constant = 2.d0 * (2.d0*PI*PME)**(1.5d0) / PH**3 * PK**(2.5d0)
		sahaScalar = constant * uu / ul * t**2.5d0 * 10.d0**(-5040.d0*chi/t) / pe

	end function sahaScalar

!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt(size(dv))
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n

		n = size(dv)
		do i = 1, n
! There is a problem for very small dampings. We force it to be 0 using a threshold
			if (da(i) < 1.d-3) da(i) = 0.d0
			z = cmplx(dv(i), da(i))
			t = cmplx(da(i), -dv(i))
			s = dabs(dv(i)) + da(i)
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
			elseif (da(i) >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
			else
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
			endif
			vecfvoigt(i) = dble(w4)
		enddo

	end function vecfvoigt
	
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt2(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt2(size(dv))
	complex :: w4(size(dv)), z(size(dv)), t(size(dv)), u(size(dv)), v4(size(dv))
	real*8 :: s(size(dv))
	integer :: i, n
                
		n = size(dv)
		
!		allocate(w4(n))
!		allocate(z(n))
!		allocate(t(n))
!		allocate(u(n))
!		allocate(v4(n))
!		allocate(s(n))
		
		where(da < 1.d-3) 
			da = 0.d0
		endwhere
		
		z = cmplx(dv,da)
		t = cmplx(da, -dv)
		s = dabs(dv) + da
		u = t*t
		
		where(s >= 15.d0)
			w4 = t * 0.5641896d0 / (0.5d0+u)
		endwhere
		where(s < 15.d0 .and. s >= 5.5d0)
			w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
		endwhere
		where(s < 15.d0 .and. s < 5.5d0 .and. da >= 0.195d0*dabs(dv)-0.176d0)
			w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
		endwhere
		where(s < 15.d0 .and. s < 5.5d0 .and. da < 0.195d0*dabs(dv)-0.176d0)
			w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
			v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
			w4 = cexp(u) - w4/v4
		endwhere
		
		vecfvoigt2 = dble(w4)
				
!		deallocate(w4)
!		deallocate(z)
!		deallocate(t)
!		deallocate(u)
!		deallocate(v4)
!		deallocate(s)

	end function vecfvoigt2
	
!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
! See Humlicek (1982) JQSRT 27, 437
!--------------------------------------------------------------
	function vecfvoigt_zeeman(da, dv)
	real*8 :: da(:)
	real*8 :: dv(:), vecfvoigt_zeeman(size(dv),2)
	complex :: w4, z, t, u, v4
	real*8 :: s
	integer :: i, n

		n = size(dv)
		do i = 1, n
! There is a problem for very small dampings. We force it to be 0 using a threshold
			if (da(i) < 1.d-3) da(i) = 0.d0
			z = cmplx(dv(i), da(i))
			t = cmplx(da(i), -dv(i))
			s = dabs(dv(i)) + da(i)
			u = t*t


			if (s >= 15.d0) then
				w4 = t * 0.5641896d0 / (0.5d0+u)
				vecfvoigt_zeeman(i,1) = dble(w4)
				vecfvoigt_zeeman(i,2) = aimag(w4)
				cycle
			elseif (s >= 5.5) then
				w4 = t*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
				vecfvoigt_zeeman(i,1) = dble(w4)
				vecfvoigt_zeeman(i,2) = aimag(w4)
				cycle
			elseif (da(i) >= 0.195d0*dabs(dv(i))-0.176d0) then
				w4 = (16.4955d0+t*(20.20933d0+t*(11.96482d0+t*(3.778987d0+t*0.5642236d0)))) / &
					(16.4955d0+t*(38.82363d0+t*(39.27121d0+t*(21.69274d0+t*(6.699398d0+t)))))
				vecfvoigt_zeeman(i,1) = dble(w4)
				vecfvoigt_zeeman(i,2) = aimag(w4)
				cycle
			else
				w4 = t*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
				v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
				w4 = cexp(u) - w4/v4
				vecfvoigt_zeeman(i,1) = dble(w4)
				vecfvoigt_zeeman(i,2) = aimag(w4)
				cycle
			endif
		enddo

	end function vecfvoigt_zeeman


!--------------------------------------------------------------
! Generates Voigt and anomalous dispersion profiles using a different damping for each frequency
! point
!*****************************************************************
! this voigt & faraday (dispersion) function is based on the     *
! paper by hui, armstrong and wray, jqsrt 19, 509 (1977). errors *
! become significant (at the < 1% level) around the knee between *
! the doppler core and the damping wings for 0.0 < a < 0.001. the*
! the normalization is such that the integral is = sqrt(pi)      *
! "a" & "vv" and output are vector of same dimension             *
! doesn't take into account case a=0 (see voigtm.f90)            *
! -------------------------------------------------------------- *
! authors: jack harvey, aake nordlund.                           *
! modified by sami solanki (1985).                               *
! modified by a.d. wittmann (1986) to include f(a,-v)            *
! transformed from function to subroutine by basilio ruiz (1993) *
! vectorized & translated to f90 by basilio ruiz (2006)          *
! -------------------------------------------------------------- *
! last update: 5-apr 06 .                                        *
!*****************************************************************
!--------------------------------------------------------------
	function voigt_faraday(a,vv)
	real(kind=8),intent(in), dimension(:):: a
	real(kind=8),intent(in), dimension(:):: vv
	real(kind=8) :: voigt_faraday(2,size(vv))
	
	complex, dimension(size(vv)) :: z

	real(kind=8) :: a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6
	data        a0,a1,a2,a3,a4,a5,a6,b0,b1,b2,b3,b4,b5,b6/         &
	 & 122.607931777104326,214.382388694706425,181.928533092181549,&
	 &  93.155580458138441, 30.180142196210589,  5.912626209773153,&
	 &    .564189583562615,122.60793177387535 ,352.730625110963558,&
	 & 457.334478783897737,348.703917719495792,170.354001821091472,&
	 &  53.992906912940207,10.479857114260399/

	real(kind=8),dimension(size(vv))    :: v,signo

	signo=1.d0
	where(vv < 0.)signo=-1.d0
	v=vv*signo

	z=cmplx(a,-v)
	z= ((((((a6*z+a5)*z+a4)*z+a3)*z+a2)*z+a1)*z+a0)/ &
	  & (((((((z+b6)*z+b5)*z+b4)*z+b3)*z+b2)*z+b1)*z+b0)

	voigt_faraday(1,:) = real(z)
	voigt_faraday(2,:) = real(.5*signo*aimag(z))
	
	end function voigt_faraday


! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary parabolically
! between points M, O and P and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)+psi(P)*S(P)
! where psi are functions of the optical distance tau(MO) and tau(OP)
! It returns the value of the three psi coefficients
! ---------------------------------------------------------
   subroutine par_sc(dtm, dtp, psim, psi0, psip)
   real(kind=8) :: short_car_parab
   real(kind=8), INTENT(IN) :: dtm, dtp
	real(kind=8), INTENT(INOUT) :: psim, psi0, psip
   real(kind=8) :: exu, u0, u1 ,u2, d2, d3, d4

         if (dtm.ge.1.d-4) then
            exu=dexp(-dtm)
            u0=1.d0-exu
            u1=dtm-1.d0+exu
            u2=dtm**2-2.d0*dtm+2.d0-2.d0*exu
         else
            d2=dtm**2
            d3=dtm**3
            d4=dtm**4
            u0=dtm-(d2/2.d0)
            u1=(d2/2.d0)-(d3/6.d0)
            u2=(d3/3.d0)-(d4/12.d0)
        endif

        if (dtm * dtp /= 0.d0 .and. dtm**2 /= 0.d0 .and. dtp**2 /= 0.d0) then
			  psim=(u2-u1*(dtp+2.d0*dtm))/(dtm*(dtm+dtp))+u0
      	  psi0=(u1*(dtm+dtp)-u2)/(dtm*dtp)
      	  psip=(u2-dtm*u1)/(dtp*(dtm+dtp))
	 	  else
		  	  psim = 0.d0
			  psi0 = 0.d0
			  psip = 0.d0
		  endif

   end subroutine par_sc

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary parabolically
! between points M, O and P and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)+psi(P)*S(P)
! where psi are functions of the optical distance tau(MO) and tau(OP)
! It returns the value of the three psi coefficients
! Works for vectorial functions
! ---------------------------------------------------------
   subroutine par_sc_vect(dtm, dtp, psim, psi0, psip)
   real(kind=8), INTENT(IN) :: dtm(:), dtp(:)
	real(kind=8), INTENT(INOUT) :: psim(size(dtm)), psi0(size(dtm)), psip(size(dtm))
   real(kind=8) :: exu(size(dtm)), u0(size(dtm)), u1(size(dtm)) ,u2(size(dtm)), d2(size(dtm))
	real(kind=8) :: d3(size(dtm)), d4(size(dtm))

		where (dtm >= 1.d-4)
			exu = dexp(-dtm)
			u0=1.d0-exu
			u1=dtm-1.d0+exu
			u2=dtm**2-2.d0*dtm+2.d0-2.d0*exu
		elsewhere
			d2=dtm**2
			d3=dtm**3
			d4=dtm**4
			u0=dtm-(d2/2.d0)
			u1=(d2/2.d0)-(d3/6.d0)
			u2=(d3/3.d0)-(d4/12.d0)
		endwhere

		psim = 0.d0
		psip = 0.d0
		psi0 = 0.d0

! 		where(dtm /= 0.d0 .or. dtp /= 0.d0)
			psim=(u2-u1*(dtp+2.d0*dtm))/(dtm*(dtm+dtp))+u0
			psi0=(u1*(dtm+dtp)-u2)/(dtm*dtp)
			psip=(u2-dtm*u1)/(dtp*(dtm+dtp))
! 		endwhere

   end subroutine par_sc_vect

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary linearly
! between points M and O and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)
! where psi are functions of the optical distance tau(MO)
! It returns the value of the two psi coefficients
! ---------------------------------------------------------
   subroutine lin_sc(dtm, psim, psi0)
   real(kind=8) :: short_car_linea
   real(kind=8), INTENT(IN) :: dtm
	real(kind=8), INTENT(INOUT) :: psim, psi0
   real(kind=8) :: exu, u0, u1, c0, cm, d2

      if (dtm.ge.1.d-4) then
         exu=dexp(-dtm)
         u0=1.d0-exu
         u1=dtm-1.d0+exu

         c0=u1/dtm
         cm=u0-c0
      else
         d2=dtm**2.d0
         c0=(dtm/2.d0)-(d2/6.d0)
         cm=(dtm/2.d0)-(d2/3.d0)
      endif
		psi0 = c0
		psim = cm

   end subroutine lin_sc

! ---------------------------------------------------------
! Solves RTE using short-characteristics. Source function is supposed to vary linearly
! between points M and O and RTE is integrated analitically, resulting in
!            I(O)=I(M)*exp(-tau(MO))+psi(M)*S(M)+psi(O)*S(O)
! where psi are functions of the optical distance tau(MO)
! It returns the value of the two psi coefficients
! Works for vectorial functions
! ---------------------------------------------------------
   subroutine lin_sc_vect(dtm, psim, psi0)
   real(kind=8), INTENT(IN) :: dtm(:)
	real(kind=8), INTENT(INOUT) :: psim(size(dtm)), psi0(size(dtm))
   real(kind=8) :: exu(size(dtm)), u0(size(dtm)), u1(size(dtm)), c0(size(dtm)), cm(size(dtm)), d2(size(dtm))

	where (dtm >= 1.d-4)
		exu=dexp(-dtm)
		u0=1.d0-exu
		u1=dtm-1.d0+exu
		c0=u1/dtm
		cm=u0-c0
	elsewhere
		d2=dtm**2.d0
		c0=(dtm/2.d0)-(d2/6.d0)
		cm=(dtm/2.d0)-(d2/3.d0)
	endwhere

	psi0 = c0
	psim = cm

   end subroutine lin_sc_vect

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel atmosphere
! for a given source function and opacity
! Opacity and source function have the index as
! opacity(depth,frequency)
! ---------------------------------------------------------
	function shortCharacteristics(height, opacity, source, mu, boundary, idir)
	real(kind=8) :: height(:), opacity(:,:), source(:,:), boundary(:), mu
	real(kind=8) :: shortCharacteristics(size(opacity(1,:)))
	integer :: k, km, kp, idir, k0, kf, ndepths
	real(kind=8), dimension(size(opacity(1,:))) :: chim, chi0, chip, sm, s0, sp, Inten, dtp, dtm, exu
	real(kind=8), dimension(size(opacity(1,:))) :: psim, psi0, psip
	real(kind=8) :: dm, dp

		ndepths = size(height)
		
! Boundary condition
		if (idir == 1) then
			k0 = 2
			kf = ndepths
			Inten = boundary
		endif
		if (idir == -1) then
			k0 = ndepths-1
			kf = 1
			Inten = boundary
		endif

		do k = k0, kf, idir

! Parabolic short-characteristics
			if (k /= kf) then
				km = k - idir
				kp = k + idir
				chim = opacity(km,:)
				chi0 = opacity(k,:)
				chip = opacity(kp,:)
				sm = source(km,:)
				s0 = source(k,:)
				sp = source(kp,:)
				dm = dabs(1.d0/mu*(height(k) - height(km)))
				dp = dabs(1.d0/mu*(height(kp) - height(k)))
			else
! Linear short-characteristics
				km = k - idir
				chim = opacity(km,:)
				chi0 = opacity(k,:)
				chip = 0.d0
				sm = source(km,:)
				s0 = source(k,:)
				sp = 0.d0
				dm = dabs(1.d0/mu*(height(k) - height(km)))
				dp = 0.d0
			endif

			where(chim <= 0.d0)
				chim = 1.d-15
			endwhere
			where(chi0 <= 0.d0)
				chi0 = 1.d-15
			endwhere
			where(chip <= 0.d0)
				chip = 1.d-15
			endwhere
			dtm = 0.5d0 * (chim + chi0) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

			where (dtm >= 1.d-4)
				exu = dexp(-dtm)
			elsewhere
				exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
			endwhere

			if (k /= kf) then
 				call par_sc_vect(dtm,dtp,psim,psi0,psip)
				Inten = Inten * exu + psim*sm + psi0*s0 + psip*sp
			else
				call lin_sc_vect(dtm,psim,psi0)
				Inten = Inten * exu + psim*sm + psi0*s0
			endif

		enddo

		shortCharacteristics = Inten

	end function shortCharacteristics

	
!---------------------------------------------------------

!                                                         robert renka
!                                                 oak ridge natl. lab.

!   this subroutine uses an order n*log(n) quick sort to sort a real (dp)
! array x into increasing order.  the algorithm is as follows.  ind is
! initialized to the ordered sequence of indices 1,...,n, and all interchanges
! are applied to ind.  x is divided into two portions by picking a central
! element t.  the first and last elements are compared with t, and
! interchanges are applied as necessary so that the three values are in
! ascending order.  interchanges are then applied so that all elements
! greater than t are in the upper portion of the array and all elements
! less than t are in the lower portion.  the upper and lower indices of one
! of the portions are saved in local arrays, and the process is repeated
! iteratively on the other portion.  when a portion is completely sorted,
! the process begins again by retrieving the indices bounding another
! unsorted portion.

! input parameters -   n - length of the array x.

!                      x - vector of length n to be sorted.

!                    ind - vector of length >= n.

! n and x are not altered by this routine.

! output parameter - ind - sequence of indices 1,...,n permuted in the same
!                          fashion as x would be.  thus, the ordering on
!                          x is defined by y(i) = x(ind(i)).

!---------------------------------------------------------

! note -- iu and il must be dimensioned >= log(n) where log has base 2.

!*********************************************************************

	subroutine qsortd(x,ind,n)
	
	implicit none
	integer, parameter  :: dp = 4
	
	real (kind=8), intent(in)  :: x(:)
	integer, intent(out)   :: ind(:)
	integer, intent(in)    :: n
	
	
	integer   :: iu(21), il(21)
	integer   :: m, i, j, k, l, ij, it, itt, indx
	real(kind=8)      :: r
	real(kind=8) :: t
	
	! local parameters -
	
	! iu,il =  temporary storage for the upper and lower
	!            indices of portions of the array x
	! m =      index for iu and il
	! i,j =    lower and upper indices of a portion of x
	! k,l =    indices in the range i,...,j
	! ij =     randomly chosen index between i and j
	! it,itt = temporary storage for interchanges in ind
	! indx =   temporary index for x
	! r =      pseudo random number for generating ij
	! t =      central element of x
	
	if (n <= 0) return
	
	! initialize ind, m, i, j, and r
	
	do  i = 1, n
	ind(i) = i
	end do
	
	m = 1
	i = 1
	j = n
	r = .375
	
	! top of loop
	
	20 if (i >= j) go to 70
	if (r <= .5898437) then
	r = r + .0390625
	else
	r = r - .21875
	end if
	
	! initialize k
	
	30 k = i
	
	! select a central element of x and save it in t
	
	ij = i + r*(j-i)
	it = ind(ij)
	t = x(it)
	
	! if the first element of the array is greater than t,
	!   interchange it with t
	
	indx = ind(i)
	if (x(indx) > t) then
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	end if
	
	! initialize l
	
	l = j
	
	! if the last element of the array is less than t,
	!   interchange it with t
	
	indx = ind(j)
	if (x(indx) >= t) go to 50
	ind(ij) = indx
	ind(j) = it
	it = indx
	t = x(it)
	
	! if the first element of the array is greater than t,
	!   interchange it with t
	
	indx = ind(i)
	if (x(indx) <= t) go to 50
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	go to 50
	
	! interchange elements k and l
	
	40 itt = ind(l)
	ind(l) = ind(k)
	ind(k) = itt
	
	! find an element in the upper part of the array which is
	!   not larger than t
	
	50 l = l - 1
	indx = ind(l)
	if (x(indx) > t) go to 50
	
	! find an element in the lower part of the array whcih is not smaller than t
	
	60 k = k + 1
	indx = ind(k)
	if (x(indx) < t) go to 60
	
	! if k <= l, interchange elements k and l
	
	if (k <= l) go to 40
	
	! save the upper and lower subscripts of the portion of the
	!   array yet to be sorted
	
	if (l-i > j-k) then
	il(m) = i
	iu(m) = l
	i = k
	m = m + 1
	go to 80
	end if
	
	il(m) = k
	iu(m) = j
	j = l
	m = m + 1
	go to 80
	
	! begin again on another unsorted portion of the array
	
	70 m = m - 1
	if (m == 0) return
	i = il(m)
	j = iu(m)
	
	80 if (j-i >= 11) go to 30
	if (i == 1) go to 20
	i = i - 1
	
	! sort elements i+1,...,j.  note that 1 <= i < j and j-i < 11.
	
	90 i = i + 1
	if (i == j) go to 70
	indx = ind(i+1)
	t = x(indx)
	it = indx
	indx = ind(i)
	if (x(indx) <= t) go to 90
	k = i
	
	100 ind(k+1) = ind(k)
	k = k - 1
	indx = ind(k)
	if (t < x(indx)) go to 100
	
	ind(k+1) = it
	go to 90
	end subroutine qsortd
	
! ---------------------------------------------------------
! Simple routine to compute the derivative dy/dx
! ---------------------------------------------------------
	subroutine deriv(x,y,dydx)
	real(kind=8) :: x(:),y(:),dydx(size(x))
	real(kind=8) :: r1,r2,r3,r4,s1,s2
	real(kind=8) :: A,B,C
	integer(kind=4) ::  i,n

    do i=2,size(x)-1
       s1 = (x(i-1)-x(i))*(x(i)**2-x(i+1)**2)
       s2 = (x(i)-x(i+1))*(x(i-1)**2-x(i)**2)
       r1 = (y(i-1)-y(i))*(x(i)-x(i+1))
       r2 = (y(i)-y(i+1))*(x(i-1)-x(i))
       r3 = (y(i-1)-y(i))*(x(i)**2-x(i+1)**2)
       r4 = (y(i)-y(i+1))*(x(i-1)**2-x(i)**2)
       A = (r1-r2)/(s2-s1)
       B = (r3-r4)/(s1-s2)
       C = y(i) - A*x(i)**2 - B*x(i)
       dydx(i) = 2*A*x(i) + B
       if(i.eq.2) dydx(i-1) = 2*A*x(i-1) + B
       if(i.eq.(size(x)-1)) dydx(i+1) = 2*A*x(i+1) + B
    enddo

    end subroutine deriv
    
!-----------------------------------------------------------------
! Calculate the damping constant a for the Voigt profile
! following Gray 1976 "The observation and analysis of stellar photospheres"
!  Natural damping: eq (11.19)
!  van der Waals damping: eq (11.35) & (11.36)
!  Stark damping: eq (11.35) & (11.36)
!  Damping: eq (11.69)
! ion_pot : ionization potential in eV
! exc_pot_lower : excitation potential of the lower level in eV
! T : temperature array
! Pg : gas pressure array
! lambda : wavelength of the transition in A
! doppler_width : Doppler width in cm/s
!-----------------------------------------------------------------
	function calculateDamping(T, nHtot, Pe, Pg, doppler_width, ion_pot, exc_pot_lower, gamma_rad, gamma_stark, gamma_vdw, lambda, alpha, sigma0, speciesMass)
	real(kind=8) :: ion_pot, exc_pot_lower, T(:), nHtot(:), Pe(:), Pg(:), lambda, doppler_width(:), speciesMass
	real(kind=8) :: calculateDamping(size(T)), gamma_rad, gamma_stark, gamma_vdw, alpha, sigma0, v0, m1, m2, reduced_mass, t1, t2, t3, t4
	integer :: n, trans
	real(kind=8) :: c6, chi_lambda, cte, a, b, alphaABO, sigmaABO
	real(kind=8), allocatable :: gamma6(:), gamma_nat(:), gamma_s(:)

		n = size(T)
		
		if (allocated(gamma6)) deallocate(gamma6)
		if (allocated(gamma_s)) deallocate(gamma_s)
		if (allocated(gamma_nat)) deallocate(gamma_nat)
		
		allocate(gamma6(n))
		allocate(gamma_s(n))
		allocate(gamma_nat(n))


! If the damping is not in the linelist, then use the formulas. In the other case, use those
! given by Kurucz

! RADIATIVE DAMPING
! Radiative damping
! 			Reference: Gray 1976, "The observation and analysis of stellar
!  		photospheres", pag 227, just after Eq. (11-19). This approximation
!  		is poor, but the radiative damping is usually negligible anyway.		

		gamma_nat = 0.22233d0 / (lambda*1.d-8)**2	

! STARK DAMPING
! Stark damping
! Formula used: gamma_4=38.8*(v**(1./3.))*(C4**(2./3.))*N , from
!   Unsold 1955, "Physik der Sternatmospharen", 2nd ed.,
!   Springer-Verlag, pp. 326 and 331. According to Gray (ref above),
!   pag 237, this is similar to
!   log (gamma_4) = 19.4 + 2./3*log C4 + log Pe - 5./6.*log T
!   The value of log C4 is approximated by -14. (see Gray, Table 11-3)

		gamma_s = 19.4 + 2.d0/3.d0*(-14.d0) + dlog10(Pe) - 5.d0/6.d0*dlog10(T)
		gamma_s = 10.d0**(gamma_s)

! HYDROGEN DAMPING
! Van der Waals damping:
! Reference: Gray (see above, pag 239, Eqs. (11-35) and (11-36))
! Formula used:
!    log (gamma_vdw) = 19.6 + 2./5.*log (C6(H)) + log Pg - 7./10.*log T
		
		if (alpha == 0 .and. sigma0 == 0) then
			chi_lambda = (PH*PC/EV_ERG) / (lambda*1.d-8)
			a = ion_pot-exc_pot_lower-chi_lambda
			b = ion_pot-exc_pot_lower
			c6 = 0.3d-30 * (1.d0/a**2 - 1.d0/b**2)
	! It seems that the ionization potential and the energy of some levels are inconsistent, and E_exc>I
			if (c6 > 0.d0) then
				gamma6 = 10.d0**(19.6d0 + 0.4d0*dlog10(c6) + dlog10(Pg) - 0.7*dlog10(T))
			else
				gamma6 = 0.d0
			endif
		else
			alphaABO = alpha
			sigmaABO = sigma0 * BOHR_R**2
			v0 = 1.d6				! 10^4 m s^-1
			m1 = speciesMass				! Species
			m2 = 1.008d0				! Hydrogen
			reduced_mass = (m1*m2) / (m1+m2) * UMA
								
			t1 = (alpha/2.d0)*dlog(4.d0/PI)
			t2 = gammln((4.d0-alpha)/2.d0)
			t3 = alpha*dlog(v0)
			t4 = dlog(sigmaABO)
			
			gamma6 = t1 + t2 + t3 + t4 + 0.5d0*(1.d0-alpha)*dlog((8.d0*PK*T)/(PI*reduced_mass))
			gamma6 = exp(gamma6)
						
! The 2 is because the formula gives the half half-width
! We only include the collisions with hydrogen, although in principle it could be possible to include all
			gamma6 = 2.d0 * gamma6 * nHtot
			
		endif

		cte = 1.d-8 / (4.d0 * PI * PC)
		calculateDamping = cte * (gamma6+gamma_nat+gamma_s) * lambda**2 / (doppler_width * lambda / PC)

		deallocate(gamma6)
		deallocate(gamma_nat)
		deallocate(gamma_s)		

	end function calculateDamping

!-----------------------------------------------------------------
! Compute the height axis from the opacity at 5000 A and the optical depth at the same wavelength
!-----------------------------------------------------------------
	subroutine computeHeight(lTau5, opacity5, height)
	real(kind=8) :: lTau5(:), opacity5(:), height(:), avgOpacity
	integer :: i, loc(1)
	
		height = 0.d0
 		do i = 2, size(lTau5)
 			avgOpacity = 0.5d0 * (opacity5(i) + opacity5(i-1))
 			height(i) = height(i-1) - (10.d0**lTau5(i) - 10.d0**lTau5(i-1)) / avgOpacity
 		enddo
				
! Put z=0 at the point where tau500=1
		loc = minloc(abs(10.d0**lTau5-1.d0))		
		height = height - height(loc(1))
		
	end subroutine computeHeight

!-----------------------------------------------------------------
! Compute the partial pressures of H, H-, H+, H2, H2+ as well as the total pressure
!-----------------------------------------------------------------
	function electronReleaseVector(nel, T, Pe)
	integer :: nel
	real(kind=8) :: T(:), Pe(:), electronReleaseVector(size(T)), u1(size(T)), u2(size(T)), u3(size(T)), ei1, ei2, weight, abundance
	real(kind=8), dimension(size(T)) :: n1overn0, n2overn1, n1overn, n2overn
	
		call partitionAtomic(T, nel, u1, u2, u3, ei1, ei2, weight, abundance)
		n1overn0 = saha(T, Pe, u1, u2, ei1)
		n2overn1 = saha(T, Pe, u2, u3, ei2)
		n1overn = n1overn0 / (1.d0 + n1overn0 * (1.d0 + n2overn1))
		n2overn = n2overn1 * n1overn	
		electronReleaseVector = abundance * (n1overn + 2 * n2overn)
		return
		
	end function electronReleaseVector
	
!-----------------------------------------------------------------
! Compute the partial pressures of H, H-, H+, H2, H2+ as well as the total pressure
!-----------------------------------------------------------------
	function electronReleaseScalar(nel, T, Pe)
	integer :: nel
	real(kind=8) :: T, Pe, electronReleaseScalar, u1, u2, u3, ei1, ei2, weight, abundance
	real(kind=8) :: n1overn0, n2overn1, n1overn, n2overn
	
		call partitionAtomic(T, nel, u1, u2, u3, ei1, ei2, weight, abundance)
		n1overn0 = saha(T, Pe, u1, u2, ei1)
		n2overn1 = saha(T, Pe, u2, u3, ei2)
		n1overn = n1overn0 / (1.d0 + n1overn0 * (1.d0 + n2overn1))
		n2overn = n2overn1 * n1overn	
		electronReleaseScalar = abundance * (n1overn + 2 * n2overn)
		return
		
	end function electronReleaseScalar

!-----------------------------------------------------------------
! Compute the partial pressures of H, H-, H+, H2, H2+ as well as the total pressure
!-----------------------------------------------------------------
	subroutine gasPressureVector(Pe, T, PH, PHminus, PHplus, PH2, PH2plus, Pg)
	real(kind=8), dimension(:) :: Pe, T, PH, PHminus, PHplus, PH2, PH2plus, Pg
	real(kind=8), parameter, dimension(4) :: ch2 = (/12.533505d0,-4.9251644d0,5.6191273d-2,-3.2687661d-3/)
	real(kind=8), parameter, dimension(4) :: ch2pl = (/1.1206998d1,-2.7942767d0,-7.9196803d-2,2.4790744d-2/)
	real(kind=8), parameter, dimension(3) :: chpl = (/-13.595d0,2.5d0,-0.4772d0/)
	real(kind=8), parameter, dimension(3) :: chmi = (/-0.747d0,2.5d0,.1249d0/)
!	Total abundance by number (it is not that of atmdat but comes from that adopted in the hsra model atmosphere; gingerich et al.: 1971, solar phys. 18, 347)
	real(kind=8), parameter :: totAbundance = 0.101d0
	real(kind=8), dimension(size(T)) :: theta, coh2, coh2pl, cohpl, cohmi, g1, g2, g3, g4, g5, a, b, c, d, e, c1, c2, c3, f1, f2, f3, f4, f5, f6
	
!
! Equilibrium constants
		theta = 5040.d0 / T

		coh2 = ch2(1)+(ch2(2)+(ch2(3)+ch2(4)*theta)*theta)*theta	
		coh2=10.d0**coh2

		coh2pl = ch2pl(1)+(ch2pl(2)+(ch2pl(3)+ch2pl(4)*theta)*theta)*theta
		coh2pl = 10.d0**coh2pl
		
!	just the saha equation
		cohpl = chpl(1)*theta+chpl(2)*log10(t)+chpl(3)
		cohpl = 10.d0**cohpl
!	just the saha equation
		cohmi = chmi(1)*theta+chmi(2)*log10(t)+chmi(3)
		cohmi = 10.d0**cohmi

! now it solves a system of six equations with six unknows.
! 	contribution to pe from elements which are not h (two ionizations)
! 	(the relevant elemets have been taken from mihalas(1967))
!	Warning, I've commented out those elements which 
!	are not used for the HSRA (table II, SPh 18, 357)

		g1 = electronRelease(2,T,Pe)  ! he
		g1 = g1 + electronRelease(6,t,pe)	! c
		g1 = g1 + electronRelease(7,t,pe)	! n
		g1 = g1 + electronRelease(8,t,pe)	! o
		g1 = g1 + electronRelease(11,t,pe) ! na
		g1 = g1 + electronRelease(12,t,pe) ! mg
		g1 = g1 + electronRelease(13,t,pe) ! al
		g1 = g1 + electronRelease(14,t,pe) ! si
		g1 = g1 + electronRelease(16,t,pe) ! s
		g1 = g1 + electronRelease(19,t,pe) ! k
		g1 = g1 + electronRelease(20,t,pe) ! ca
		g1 = g1 + electronRelease(24,t,pe) ! cr
		g1 = g1 + electronRelease(26,t,pe) ! fe
	
		g2 = cohpl / pe
		g3 = pe/cohmi
		g4 = pe/coh2pl
		g5 = pe/coh2

		a = 1+g2+g3
		b = 2*(1+g2*g4/g5)
		c = g5
		d = g2-g3
		e = g2*g4/g5

		c1 = c*b*b+a*d*b-e*a*a
		c2 = 2*a*e-d*b+a*b*g1
		c3 = (-1)*(e+b*g1)
		
		f1 = (-1)*c2/2/c1
		f2 = sqrt((c2/2/c1)**2.d0-c3/c1)
		
		where(c1 >= 0.d0)
			f1 = f1 + f2
		elsewhere		
			f1 = f1 - f2
		endwhere
				
		f2 = g2*f1
		f3 = g3*f1
		f5 = (1-a*f1)/b
		f4 = e*f5
		f6 = pe/(f2-f3+g1+f4)

		PH = f1*f6
		PHplus = f2*f6
		PHminus = f3*f6
		PH2plus = f4*f6
		PH2 = f5*f6
		pg = pe + f6*(totAbundance+(f1+f2+f3+f4+f5))	
		
	end subroutine gasPressureVector
	
!-----------------------------------------------------------------
! Compute the partial pressures of H, H-, H+, H2, H2+ as well as the total pressure
!-----------------------------------------------------------------
	subroutine gasPressureScalar(Pe, T, PH, PHminus, PHplus, PH2, PH2plus, Pg)
	real(kind=8) :: Pe, T, PH, PHminus, PHplus, PH2, PH2plus, Pg
	real(kind=8), parameter, dimension(4) :: ch2 = (/12.533505d0,-4.9251644d0,5.6191273d-2,-3.2687661d-3/)
	real(kind=8), parameter, dimension(4) :: ch2pl = (/1.1206998d1,-2.7942767d0,-7.9196803d-2,2.4790744d-2/)
	real(kind=8), parameter, dimension(3) :: chpl = (/-13.595d0,2.5d0,-0.4772d0/)
	real(kind=8), parameter, dimension(3) :: chmi = (/-0.747d0,2.5d0,.1249d0/)
!	Total abundance by number (it is not that of atmdat but comes from that adopted in the hsra model atmosphere; gingerich et al.: 1971, solar phys. 18, 347)
	real(kind=8), parameter :: totAbundance = 0.101d0
	real(kind=8) :: theta, coh2, coh2pl, cohpl, cohmi, g1, g2, g3, g4, g5, a, b, c, d, e, c1, c2, c3, f1, f2, f3, f4, f5, f6
	
!
! Equilibrium constants
		theta = 5040.d0 / T

		coh2 = ch2(1)+(ch2(2)+(ch2(3)+ch2(4)*theta)*theta)*theta	
		coh2=10.d0**coh2

		coh2pl = ch2pl(1)+(ch2pl(2)+(ch2pl(3)+ch2pl(4)*theta)*theta)*theta
		coh2pl = 10.d0**coh2pl
		
!	just the saha equation
		cohpl = chpl(1)*theta+chpl(2)*log10(t)+chpl(3)
		cohpl = 10.d0**cohpl
!	just the saha equation
		cohmi = chmi(1)*theta+chmi(2)*log10(t)+chmi(3)
		cohmi = 10.d0**cohmi

! now it solves a system of six equations with six unknows.
! 	contribution to pe from elements which are not h (two ionizations)
! 	(the relevant elemets have been taken from mihalas(1967))
!	Warning, I've commented out those elements which 
!	are not used for the HSRA (table II, SPh 18, 357)

		g1 = electronRelease(2,T,Pe)  ! he
		g1 = g1 + electronRelease(6,t,pe)	! c
		g1 = g1 + electronRelease(7,t,pe)	! n
		g1 = g1 + electronRelease(8,t,pe)	! o
		g1 = g1 + electronRelease(11,t,pe) ! na
		g1 = g1 + electronRelease(12,t,pe) ! mg
		g1 = g1 + electronRelease(13,t,pe) ! al
		g1 = g1 + electronRelease(14,t,pe) ! si
		g1 = g1 + electronRelease(16,t,pe) ! s
		g1 = g1 + electronRelease(19,t,pe) ! k
		g1 = g1 + electronRelease(20,t,pe) ! ca
		g1 = g1 + electronRelease(24,t,pe) ! cr
		g1 = g1 + electronRelease(26,t,pe) ! fe
	
		g2 = cohpl / pe
		g3 = pe/cohmi
		g4 = pe/coh2pl
		g5 = pe/coh2

		a = 1+g2+g3
		b = 2*(1+g2*g4/g5)
		c = g5
		d = g2-g3
		e = g2*g4/g5

		c1 = c*b*b+a*d*b-e*a*a
		c2 = 2*a*e-d*b+a*b*g1
		c3 = (-1)*(e+b*g1)
		
		f1 = (-1)*c2/2/c1
		f2 = sqrt((c2/2/c1)**2.d0-c3/c1)
		
		if (c1 >= 0.d0) then
			f1 = f1 + f2
		else
			f1 = f1 - f2
		endif
				
		f2 = g2*f1
		f3 = g3*f1
		f5 = (1-a*f1)/b
		f4 = e*f5
		f6 = pe/(f2-f3+g1+f4)

		PH = f1*f6
		PHplus = f2*f6
		PHminus = f3*f6
		PH2plus = f4*f6
		PH2 = f5*f6
		pg = pe + f6*(totAbundance+(f1+f2+f3+f4+f5))	
		
	end subroutine gasPressureScalar
	

! ---------------------------------------------------------
! This routine performs an interpolation of the curve that
! passes through x,y to xx,yy. The interpolation is performed joining
! adjacent points by a straight line and then smoothing the result 
! by convolving it with a Gaussian whose width is a half of the distance
! between the two closest nodes.
! ---------------------------------------------------------
	subroutine smoothed_lines(n, nn, x, y, xx, yy)
	integer :: n, nn, inode, ipoint, jpoint, i
	real(kind=8), dimension (n) :: x, y
	real(kind=8), dimension (nn) :: xx, yy, yy2
	real(kind=8), dimension (n) :: xdouble, ydouble
	real(kind=8), dimension (nn) :: xxdouble, yydouble
	real :: x1, x2, y1, y2, closest, expo, norm

		xdouble = x
		ydouble = y
		xxdouble = xx
		call spline(xdouble,ydouble,xxdouble,yydouble)
		yy = real(yydouble)

		return

! compute segments
	yy=0.
	do inode=1, n-1
		x1=x(inode)
		x2=x(inode+1)
		y1=y(inode)
		y2=y(inode+1)
		do ipoint=1, nn
			if (xx(ipoint) .ge. x1 .and. xx(ipoint) .le. x2) then
				yy(ipoint)=y1+(y2-y1)/(x2-x1)*(xx(ipoint)-x1)
			endif
		enddo
	enddo

! Smooth the result
! Find the closest nodes

	closest=abs(x(1)-x(2))
	do inode=2, n-1
		if (abs(x(inode)-x(inode+1)) .lt. closest) &
				closest=abs(x(inode)-x(inode+1))
	end do

! convolve (note that the model is not necessarily equispaced)
	closest=closest/2.
	yy2=0.
	do ipoint=1, nn
		norm=0.
		do jpoint=1,nn
			x1=xx(jpoint)-xx(ipoint)
			if (x1*x1/closest/closest .lt. 25.) then
				expo=exp(-x1*x1/(closest*closest))
				yy2(ipoint)=yy2(ipoint)+yy(jpoint)*expo
				norm=norm+expo
			end if
		end do
		yy2(ipoint)=yy2(ipoint)/norm
	end do
	yy=yy2
		
	end subroutine smoothed_lines

!--------------------------------------------------------------
! Inversion of a 4x4 matrix
!--------------------------------------------------------------
	subroutine invert(a)
	real(kind=8) :: a(4,4)
	real(kind=8) :: b(4,4), det, maxim, fabsmax
	! First some tests of singularity
		b = dabs(a)
		maxim = maxval(b)
		fabsmax = 1.d0 / maxim
		if (maxim == 0.d0) then
			! print *, 'Singularity in the inversion'
			stop
		endif

		a = a * fabsmax

   	b(1,1) = a(2,2) * a(3,3) * a(4,4) + a(2,3) * a(3,4) * a(4,2)&
      	+ a(2,4) * a(3,2) * a(4,3) - a(2,2) * a(3,4) * a(4,3)&
	   	- a(2,3) * a(3,2) * a(4,4) - a(2,4) * a(3,3) * a(4,2)
   	b(2,1) = a(2,3) * a(3,1) * a(4,4) + a(2,4) * a(3,3) * a(4,1)&
	   	+ a(2,1) * a(3,4) * a(4,3) - a(2,3) * a(3,4) * a(4,1)&
   		- a(2,4) * a(3,1) * a(4,3) - a(2,1) * a(3,3) * a(4,4)
   	b(3,1) = a(2,4) * a(3,1) * a(4,2) + a(2,1) * a(3,2) * a(4,4)&
	   	+ a(2,2) * a(3,4) * a(4,1) - a(2,4) * a(3,2) * a(4,1)&
   		- a(2,1) * a(3,4) * a(4,2) - a(2,2) * a(3,1) * a(4,4)
   	b(4,1) = a(2,1) * a(3,3) * a(4,2) + a(2,2) * a(3,1) * a(4,3)&
	   	+ a(2,3) * a(3,2) * a(4,1) - a(2,1) * a(3,2) * a(4,3)&
   		- a(2,2) * a(3,3) * a(4,1) - a(2,3) * a(3,1) * a(4,2)
   	b(1,2) = a(3,2) * a(4,4) * a(1,3) + a(3,3) * a(4,2) * a(1,4)&
	   	+ a(3,4) * a(4,3) * a(1,2) - a(3,2) * a(4,3) * a(1,4)&
   		- a(3,3) * a(4,4) * a(1,2) - a(3,4) * a(4,2) * a(1,3)
   	b(2,2) = a(3,3) * a(4,4) * a(1,1) + a(3,4) * a(4,1) * a(1,3)&
	   	+ a(3,1) * a(4,3) * a(1,4) - a(3,3) * a(4,1) * a(1,4)&
   		- a(3,4) * a(4,3) * a(1,1) - a(3,1) * a(4,4) * a(1,3)
   	b(3,2) = a(3,4) * a(4,2) * a(1,1) + a(3,1) * a(4,4) * a(1,2)&
	   	+ a(3,2) * a(4,1) * a(1,4) - a(3,4) * a(4,1) * a(1,2)&
   		- a(3,1) * a(4,2) * a(1,4) - a(3,2) * a(4,4) * a(1,1)
   	b(4,2) = a(3,1) * a(4,2) * a(1,3) + a(3,2) * a(4,3) * a(1,1)&
	   	+ a(3,3) * a(4,1) * a(1,2) - a(3,1) * a(4,3) * a(1,2)&
   		- a(3,2) * a(4,1) * a(1,3) - a(3,3) * a(4,2) * a(1,1)
   	b(1,3) = a(4,2) * a(1,3) * a(2,4) + a(4,3) * a(1,4) * a(2,2)&
   		+ a(4,4) * a(1,2) * a(2,3) - a(4,2) * a(1,4) * a(2,3)&
	   	- a(4,3) * a(1,2) * a(2,4) - a(4,4) * a(1,3) * a(2,2)
   	b(2,3) = a(4,3) * a(1,1) * a(2,4) + a(4,4) * a(1,3) * a(2,1)&
	   	+ a(4,1) * a(1,4) * a(2,3) - a(4,3) * a(1,4) * a(2,1)&
   		- a(4,4) * a(1,1) * a(2,3) - a(4,1) * a(1,3) * a(2,4)
   	b(3,3) = a(4,4) * a(1,1) * a(2,2) + a(4,1) * a(1,2) * a(2,4)&
	   	+ a(4,2) * a(1,4) * a(2,1) - a(4,4) * a(1,2) * a(2,1)&
   		- a(4,1) * a(1,4) * a(2,2) - a(4,2) * a(1,1) * a(2,4)
   	b(4,3) = a(4,1) * a(1,3) * a(2,2) + a(4,2) * a(1,1) * a(2,3)&
	   	+ a(4,3) * a(1,2) * a(2,1) - a(4,1) * a(1,2) * a(2,3)&
   		- a(4,2) * a(1,3) * a(2,1) - a(4,3) * a(1,1) * a(2,2)
   	b(1,4) = a(1,2) * a(2,4) * a(3,3) + a(1,3) * a(2,2) * a(3,4)&
	   	+ a(1,4) * a(2,3) * a(3,2) - a(1,2) * a(2,3) * a(3,4)&
   		- a(1,3) * a(2,4) * a(3,2) - a(1,4) * a(2,2) * a(3,3)
   	b(2,4) = a(1,3) * a(2,4) * a(3,1) + a(1,4) * a(2,1) * a(3,3)&
	   	+ a(1,1) * a(2,3) * a(3,4) - a(1,3) * a(2,1) * a(3,4)&
   		- a(1,4) * a(2,3) * a(3,1) - a(1,1) * a(2,4) * a(3,3)
   	b(3,4) = a(1,4) * a(2,2) * a(3,1) + a(1,1) * a(2,4) * a(3,2)&
	   	+ a(1,2) * a(2,1) * a(3,4) - a(1,4) * a(2,1) * a(3,2)&
   		- a(1,1) * a(2,2) * a(3,4) - a(1,2) * a(2,4) * a(3,1)
   	b(4,4) = a(1,1) * a(2,2) * a(3,3) + a(1,2) * a(2,3) * a(3,1)&
	   	+ a(1,3) * a(2,1) * a(3,2) - a(1,1) * a(2,3) * a(3,2)&
   		- a(1,2) * a(2,1) * a(3,3) - a(1,3) * a(2,2) * a(3,1)

		det = a(1,1) * b(1,1) + a(1,2) * b(2,1) + a(1,3) * b(3,1) + a(1,4) * b(4,1)

		a = b * (fabsmax / det)	

	end subroutine invert

!----------------------------------------------------------------
! Given the values of the Zeeman opacities, fill the absorption matrix
! for each depth point
!----------------------------------------------------------------
		subroutine fill_matrix(opacity,ab_matrix)
		real(kind=8) :: opacity(:,:), ab_matrix(:,:,:)

			ab_matrix(1,1,:) = opacity(1,:)   !eta_I
			ab_matrix(2,2,:) = opacity(1,:)   !eta_I
			ab_matrix(3,3,:) = opacity(1,:)   !eta_I
			ab_matrix(4,4,:) = opacity(1,:)   !eta_I
			ab_matrix(1,2,:) = opacity(2,:)   !eta_Q
			ab_matrix(2,1,:) = opacity(2,:)   !eta_Q
			ab_matrix(1,3,:) = opacity(3,:)   !eta_U
			ab_matrix(3,1,:) = opacity(3,:)   !eta_U
			ab_matrix(1,4,:) = opacity(4,:)   !eta_V
			ab_matrix(4,1,:) = opacity(4,:)   !eta_V
			ab_matrix(2,3,:) = opacity(7,:)   !rho_V
			ab_matrix(3,2,:) = -opacity(7,:)  !-rho_V
			ab_matrix(2,4,:) = -opacity(6,:)  !-rho_U
			ab_matrix(4,2,:) = opacity(6,:)   !rho_U
			ab_matrix(3,4,:) = opacity(5,:)   !rho_Q
			ab_matrix(4,3,:) = -opacity(5,:)  !-rho_Q

		end subroutine fill_matrix

! ---------------------------------------------------------
! Performs the formal solution of the RT equation for a plane-parallel magnetized atmosphere
! for a given source function and opacity
! ---------------------------------------------------------
	function formal_sol_polarized(height, opacity, source, mu, boundary, idir)
	real(kind=8) :: height(:), opacity(:,:), source(:), mu, boundary(:)
	real(kind=8) :: formal_sol_polarized(4), Inten(4)
	integer :: k, km, kp, idir, k0, kf
	real(kind=8) :: chim, chi0, chip, dtp, dtm, exu
	real(kind=8) :: psim, psi0, psip, psim_lin, psi0_lin, dm, dp

	integer :: i, j, n
	real(kind=8), allocatable :: ab_matrix(:,:,:), source_vector(:,:)
	real(kind=8) :: sm(4), s0(4), sp(4), mat1(4,4), mat2(4,4), mat3(4,4)

		n = size(height)

		allocate(ab_matrix(4,4,n))
		allocate(source_vector(4,n))

		call fill_matrix(opacity,ab_matrix)

! Transform K into K* and then into K'
		do i = 1, 4
			do j = 1, 4
				ab_matrix(i,j,:) = ab_matrix(i,j,:) / opacity(1,:)
			enddo
			ab_matrix(i,i,:) = ab_matrix(i,i,:) - 1.d0
			source_vector(i,:) = opacity(i,:) * source / opacity(1,:)
		enddo

! Boundary condition
		if (idir == 1) then
			k0 = 2
			kf = n
			Inten = boundary
		endif
		if (idir == -1) then
			k0 = n-1
			kf = 1
			Inten = boundary
		endif

		do k = k0, kf, idir

! Parabolic short-characteristics
			if (k /= kf) then
				km = k - idir
				kp = k + idir
				chim = opacity(1,km)
				chi0 = opacity(1,k)
				chip = opacity(1,kp)
				sm = source_vector(:,km)
				s0 = source_vector(:,k)
				sp = source_vector(:,kp)
				dm = dabs(1.d0/mu*(height(k) - height(km)))
				dp = dabs(1.d0/mu*(height(kp) - height(k)))
			else
! Linear short-characteristics
				km = k - idir
				chim = opacity(1,km)
				chi0 = opacity(1,k)
				chip = 0.d0
				sm = source_vector(:,km)
				s0 = source_vector(:,k)
				sp = 0.d0
				dm = dabs(1.d0/mu*(height(k) - height(km)))
				dp = 0.d0
			endif

			dtm = 0.5d0 * (chim + chi0) * dm
			dtp = 0.5d0 * (chi0 + chip) * dp

			if (dtm >= 1.d-4) then
				exu = dexp(-dtm)
			else
				exu = 1.d0 - dtm + 0.5d0 * dtm**2.d0
			endif

			call lin_sc(dtm,psim_lin,psi0_lin)
			mat1 = exu * identity - psim_lin*ab_matrix(:,:,km)
			mat2 = identity + psi0_lin * ab_matrix(:,:,k)
			call invert(mat2)		

			if (k /= kf) then
				call par_sc(dtm,dtp,psim,psi0,psip)
				Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0 + psip*sp)
			else
				call lin_sc(dtm,psim,psi0)
				Inten = matmul(mat2,matmul(mat1,Inten) + psim*sm + psi0*s0)
			endif

		enddo

		deallocate(ab_matrix)
		deallocate(source_vector)

		formal_sol_polarized = Inten

	end function formal_sol_polarized


end module mathsModule