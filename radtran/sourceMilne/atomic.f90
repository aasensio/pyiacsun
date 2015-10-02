module atomic_functions
use math_functions
use vars
implicit none
contains

!-----------------------------------------------------------------
! Return the profiles weighted by the strength of the components for a given frequency
! It returns zeeman_profile(q,n_depths) with
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------	
	subroutine zeeman_profile(Stokes_Syn,model,linea,zeeman_voigt,zeeman_faraday)
	type(stokes_type) :: Stokes_Syn
	type(modelo_type) :: model
	type(line_type) :: linea
	real(kind=8) :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	integer :: n, nlow, nup, iup, ilow, i_pi, i_blue, i_red, cual
	real(kind=8) :: Mup, Mlow, strength, va, vb, splitting
	
		nup = 2*linea%Jup+1
		nlow = 2*linea%Jlow+1
		
		zeeman_voigt = 0.d0
		zeeman_faraday = 0.d0
		
		i_red = 0
		i_pi = 0
		i_blue = 0
		
		v = (Stokes_Syn%lambda-linea%wave0) / model%doppler
		va = 1.d5 * linea%wave0 * model%vmac / (PC*model%doppler)
		vb = linea%wave0**2 * model%Bfield * 4.6686d-13 / model%doppler
						
		do iup = 1, nup
			Mup = linea%Jup + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= linea%Jlow) then
										
					strength = strength_zeeman(linea%Jup,linea%Jlow,Mup,Mlow)
					splitting = linea%gup*Mup - linea%glow*Mlow
					
					profile = fvoigt_zeeman(model%damping,v-va+vb*splitting)
					zeeman_voigt(ilow,:) = zeeman_voigt(ilow,:) + strength * profile(1,:) / sqrt(PI)
					zeeman_faraday(ilow,:) = zeeman_faraday(ilow,:) + strength * profile(2,:) / sqrt(PI)		
					
				endif
			enddo
		enddo	
				
	end subroutine zeeman_profile
		
!-----------------------------------------------------------------
! Return the seven independent elements of the absorption matrix
! Remember that, zeeman_voigt(q,:) and zeeman_faraday(q,:) have
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------	
	subroutine zeeman_opacity(model,zeeman_voigt,zeeman_faraday,ki,kq,ku,kv,fq,fu,fv)
	type(modelo_type) :: model
	integer :: line
	real(kind=8) :: zeeman_voigt(:,:), zeeman_faraday(:,:)
	real(kind=8) :: ki(:), kq(:), ku(:), kv(:), fq(:), fu(:), fv(:)
	real(kind=8) :: sin_theta, cos_theta, sin_2chi, cos_2chi

		sin_theta = sin(model%theta * PI / 180.d0)
		cos_theta = cos(model%theta * PI / 180.d0)
		
		sin_2chi = sin(2.d0 * model%chi * PI / 180.d0)
		cos_2chi = cos(2.d0 * model%chi * PI / 180.d0)		

! Classical absorption coefficients		
		ki = ki + model%kl * 0.5d0 * (zeeman_voigt(2,:)*sin_theta**2 + 0.5d0*(zeeman_voigt(1,:)+zeeman_voigt(3,:))*(1.d0+cos_theta**2))  ! eta_I
		kq = kq + model%kl * 0.5d0 * (zeeman_voigt(2,:) - 0.5d0*(zeeman_voigt(1,:)+zeeman_voigt(3,:))) * sin_theta**2*cos_2chi  ! eta_Q
		ku = ku + model%kl * 0.5d0 * (zeeman_voigt(2,:) - 0.5d0*(zeeman_voigt(1,:)+zeeman_voigt(3,:))) * sin_theta**2*sin_2chi  ! eta_U
		kv = kv + model%kl * 0.5d0 * (zeeman_voigt(3,:)-zeeman_voigt(1,:)) * cos_theta  ! eta_V

! Magneto-optical coefficients		
		fq = fq + model%kl * 0.5d0 * (zeeman_faraday(2,:) - 0.5d0*(zeeman_faraday(1,:)+zeeman_faraday(3,:))) * sin_theta**2*cos_2chi  ! rho_Q
		fu = fu + model%kl * 0.5d0 * (zeeman_faraday(2,:) - 0.5d0*(zeeman_faraday(1,:)+zeeman_faraday(3,:))) * sin_theta**2*sin_2chi  ! rho_U
		fv = fv + model%kl * 0.5d0 * (zeeman_faraday(3,:)-zeeman_faraday(1,:)) * cos_theta  ! rho_V
				
	end subroutine zeeman_opacity	
	
!-----------------------------------------------------------------
! Add the atomic opacity to the opacity including the effect of a magnetic field
!-----------------------------------------------------------------	
	subroutine synthesize(model,linea,Stokes_Syn)
	type(modelo_type) :: model
	type(stokes_type) :: Stokes_Syn
	type(line_type) :: linea(:)
	integer :: i, j, k, n, n_lineas	
	real(kind=8) :: factor1 ,factor2, lmax, lmin, lstep
	character(len=1) :: str					
					
		ki = 0.d0
		kq = 0.d0
		ku = 0.d0
		kv = 0.d0
		fq = 0.d0
		fu = 0.d0
		fv = 0.d0
			
		do i = 1, nLines
			call zeeman_profile(Stokes_Syn,model,linea(i),zeeman_voigt,zeeman_faraday)
							
			call zeeman_opacity(model,zeeman_voigt,zeeman_faraday,ki,kq,ku,kv,fq,fu,fv)
		enddo
																				
		delta = (1.d0+ki)**4 + (1.d0+ki)**2 * (fq**2+fu**2+fv**2-kq**2-ku**2-kv**2) - &
			(kq*fq+ku*fu+kv*fv)**2
		
		stokes(1,:) = model%B0 * (1.d0+model%B1/model%B0*model%mu*(1.d0+ki) / delta * &
			((1.d0+ki)**2 + fq**2 + fu**2 + fv**2))
		stokes(2,:) = -model%B1*model%mu / delta * ((1.d0+ki)**2*kq - (1.d0+ki)*(ku*fv-kv*fu) + &
			fq*(kq*fq+ku*fu+kv*fv))
		stokes(3,:) = -model%B1*model%mu / delta * ((1.d0+ki)**2*ku - (1.d0+ki)*(kv*fq-kq*fv) + &
			fu*(kq*fq+ku*fu+kv*fv))
		stokes(4,:) = -model%B1*model%mu / delta * ((1.d0+ki)**2*kv - (1.d0+ki)*(kq*fu-ku*fq) + &
			fv*(kq*fq+ku*fu+kv*fv))
			
		Stokes_Syn%stokes = stokes
			
	end subroutine synthesize
		
end module atomic_functions
