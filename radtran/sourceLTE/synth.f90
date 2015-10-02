module synthModule
use globalModule, only : OPA, PH, PC, PK, PHK, UMA, EV_ERG, SQRTPI, LARMOR, atmosphereType, lineListType, lineList, transitionType
use atomicPartitionModule, only : partitionAtomic
use mathsModule, only : saha, vecfvoigt, calculateDamping, planckFrequency, shortCharacteristics, computeHeight, gasPressure, vecfvoigt_zeeman, strength_zeeman, formal_sol_polarized
use backgroundOpacityModule, only : backgroundOpacity
use hydrostaticModule, only : hydrostaticEquilibrium
implicit none
contains

!-----------------------------------------------------------------
! Fill the strength and splitting variables of the atomic_transition structure
!-----------------------------------------------------------------
	subroutine generate_atomic_zeeman_components(transition)
	type(transitionType) :: transition
	integer :: i, n, nlow, nup, iup, ilow, i_pi, i_red, i_blue, cual
	real(kind=8) :: Mup, Mlow, strength

		nup = 2*transition%Ju+1
		nlow = 2*transition%Jl+1

		i_pi = 0
		i_blue = 0
		i_red = 0

! First count the number of components of each type (
		transition%nComponents = 0
		do iup = 1, nup
			Mup = transition%Ju + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= transition%Jl) then
					transition%nComponents(ilow) = transition%nComponents(ilow)+1
				endif
			enddo
		enddo

		do i = 1, 3
			allocate(transition%ZeemanComponent(i)%splitting(transition%nComponents(i)))
			allocate(transition%ZeemanComponent(i)%strength(transition%nComponents(i)))
		enddo

! Now generate all the data
		do iup = 1, nup
			Mup = transition%Ju + 1 - iup  ! Mup=J...-J
			do ilow = 1, 3
				Mlow = Mup-2+ilow			! Mlow=Mup-1,Mup,Mup+1 (the allowed transitions)
				if (abs(Mlow) <= transition%Jl) then
					if (ilow == 1) then
						i_blue = i_blue + 1
						cual = i_blue
					endif
					if (ilow == 2) then
						i_pi = i_pi + 1
						cual = i_pi
					endif
					if (ilow == 3) then
						i_red = i_red + 1
						cual = i_red
					endif

					transition%ZeemanComponent(ilow)%strength(cual) = strength_zeeman(transition%Ju,transition%Jl,Mup,Mlow)
					transition%ZeemanComponent(ilow)%splitting(cual) = (transition%gu*Mup - transition%gl*Mlow)

				endif
			enddo
		enddo

	end subroutine generate_atomic_zeeman_components

!-----------------------------------------------------------------
! Return the seven independent elements of the absorption matrix
! Remember that, zeeman_voigt(q,:) and zeeman_faraday(q,:) have
!  q=1  Mlow=Mup-1  (sigma blue)
!  q=2  Mlow=Mup    (sigma pi)
!  q=3  Mlow=Mup+1  (sigma red)
!-----------------------------------------------------------------
	subroutine zeeman_opacity(transition, atmos, frequency)	
	real(kind=8) :: frequency
	type(atmosphereType) :: atmos
	type(transitionType) :: transition	
	integer :: n, nlow, nup, iup, ilow, i_pi, i_blue, i_red, cual
	real(kind=8) :: Mup, Mlow, strength
	
		nup = 2*transition%Ju+1
		nlow = 2*transition%Jl+1

		atmos%zeeman_voigt = 0.d0
		atmos%zeeman_faraday = 0.d0
	

		do ilow = 1, 3
			do cual = 1, transition%nComponents(ilow)
				atmos%splitting = LARMOR * atmos%B * transition%ZeemanComponent(ilow)%splitting(cual)
				atmos%profile = vecfvoigt_zeeman(transition%damping, (transition%frequency0 - frequency - transition%frequency0 * 1e5 * atmos%vmac / PC + atmos%splitting) / transition%deltaNu)

				atmos%zeeman_voigt(:,ilow) = atmos%zeeman_voigt(:,ilow) + transition%ZeemanComponent(ilow)%strength(cual) * atmos%profile(:,1)
				atmos%zeeman_faraday(:,ilow) = atmos%zeeman_faraday(:,ilow) + transition%ZeemanComponent(ilow)%strength(cual) * atmos%profile(:,2)

			enddo
		enddo	

! Classical absorption coefficients
		atmos%coefficients(:,1) = 0.5d0 * (atmos%zeeman_voigt(:,2)*atmos%sinthetaB2 + &
			0.5d0*(atmos%zeeman_voigt(:,1)+atmos%zeeman_voigt(:,3))*(1.d0+atmos%costhetaB2))  ! eta_I
 		atmos%coefficients(:,2) = 0.5d0 * (atmos%zeeman_voigt(:,2) - &
 			0.5d0*(atmos%zeeman_voigt(:,1)+atmos%zeeman_voigt(:,3))) * atmos%sinthetaB2 * atmos%cos2chiB  ! eta_Q
 		atmos%coefficients(:,3) = 0.5d0 * (atmos%zeeman_voigt(:,2) - &
 			0.5d0*(atmos%zeeman_voigt(:,1)+atmos%zeeman_voigt(:,3))) * atmos%sinthetaB2 * atmos%sin2chiB  ! eta_U
 		atmos%coefficients(:,4) = 0.5d0 * (atmos%zeeman_voigt(:,3)-atmos%zeeman_voigt(:,1)) * atmos%costhetaB         ! eta_V
! 
! ! Magneto-optical coefficients
 		atmos%coefficients(:,5) = 0.5d0 * (atmos%zeeman_faraday(:,2) - &
 			0.5d0*(atmos%zeeman_faraday(:,1)+atmos%zeeman_faraday(:,3))) * atmos%sinthetaB2 * atmos%cos2chiB  ! rho_Q
 		atmos%coefficients(:,6) = 0.5d0 * (atmos%zeeman_faraday(:,2) - &
 			0.5d0*(atmos%zeeman_faraday(:,1)+atmos%zeeman_faraday(:,3))) * atmos%sinthetaB2 * atmos%sin2chiB  ! rho_U
 		atmos%coefficients(:,7) = 0.5d0 * (atmos%zeeman_faraday(:,3)-atmos%zeeman_faraday(:,1)) * atmos%costhetaB  ! rho_V

	end subroutine zeeman_opacity
	
!------------------------------------------------
! Synthesize all lines
!------------------------------------------------
	subroutine synthLines(atmosphere, stokesOut)
	type(atmosphereType) :: atmosphere	
	real(kind=8) :: n1overn0, n2overn1, ei1, ei2, weight, abundance, stokesOut(:,:)
	real(kind=8), allocatable :: voigtProfile(:)
	integer :: i, j, loop, fromLambda, toLambda, k
				
		if (.not.allocated(voigtProfile)) allocate(voigtProfile(atmosphere%nDepths))
	
! Put the atmosphere in hydrostatic equilibrium to compute the electron pressure		
 		call hydrostaticEquilibrium(atmosphere)
 		 		 		
! Compute the chemical equilibrium
		call gasPressure(atmosphere%Pe, atmosphere%T, atmosphere%PH, atmosphere%PHminus, atmosphere%PHplus, atmosphere%PH2, atmosphere%PH2plus, atmosphere%PTotal)							

! Total hydrogen density
		atmosphere%nhtot = (atmosphere%PH + atmosphere%PHminus + atmosphere%PHplus + atmosphere%PH2 + atmosphere%PH2plus) / (PK * atmosphere%T)
		
! Generate height axis
		do j = 1, atmosphere%nDepths
			atmosphere%opacity500(j) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
				atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), 5000.d0)		
		enddo
		
! Compute the height depth scale
		call computeHeight(atmosphere%lTau500, atmosphere%opacity500, atmosphere%height)
		
		fromLambda = 1
		loop = 1
		
		lineList%opacity = 0.d0
		lineList%boundary = 0.d0
		lineList%source = 0.d0
		
		do i = 1, lineList%nLines
			do j = 1, atmosphere%nDepths
				call partitionAtomic(atmosphere%T(j), lineList%transition(i)%element, atmosphere%u1(j), atmosphere%u2(j), &
					atmosphere%u3(j), ei1, ei2, weight, atmosphere%abundance)
				lineList%transition(i)%backOpacity(j) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
					atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), lineList%transition(i)%lambda0)
			enddo
			
			atmosphere%n1overn0 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u1, atmosphere%u2, ei1)
			atmosphere%n2overn1 = saha(atmosphere%T, atmosphere%Pe, atmosphere%u2, atmosphere%u3, ei2)
			atmosphere%niovern = 1.d0 / (1.d0 + atmosphere%n1overn0 + atmosphere%n2overn1 * atmosphere%n1overn0)						
						
			select case (lineList%transition(i)%ionization)
				case(1)
					atmosphere%ui = atmosphere%u1
				case(2)
					atmosphere%ui = atmosphere%u2
				case(3)
					atmosphere%ui = atmosphere%u3
			end select
		
! Compute line opacity
			lineList%transition(i)%lineOpacity = OPA * lineList%transition(i)%gf / atmosphere%ui * dexp(-lineList%transition(i)%Elow / (PK * atmosphere%T)) *&
				(1.d0 - dexp(-PHK * lineList%transition(i)%frequency0 / atmosphere%T)) * atmosphere%niovern * &
				(atmosphere%nhtot * atmosphere%abundance)
												
! Doppler width
			lineList%transition(i)%dopplerWidth = dsqrt(atmosphere%microturbulence**2.d0 + 2.d0 * PK * atmosphere%T / (weight * UMA))
			lineList%transition(i)%deltaNu = lineList%transition(i)%dopplerWidth * lineList%transition(i)%frequency0 / PC
			
! Background opacity
			do j = 1, atmosphere%nDepths
				lineList%transition(i)%backOpacity(j) = backgroundOpacity(atmosphere%T(j), atmosphere%Pe(j), atmosphere%PH(j), atmosphere%PHminus(j), &
					atmosphere%PHplus(j), atmosphere%PH2(j), atmosphere%PH2plus(j), lineList%transition(i)%lambda0)
			enddo
			
! Compute the damping
			lineList%transition(i)%damping = calculateDamping(atmosphere%T, atmosphere%nHtot, atmosphere%Pe, atmosphere%Pg, &
				lineList%transition(i)%dopplerWidth, ei1, lineList%transition(i)%Elow / EV_ERG, 0.d0, 0.d0, 0.d0, lineList%transition(i)%lambda0, &
				lineList%transition(i)%alphaABO, lineList%transition(i)%sigmaABO, weight)
								
! Compute frequency dependency of line opacity, source function and boundary condition			
			lineList%transition(i)%opacity = 0.d0
			do j = 1, lineList%transition(i)%nLambda
				call zeeman_opacity(lineList%transition(i), atmosphere, lineList%transition(i)%frequency(j))
				do k = 1, 7
					lineList%transition(i)%opacity(k,:,j) = lineList%transition(i)%opacity(k,:,j) + lineList%transition(i)%lineOpacity * atmosphere%coefficients(:,k) / (lineList%transition(i)%deltaNu * SQRTPI)
				enddo
				lineList%transition(i)%opacityContinuum(:,j) = lineList%transition(i)%backOpacity								
				lineList%transition(i)%source(:,j) = planckFrequency(lineList%transition(i)%frequency(j), atmosphere%T)				
				lineList%transition(i)%boundary = lineList%transition(i)%source(atmosphere%nDepths,j)
			enddo

! Add to total opacity
			lineList%opacity = lineList%opacity + lineList%transition(i)%opacity
			lineList%opacityContinuum = lineList%transition(i)%opacityContinuum
			lineList%source = lineList%transition(i)%source
			lineList%boundary(1,:) = lineList%transition(i)%boundary

		enddo
		
! Add continuum opacity
		lineList%opacity(1,:,:) = lineList%opacity(1,:,:) + lineList%opacityContinuum
						
! Solve RT equation
		do i = 1, lineList%nLambdaTotal
			stokesOut(1:4,i) = formal_sol_polarized(atmosphere%height, lineList%opacity(:,:,i), lineList%source(:,i), 1.d0, lineList%boundary(:,i), -1)			
		enddo
		
  		stokesOut(5,:) = shortCharacteristics(atmosphere%height, lineList%opacityContinuum, lineList%source, 1.d0, lineList%boundary(1,:), -1)
				
	end subroutine synthLines
end module synthModule