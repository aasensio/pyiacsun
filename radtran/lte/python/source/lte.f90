module lteMod
use globalModule, only : atmosphere, lineList, PH, PC, PK, PHK, UMA, EV_ERG, SQRTPI
use synthModule, only : synthLines
use mathsModule, only: qsortd
use iso_c_binding, only: c_int, c_double
implicit none

contains

	subroutine c_initatmosphere(nDepths) bind(c)
	integer(c_int), intent(in) :: nDepths
   
		atmosphere%nDepths = nDepths
		
		print *, 'N. depth points = ', atmosphere%nDepths
		
		if (associated(atmosphere%lTau500)) deallocate(atmosphere%lTau500)
		allocate(atmosphere%lTau500(atmosphere%nDepths))
		
		if (associated(atmosphere%lTau500Ordered)) deallocate(atmosphere%lTau500Ordered)
		allocate(atmosphere%lTau500Ordered(atmosphere%nDepths))
		
		if (associated(atmosphere%height)) deallocate(atmosphere%height)
		allocate(atmosphere%height(atmosphere%nDepths))
		
		if (associated(atmosphere%T)) deallocate(atmosphere%T)
		allocate(atmosphere%T(atmosphere%nDepths))
		
		if (associated(atmosphere%TOrdered)) deallocate(atmosphere%TOrdered)
		allocate(atmosphere%TOrdered(atmosphere%nDepths))
		
		if (associated(atmosphere%microturbulence)) deallocate(atmosphere%microturbulence)
		allocate(atmosphere%microturbulence(atmosphere%nDepths))
		
		if (associated(atmosphere%vmac)) deallocate(atmosphere%vmac)
		allocate(atmosphere%vmac(atmosphere%nDepths))
		
		if (associated(atmosphere%niovern)) deallocate(atmosphere%niovern)
		allocate(atmosphere%niovern(atmosphere%nDepths))
		
		if (associated(atmosphere%ui)) deallocate(atmosphere%ui)
		allocate(atmosphere%ui(atmosphere%nDepths))
		
		if (associated(atmosphere%u1)) deallocate(atmosphere%u1)
		allocate(atmosphere%u1(atmosphere%nDepths))
		
		if (associated(atmosphere%u2)) deallocate(atmosphere%u2)
		allocate(atmosphere%u2(atmosphere%nDepths))
		
		if (associated(atmosphere%u3)) deallocate(atmosphere%u3)
		allocate(atmosphere%u3(atmosphere%nDepths))
		
		if (associated(atmosphere%n1overn0)) deallocate(atmosphere%n1overn0)
		allocate(atmosphere%n1overn0(atmosphere%nDepths))
		
		if (associated(atmosphere%n2overn1)) deallocate(atmosphere%n2overn1)
		allocate(atmosphere%n2overn1(atmosphere%nDepths))
		
		if (associated(atmosphere%Pe)) deallocate(atmosphere%Pe)
		allocate(atmosphere%Pe(atmosphere%nDepths))
		
		if (associated(atmosphere%Pg)) deallocate(atmosphere%Pg)
		allocate(atmosphere%Pg(atmosphere%nDepths))
		
		if (associated(atmosphere%PH)) deallocate(atmosphere%PH)
		allocate(atmosphere%PH(atmosphere%nDepths))
		
		if (associated(atmosphere%PHminus)) deallocate(atmosphere%PHminus)
		allocate(atmosphere%PHminus(atmosphere%nDepths))
		
		if (associated(atmosphere%PHplus)) deallocate(atmosphere%PHplus)
		allocate(atmosphere%PHplus(atmosphere%nDepths))
		
		if (associated(atmosphere%PH2)) deallocate(atmosphere%PH2)
		allocate(atmosphere%PH2(atmosphere%nDepths))
		
		if (associated(atmosphere%PH2plus)) deallocate(atmosphere%PH2plus)
		allocate(atmosphere%PH2plus(atmosphere%nDepths))
		
		if (associated(atmosphere%PTotal)) deallocate(atmosphere%PTotal)
		allocate(atmosphere%PTotal(atmosphere%nDepths))
		
		if (associated(atmosphere%nHtot)) deallocate(atmosphere%nHtot)
		allocate(atmosphere%nHtot(atmosphere%nDepths))
		
		if (associated(atmosphere%order)) deallocate(atmosphere%order)
		allocate(atmosphere%order(atmosphere%nDepths))		
		
		if (associated(atmosphere%opacity500)) deallocate(atmosphere%opacity500)
		allocate(atmosphere%opacity500(atmosphere%nDepths))
				
	end subroutine c_initatmosphere
   
	subroutine c_initlines(nLines, lineListIn, nLambda, lambdaAxis) bind(c)
	integer(c_int), intent(in) :: nLines, nLambda
	real(c_double), intent(in) :: lineListIn(11,nLines)
	real(c_double), intent(in) :: lambdaAxis(nLambda)
	integer :: activeLines, i
   
		lineList%nLines = nLines
		
		if (associated(linelist%transition)) deallocate(lineList%transition)
		allocate(lineList%transition(lineList%nLines))
		
		activeLines = 0
		
! 		lineList%transition(i)%lambdaLeft, lineList%transition(i)%lambdaRight
		
		do i = 1, lineList%nLines
			lineList%transition(i)%nLambda = nLambda			
			lineList%transition(i)%lambda0 = lineListIn(1,i)
			lineList%transition(i)%element = lineListIn(2,i)
			lineList%transition(i)%ionization = lineListIn(3,i)
			lineList%transition(i)%gf = lineListIn(4,i)
			lineList%transition(i)%Elow = lineListIn(5,i)
			lineList%transition(i)%gu = lineListIn(6,i)
			lineList%transition(i)%gl = lineListIn(7,i)
			lineList%transition(i)%Ju = lineListIn(8,i)
			lineList%transition(i)%Jl = lineListIn(9,i)
			lineList%transition(i)%sigmaABO = lineListIn(10,i)
			lineList%transition(i)%alphaABO = lineListIn(11,i)
						
			activeLines = activeLines + 1
			
			lineList%transition(i)%gf = 10.d0**lineList%transition(i)%gf
			lineList%transition(i)%Elow = lineList%transition(i)%Elow * PH * PC         ! Transform to erg
			lineList%transition(i)%frequency0 = PC / (lineList%transition(i)%lambda0 * 1.d-8)
						
			allocate(lineList%transition(i)%lineOpacity(atmosphere%nDepths))
			allocate(lineList%transition(i)%backOpacity(atmosphere%nDepths))
			allocate(lineList%transition(i)%dopplerWidth(atmosphere%nDepths))
			allocate(lineList%transition(i)%deltaNu(atmosphere%nDepths))
			allocate(lineList%transition(i)%damping(atmosphere%nDepths))						
			
			allocate(lineList%transition(i)%lambda(lineList%transition(i)%nLambda))
			allocate(lineList%transition(i)%frequency(lineList%transition(i)%nLambda))
			allocate(lineList%transition(i)%intensity(lineList%transition(i)%nLambda))
			allocate(lineList%transition(i)%continuumIntensity(lineList%transition(i)%nLambda))
			allocate(lineList%transition(i)%observed(lineList%transition(i)%nLambda))
				
			allocate(lineList%transition(i)%opacity(atmosphere%nDepths,lineList%transition(i)%nLambda))
			allocate(lineList%transition(i)%opacityContinuum(atmosphere%nDepths,lineList%transition(i)%nLambda))
			allocate(lineList%transition(i)%source(atmosphere%nDepths,lineList%transition(i)%nLambda))
			
			lineList%transition(i)%lambda = lambdaAxis
			lineList%transition(i)%frequency = PC / (lineList%transition(i)%lambda * 1.d-8)
							
		enddo
		
		allocate(lineList%opacity(atmosphere%nDepths,nLambda))
		allocate(lineList%source(atmosphere%nDepths,nLambda))
		allocate(lineList%boundary(nLambda))
		allocate(lineList%opacityContinuum(atmosphere%nDepths,nLambda))
		allocate(lineList%intensity(nLambda))
		allocate(lineList%intensityContinuum(nLambda))
				
		write(*,*) 'Number of active lines : ', activeLines
	end subroutine c_initlines
   
	subroutine c_synthlines(nDepths, atmosphereIn, nLambda, stokesOut) bind(c)
	integer(c_int), intent(in) :: nDepths, nLambda
	real(c_double), intent(in) :: atmosphereIn(4,nDepths)
	real(c_double), intent(inout) :: stokesOut(5,nLambda)
	integer :: i
   		
		atmosphere%lTau500 = atmosphereIn(1,:)
		atmosphere%T = atmosphereIn(2,:)
		atmosphere%microturbulence= atmosphereIn(3,:)
		atmosphere%vmac = atmosphereIn(4,:)
						
! Order the optical depth
		! call qsortd(atmosphere%lTau500,atmosphere%order,atmosphere%nDepths)
		
! Synthesize lines				
		call synthLines(atmosphere, stokesOut)
		
	end subroutine c_synthlines
   
end module lteMod