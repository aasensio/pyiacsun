module ioModule
use globalModule, only : atmosphereType, atlasType, lineListType, PH, PC
use mathsModule, only: extremes, qsortd
use hydrostaticModule, only : hydrostaticEquilibrium
use atomicPartitionModule, only : partitionAtomic
implicit none
contains

!------------------------------------------------
! Read input data
!------------------------------------------------
	subroutine readInput(atmosphereIn, lineListIn)
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList
	type(atlasType) :: atlas
	integer :: j
	
		open(unit=12,file='turb.conf',action='read',status='old')
		
		read(12,*)
		read(12,*) atmosphere%file

! Read the atmosphere model
		call readAtmosphere(atmosphere)
		
		read(12,*)
		read(12,*)
		read(12,*) lineList%file
		
! Read the line list
		call readLineList(lineList)
		
		read(12,*)
		read(12,*)
		read(12,*) atlas%file
		
! Read the atlas
		call readAtlas(atlas)
		
! Read the nodes for the temperature
		read(12,*)
		read(12,*)
		read(12,*) atmosphere%nNodes		
		allocate(atmosphere%nodes(atmosphere%nNodes))		
		read(12,*) (atmosphere%nodes(j),j=1,atmosphere%nNodes)
		
		close(12)
		
! Set-up line information
 		call setupLines(atmosphere,lineList)
		
	end subroutine readInput
	
!------------------------------------------------
! Read the model atmosphere
!------------------------------------------------
	subroutine readAtmosphere(nDepths, atmosphereIn, atmosphere)
	real(c_double), intent(in) :: atmosphereIn(3,nDepths)
	type(atmosphereType) :: atmosphere
	integer :: i
	real (kind=8) :: u1, u2, u3, ei1, ei2, weight	
	
		open(unit=13,file=trim(adjustl(atmosphere%file)),action='read',status='old')
		read(13,*) atmosphere%nDepths
		
		allocate(atmosphere%lTau500(atmosphere%nDepths))
		allocate(atmosphere%lTau500Ordered(atmosphere%nDepths))
		allocate(atmosphere%height(atmosphere%nDepths))
		allocate(atmosphere%T(atmosphere%nDepths))
		allocate(atmosphere%TOrdered(atmosphere%nDepths))
		allocate(atmosphere%microturbulence(atmosphere%nDepths))
		allocate(atmosphere%niovern(atmosphere%nDepths))
		allocate(atmosphere%ui(atmosphere%nDepths))
		allocate(atmosphere%u1(atmosphere%nDepths))
		allocate(atmosphere%u2(atmosphere%nDepths))
		allocate(atmosphere%u3(atmosphere%nDepths))
		allocate(atmosphere%n1overn0(atmosphere%nDepths))
		allocate(atmosphere%n2overn1(atmosphere%nDepths))
		
		
		allocate(atmosphere%Pe(atmosphere%nDepths))
		allocate(atmosphere%Pg(atmosphere%nDepths))
		allocate(atmosphere%PH(atmosphere%nDepths))
		allocate(atmosphere%PHminus(atmosphere%nDepths))
		allocate(atmosphere%PHplus(atmosphere%nDepths))
		allocate(atmosphere%PH2(atmosphere%nDepths))
		allocate(atmosphere%PH2plus(atmosphere%nDepths))
		allocate(atmosphere%PTotal(atmosphere%nDepths))
		allocate(atmosphere%nHtot(atmosphere%nDepths))
		allocate(atmosphere%order(atmosphere%nDepths))		
		allocate(atmosphere%opacity500(atmosphere%nDepths))
		
		do i = 1, atmosphere%nDepths
			read(13,*) atmosphere%lTau500(i), atmosphere%T(i), atmosphere%microturbulence(i)
		enddo
		
		close(13)
		
! Iron abundance
		call partitionAtomic(atmosphere%T(1), 26, u1, u2, u3, ei1, ei2, weight, atmosphere%abundance)
				
! Order the optical depth
		call qsortd(atmosphere%lTau500,atmosphere%order,atmosphere%nDepths)
		
! Put the atmosphere in hydrostatic equilibrium		
!  		call hydrostaticEquilibrium(atmosphere)
		
	end subroutine readAtmosphere
	
!------------------------------------------------
! Read the model atmosphere
!------------------------------------------------
	subroutine readLineList(lineList)
	integer :: i, active, activeLines
	type(lineListType) :: lineList
	
		open(unit=13,file=trim(adjustl(lineList%file)),action='read',status='old')
		read(13,*) lineList%nLines
		
		allocate(lineList%transition(lineList%nLines))
		
		activeLines = 0
				
		do i = 1, lineList%nLines
			read(13,*) active, lineList%transition(i)%lambda0, lineList%transition(i)%lambdaNew, lineList%transition(i)%lambdaKurucz, lineList%transition(i)%gf, lineList%transition(i)%Elow, &
				lineList%transition(i)%geff, lineList%transition(i)%Gt, lineList%transition(i)%sigmaABO,&
				lineList%transition(i)%alphaABO, lineList%transition(i)%lambdaLeft, lineList%transition(i)%lambdaRight, lineList%transition(i)%lambdaContLeft, lineList%transition(i)%lambdaContRight
			
			if (active == 1 .and. lineList%transition(i)%lambda0 /= 0.d0) then
				activeLines = activeLines + 1
			else						
				lineList%transition(i)%lambda0 = 0.d0
			endif
			
			lineList%transition(i)%gf = 10.d0**lineList%transition(i)%gf
			lineList%transition(i)%Elow = lineList%transition(i)%Elow * PH * PC         ! Transform to erg
			lineList%transition(i)%frequency0 = PC / (lineList%transition(i)%lambdaNew * 1.d-8)
		enddo
		
		close(13)
		
		write(*,*) 'Number of active lines : ', activeLines
		
		
	end subroutine readLineList
		
!------------------------------------------------
! Find indices
!------------------------------------------------
	subroutine setupLines(atmosphere,lineList,atlas)
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList
	type(atlasType) :: atlas
	integer :: i, limits(2), nTotal
	real(kind=8) :: a, b, x1, y1, x2, y2
		
		lineList%nLambdaTotal = 0
		nTotal = 0
		
		do i = 1, lineList%nLines
			
			lineList%transition(i)%active = .FALSE.
			
			if (lineList%transition(i)%lambda0 /= 0.d0) then
			
				nTotal = nTotal + 1
! Get ranges of the line
				lineList%transition(i)%active = .TRUE.
				limits = extremes(atlas%lambda, lineList%transition(i)%lambdaLeft, lineList%transition(i)%lambdaRight)
				write(*,FMT='(A,I3,A,I3,A,I7,A,I7,A,F8.3,A,F8.3)') 'Line ', i, '(', nTotal, ') - Range : ', limits(1), ' -> ', limits(2), ' - lambda: ', atlas%lambda(limits(1)), ' -> ', atlas%lambda(limits(2))
				
				lineList%transition(i)%nLambda = limits(2) - limits(1) + 1
				allocate(lineList%transition(i)%lambda(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%frequency(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%intensity(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%continuumIntensity(lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%observed(lineList%transition(i)%nLambda))
				
				allocate(lineList%transition(i)%opacity(atmosphere%nDepths,lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%opacityContinuum(atmosphere%nDepths,lineList%transition(i)%nLambda))
				allocate(lineList%transition(i)%source(atmosphere%nDepths,lineList%transition(i)%nLambda))
				
				allocate(lineList%transition(i)%lineOpacity(atmosphere%nDepths))
				allocate(lineList%transition(i)%backOpacity(atmosphere%nDepths))
				allocate(lineList%transition(i)%dopplerWidth(atmosphere%nDepths))
				allocate(lineList%transition(i)%deltaNu(atmosphere%nDepths))
				allocate(lineList%transition(i)%damping(atmosphere%nDepths))
				
				lineList%transition(i)%lambda = atlas%lambda(limits(1):limits(2))
				lineList%transition(i)%observed = atlas%intensity(limits(1):limits(2))
				lineList%transition(i)%frequency = PC / (lineList%transition(i)%lambda * 1.d-8)
				
				lineList%nLambdaTotal = lineList%nLambdaTotal + lineList%transition(i)%nLambda				
				
			endif
		
		enddo
		
		print *, 'Total number of wavelengths : ', lineList%nLambdaTotal
		
		allocate(lineList%opacity(atmosphere%nDepths,lineList%nLambdaTotal))
		allocate(lineList%source(atmosphere%nDepths,lineList%nLambdaTotal))
		allocate(lineList%boundary(lineList%nLambdaTotal))
		allocate(lineList%opacityContinuum(atmosphere%nDepths,lineList%nLambdaTotal))
		allocate(lineList%intensity(lineList%nLambdaTotal))
		allocate(lineList%intensityContinuum(lineList%nLambdaTotal))
				
	end subroutine setupLines

	
end module ioModule