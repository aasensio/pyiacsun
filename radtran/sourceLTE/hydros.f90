module hydrostaticModule
use globalModule, only : atmosphereType, GRAV
use atomicPartitionModule, only : partitionAtomic
use mathsModule, only : saha, deriv, gasPressure
use backgroundOpacityModule, only : backgroundOpacityPerHParticle
implicit none
contains

!------------------------------------------------
! Hydrostatic equilibrium
!------------------------------------------------
	subroutine hydrostaticEquilibrium(atmosphere)
	type(atmosphereType) :: atmosphere
	real(kind=8) :: step
	real(kind=8), allocatable :: phi(:,:), tdphi(:,:), u1(:), u2(:), u3(:), one(:), dtdx(:), y(:), dy(:)
	integer :: electronRelease(10) = (/1,2,6,8,11,12,13,14,16,26/)
	real(kind=8) :: abundances(10)
	real(kind=8) :: abundancesRatio(10), ei1, ei2, weight
	real(kind=8) :: PH, PHminus, PHplus, PH2, PH2plus, PTotal, pe0, k_c, delta, pe2, dpdgt, der, mu, Pg
	integer :: i, loop
	
		atmosphere%TOrdered = atmosphere%T(atmosphere%order)
		atmosphere%lTau500Ordered = atmosphere%lTau500(atmosphere%order)
		step = (atmosphere%lTau500Ordered(2) - atmosphere%lTau500Ordered(1))
		
		if (allocated(phi)) deallocate(phi)
		if (allocated(tdphi)) deallocate(tdphi)
		if (allocated(u1)) deallocate(u1)
		if (allocated(u2)) deallocate(u2)
		if (allocated(u3)) deallocate(u3)
		if (allocated(one)) deallocate(one)
		if (allocated(dtdx)) deallocate(dtdx)
		if (allocated(dy)) deallocate(dy)
		if (allocated(y)) deallocate(y)
		
		allocate(phi(10,atmosphere%nDepths))
		allocate(tdphi(10,atmosphere%nDepths))
		allocate(u1(atmosphere%nDepths))
		allocate(u2(atmosphere%nDepths))
		allocate(u3(atmosphere%nDepths))
		allocate(one(atmosphere%nDepths))
		allocate(dtdx(atmosphere%nDepths))
		allocate(dy(atmosphere%nDepths))
		allocate(y(atmosphere%nDepths))
		
		one = 1.d0
		
		
		do i = 1, 10
			call partitionAtomic(atmosphere%TOrdered, electronRelease(i), u1, u2, u3, ei1, ei2, weight, abundances(i))			
			phi(i,:) = saha(atmosphere%TOrdered, one, u1, u2, ei1)
			tdphi(i,:) = saha(1.001d0*atmosphere%TOrdered, one, u1, u2, ei1)
		enddo		
		
		tdphi = (tdphi - phi) / 0.001d0
		abundancesRatio = abundances / 1.10165d0
						
		call deriv(atmosphere%lTau500Ordered,log10(atmosphere%TOrdered),dtdx)			
				
! Integrate the hydrostatic equilibrium equation
		pe0 = 0.001d0
		mu = 2.3848d-24  ! Mass per H particle
		dy = 0.d0
		
		call gasPressure(pe0, atmosphere%T(1), PH, PHminus, PHplus, PH2, PH2plus, PTotal)
		k_c = backgroundOpacityPerHParticle(atmosphere%T(1), pe0, PH, PHminus, PHplus, PH2, PH2plus, 5000.d0) / mu		
		atmosphere%Pe(1) = sqrt(GRAV * 10.d0**atmosphere%lTau500Ordered(1) / k_c * pe0)		
		delta = abs(atmosphere%Pe(1) / pe0 - 1.d0)
		pe0 = atmosphere%Pe(1)
						
		loop = 0
		do while (delta > 0.1d0)
			pe2 = atmosphere%Pe(1) / 2.d0
			atmosphere%Pg(1) = sum(abundancesRatio * phi(:,1) / (phi(:,1) + pe2))
			atmosphere%Pg(1) = pe2 * (1.d0 + 1.d0 / atmosphere%Pg(1))
			dpdgt = sum(abundancesRatio * tdphi(:,1) / (phi(:,1) + pe2)**2)
			dpdgt = (atmosphere%Pg(1) - pe2)**2 * dpdgt * dtdx(1)
			der = sum(abundancesRatio * phi(:,1) / (phi(:,1) + pe2)**2)
			der = atmosphere%Pg(1) + (atmosphere%Pg(1) - pe2)**2 * der
			
			call gasPressure(pe2, atmosphere%T(1), PH, PHminus, PHplus, PH2, PH2plus, PTotal)
			k_c = backgroundOpacityPerHParticle(atmosphere%T(1), pe2, PH, PHminus, PHplus, PH2, PH2plus, 5000.d0) / mu
			
			dy(1) = (GRAV * 10.d0**atmosphere%lTau500Ordered(1)/k_c+dpdgt)/der
			pe0 = atmosphere%Pe(1)
			if (dy(1) <= 0.d0) then
				atmosphere%Pe(1) = pe0 / 2.d0
			else
				atmosphere%Pe(1) = (dy(1)+1.d0)*pe2
			endif
			delta = abs(atmosphere%Pe(1)/pe0-1.d0)
			loop = loop + 1
			if (loop == 20) then
				delta = 0.d0
			endif
		enddo
								
		y(1) = log10(atmosphere%Pe(1))
		if (loop == 20) then
			! print *, 'Hydrostatic equilibrium: Pe at boundary might be wrong'
			stop
		endif
		
		do i = 2, atmosphere%nDepths
			y(i) = y(i-1) + (atmosphere%lTau500Ordered(i) - atmosphere%lTau500Ordered(i-1)) * dy(i-1)
			atmosphere%Pe(i) = 10.d0**y(i)
			atmosphere%Pg(i) = sum(abundancesRatio * phi(:,i) / (phi(:,i) + atmosphere%Pe(i)))
			atmosphere%Pg(i) = atmosphere%Pe(i) * (1.d0 + 1.d0 / atmosphere%Pg(i))
			dpdgt = sum(abundancesRatio * tdphi(:,i) / (phi(:,i) + atmosphere%Pe(i))**2)
			dpdgt = (atmosphere%Pg(i) - atmosphere%Pe(i))**2 * dpdgt * dtdx(i)
			der = sum(abundancesRatio * phi(:,i) / (phi(:,i) + atmosphere%Pe(i))**2)			
			der = atmosphere%Pg(i) + (atmosphere%Pg(i) - atmosphere%Pe(i))**2 * der

			call gasPressure(atmosphere%Pe(i), atmosphere%T(i), PH, PHminus, PHplus, PH2, PH2plus, PTotal)
			k_c = backgroundOpacityPerHParticle(atmosphere%T(i), atmosphere%Pe(i), PH, PHminus, PHplus, PH2, PH2plus, 5000.d0) / mu
						
			dy(i) = (GRAV * 10.d0**atmosphere%lTau500Ordered(i)/k_c+dpdgt)/der
		enddo
				
		atmosphere%Pe = atmosphere%Pe(atmosphere%order)		
		atmosphere%Pg = atmosphere%Pg(atmosphere%order)			

	end subroutine

end module hydrostaticModule