module globalModule
implicit none

	real(kind=8), parameter :: PK = 1.3806503d-16, UMA = 1.66053873d-24, PC = 2.99792458d10
	real(kind=8), parameter :: PH = 6.62606876d-27, PHK = PH / PK, PHC = 2.d0 * PH / PC**2
	real(kind=8), parameter :: PME = 9.10938188d-28, PE = 4.8032d-10, PI = 3.14159265359d0
	real(kind=8), parameter :: PHC2 = 2.d0 * PH * PC**2
	real(kind=8), parameter :: SQRTPI = 1.77245385091d0, EARTH_SUN_DISTANCE = 1.495979d13
	real(kind=8), parameter :: LARMOR = PE/(4.d0*PI*PME*PC), PMP = 1.67262158d-24
	real(kind=8), parameter :: NUCLEAR_MAGNETON = PE*PH/(2.d0*PI*PMP)
	real(kind=8), parameter :: BOHR_MAGNETON = PE*PH/(2.d0*PI*PME)
	real(kind=8), parameter :: PH2C3 = 2.d0 * PH**2 * PC**3, PERTURBATION = 0.005d0
	real(kind=8), parameter :: OPA = PI * PE**2 / (PME * PC), EV_ERG = 1.60217646d-12, GRAV = 2.7398d4
	real(kind=8), parameter :: BOHR_R = (PH/(2.d0*PI))**2/(PME*PE**2), CEV_ERG = 1.60217646d-12


	type atmosphereType		
		character(len=256) :: file
		integer :: nDepths, nNodes
		real(kind=8) :: abundance
		real(kind=8), dimension(:), pointer :: lTau500, height, T, nodes, niovern, ui, u1, u2, u3, n2overn1, n1overn0, microturbulence, vmac
		real(kind=8), dimension(:), pointer :: opacity500
 		real(kind=8), dimension(:), pointer :: Pe, Pg, PH, PHminus, PHplus, PH2, PH2plus, PTotal, nHtot, TOrdered, lTau500Ordered
 		integer, dimension(:), pointer :: order 		
	end type atmosphereType

	type transitionType
		integer :: nLambda, element, ionization
		logical :: active
		real(kind=8) :: lambda0, gf, Elow, gu, Ju, gl, Jl, sigmaABO, alphaABO, lambdaLeft, lambdaRight
		real(kind=8) :: frequency0, boundary
		real(kind=8), dimension(:), pointer :: lambda, intensity, observed, continuumIntensity
		real(kind=8), dimension(:), pointer :: lineOpacity, backOpacity, dopplerWidth, deltaNu, damping, frequency
		real(kind=8), dimension(:,:), pointer :: opacity, source, opacityContinuum
	end type transitionType
	
	type lineListType
		character(len=256) :: file
		integer :: nLines, nLambdaTotal
		type(transitionType), dimension(:), pointer :: transition
		real(kind=8), dimension(:,:), pointer :: opacity, source, opacityContinuum
		real(kind=8), dimension(:), pointer :: intensity, boundary, intensityContinuum
	end type lineListType
	
	
	type atlasType
		character(len=256) :: file
		integer :: nLambda
		real(kind=8), dimension(:), pointer :: lambda, intensity
	end type atlasType
	
	type(atmosphereType) :: atmosphere
	type(lineListType) :: lineList
	type(atlasType) :: atlas

end module globalModule