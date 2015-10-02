module dualClass
use adMaths
implicit none

	type, public :: dual	
		real(kind=8) :: rp
		real(kind=8), dimension(dualSize) :: ip   ! Real and infinitesimal part
	contains
		procedure, pass :: init
	end type dual
	
! Operator overload
	interface operator(+)
		module procedure additionDualDual, additionDualReal, additionRealDual
	end interface
	
	interface operator(-)
		module procedure substractionDualDual, substractionDualReal, substractionRealDual, signChangeDual
	end interface
	
	interface operator(*)
		module procedure productDualDual, productDualReal, productRealDual
	end interface
	
	interface operator(/)
		module procedure divisionDualDual, divisionDualReal, divisionRealDual
	end interface
	
	interface operator(**)
		module procedure powerDual, powerIntDual
	end interface
	
	interface exp
		module procedure expDual
	end interface
		
	interface sin
		module procedure sinDual
	end interface
	
	interface cos
		module procedure cosDual
	end interface
	
	interface log
		module procedure logDual
	end interface
	
	interface tanh
		module procedure tanhDual
	end interface
	
	interface sqrt
		module procedure sqrtDual
	end interface
	
	interface sum
		module procedure sumDual
	end interface
	
	interface lnGamma
		module procedure lnGammaDual
	end interface
	
	interface voigt
		module procedure voigtDual
	end interface
	
	contains

!---------------------------
! Creation of the dual number
! value and whichIndependent are optional arguments
! value is used to initialize the value of the number
! whichIndependent is used to set independent variables. If not present, the variable is dependent
!---------------------------
	elemental subroutine init(this, value, whichIndependent)
	class(dual), intent(inout) :: this
	real(kind=8), intent(in), optional :: value
	integer, intent(in), optional :: whichIndependent
				
		this%ip = 0.d0
		
		this%rp = 0.d0
		if (present(value)) then
			this%rp = value
		endif
		
		this%ip = 0.d0
		if (present(whichIndependent)) then
			this%ip(whichIndependent) = 1.d0
		endif
		
	end subroutine init
	
! Functions for defining the basic operations

! ---------------
! MULTIPLICATION
!----------------
	elemental function productDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp * b%rp
		c%ip = a%rp * b%ip + b%rp * a%ip
	end function productDualDual
	
	elemental function productDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp * b
		c%ip = a%ip * b
	end function productDualReal
	
	elemental function productRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp * a
		c%ip = b%ip * a
	end function productRealDual

! ---------------
! DIVISION
!----------------
	elemental function divisionDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp / b%rp
		c%ip = (b%rp * a%ip - a%rp * b%ip) / b%rp**2
	end function divisionDualDual
	
	elemental function divisionDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp / b
		c%ip = a%ip / b
	end function divisionDualReal
	
	elemental function divisionRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = a / b%rp
		c%ip = -a * b%ip / b%rp**2
	end function divisionRealDual

!---------------
! ADDITION
!---------------
	elemental function additionDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp + b%rp
		c%ip = a%ip + b%ip
	end function additionDualDual
		
	elemental function additionDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp + b
		c%ip = a%ip
	end function additionDualReal
	
	elemental function additionRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp + a
		c%ip = b%ip
	end function additionRealDual

!---------------
! SUBSTRACTION
!---------------
	elemental function substractionDualDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		type (dual) :: c
		c%rp = a%rp - b%rp
		c%ip = a%ip - b%ip
	end function substractionDualDual
	
	elemental function substractionDualReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c
		c%rp = a%rp - b
		c%ip = a%ip
	end function substractionDualReal
	
	elemental function substractionRealDual(a,b) result(c)
		real(kind=8), intent(in) :: a
		type (dual), intent(in) :: b		
		type (dual) :: c
		c%rp = b%rp - a
		c%ip = b%ip
	end function substractionRealDual
	
	elemental function signChangeDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = -a%rp
		c%ip = -a%ip
	end function signChangeDual

!---------------
! POWER
!---------------	
	elemental function powerDual(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		type (dual) :: c		
		c%rp = a%rp**b
		c%ip = b * a%rp**(b-1.d0) * a%ip		
	end function powerDual
	
	elemental function powerIntDual(a,b) result(c)
		type (dual), intent(in) :: a
		integer, intent(in) :: b
		type (dual) :: c		
		c%rp = a%rp**b
		c%ip = b * a%rp**(b-1.d0) * a%ip		
	end function powerIntDual
	
!---------------
! EXP
!---------------
	elemental function expDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c				
		c%rp = exp(a%rp)
		c%ip = a%ip * exp(a%rp)
	end function expDual

!---------------
! SIN
!---------------
	elemental function sinDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = sin(a%rp)
		c%ip = a%ip * cos(a%rp)
	end function sinDual

!---------------
! COS
!---------------
	elemental function cosDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = cos(a%rp)
		c%ip = -a%ip * sin(a%rp)
	end function cosDual

!---------------
! LN
!---------------
	elemental function logDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		c%rp = log(a%rp)
		c%ip = a%ip / a%rp
	end function logDual
		
!---------------
! TANH
!---------------
	elemental function tanhDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		real(kind=8) :: t
		c%rp = tanh(a%rp)
		c%ip = 1.d0/cosh(a%rp)**2 * a%ip
	end function tanhDual
	
!---------------
! SQRT
!---------------
	elemental function sqrtDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		
		c%rp = sqrt(a%rp)
		c%ip = 0.5d0/sqrt(a%rp) * a%ip
	end function sqrtDual
	
!---------------
! ADDITION OF A VECTOR
!---------------
	function sumDual(a) result(c)
		type (dual), intent(in), dimension(:) :: a
		type (dual), dimension(size(a)) :: c
		integer :: i
		c%rp = sum(a%rp)
 		do i = 1, size(a(1)%ip)
 			c%ip(i) = sum(a%ip(i))
 		enddo
	end function sumDual

!---------------
! lnGamma
!---------------
	function lnGammaDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		integer :: ierr		
		
		c%rp = alngam(a%rp, ierr)
		c%ip = digama(a%rp, ierr) * a%ip
		
	end function lnGammaDual
	
!---------------
! Voigt profile
!---------------
	function voigtDual(a, v) result(c)
		type (dual), intent(in) :: a, v
		type (dual) :: c
		integer :: ierr
		real(kind=8) :: H(2), ga, gv, PI = 3.14159265359d0
		
			H = voigtAutomaticDiff(a%rp, v%rp)
			ga = 2.d0 * (-1.d0 / sqrt(PI) + a%rp * H(1) + v%rp * H(2))
			gv = 2.d0 * (-v%rp * H(1) + a%rp * H(2))
						
			c%rp = H(1)
			c%ip = ga * a%ip + gv * v%ip
		
	end function voigtDual
	
end module dualClass