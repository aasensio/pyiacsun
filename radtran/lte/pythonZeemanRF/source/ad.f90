module dualClass
implicit none

	type, public :: dual	
		real(kind=8) :: rp
		real(kind=8), dimension(64) :: ip   ! Real and infinitesimal part
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

	interface assignment(=)
		module procedure equalDual, equalReal
	end interface

	interface operator(==)
		module procedure equalToDual, equalToReal
	end interface

	interface operator(/=)
		module procedure notequalToDual, notequalToReal
	end interface

	interface operator(>)
		module procedure greaterThanDual, greaterThanReal
	end interface

	interface operator(>=)
		module procedure greaterEqualThanDual, greaterEqualThanReal
	end interface

	interface operator(<)
		module procedure smallerThanDual, smallerThanReal
	end interface

	interface operator(<=)
		module procedure smallerEqualThanDual, smallerEqualThanReal
	end interface

	interface maxval
		module procedure maxvalDual
	end interface

	! interface matmul
	! 	module procedure matmulDual
	! end interface
	
	interface exp
		module procedure expDual
	end interface

	interface abs
		module procedure absDual
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
		c%rp = a * b%rp
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
		c%rp = a + b%rp
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
		c%rp = a - b%rp
		c%ip = -b%ip
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
! Equal
!---------------
	elemental subroutine equalDual(b,a)
		type (dual), intent(in) :: a
		type (dual), intent(out) :: b			
		b%rp = a%rp
		b%ip = a%ip
	end subroutine equalDual

	elemental subroutine equalReal(b, a)
		real(kind=8), intent(in) :: a		
		type (dual), intent(out) :: b
		b%rp = a
		b%ip = 0.d0
	end subroutine equalReal

!---------------
! Equal to
!---------------
	elemental function equalToDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		logical :: c
		if (a%rp == b%rp) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function equalToDual

	elemental function equalToReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		logical :: c
		if (a%rp == b) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function equalToReal

!---------------
! Equal to
!---------------
	elemental function notequalToDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		logical :: c
		if (a%rp /= b%rp) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function notequalToDual

	elemental function notequalToReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		logical :: c
		if (a%rp /= b) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function notequalToReal

!---------------
! Greater than
!---------------
	elemental function greaterThanDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		logical :: c
		if (a%rp > b%rp) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function greaterThanDual

	elemental function greaterThanReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		logical :: c
		if (a%rp > b) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function greaterThanReal

!---------------
! Smaller than
!---------------
	elemental function smallerThanDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		logical :: c
		if (a%rp < b%rp) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function smallerThanDual

	elemental function smallerThanReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		logical :: c
		if (a%rp < b) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function smallerThanReal

!---------------
! Greater Equal than
!---------------
	elemental function greaterEqualThanDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		logical :: c
		if (a%rp >= b%rp) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function greaterEqualThanDual

	elemental function greaterEqualThanReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		logical :: c
		if (a%rp >= b) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function greaterEqualThanReal

!---------------
! Smaller Equal than
!---------------
	elemental function smallerEqualThanDual(a,b) result(c)
		type (dual), intent(in) :: a, b
		logical :: c
		if (a%rp <= b%rp) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function smallerEqualThanDual

	elemental function smallerEqualThanReal(a,b) result(c)
		type (dual), intent(in) :: a
		real(kind=8), intent(in) :: b
		logical :: c
		if (a%rp <= b) then
			c = .TRUE.
		else
			c = .FALSE.
		endif
	end function smallerEqualThanReal
	
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
! ABS
!---------------
	elemental function absDual(a) result(c)
		type (dual), intent(in) :: a
		type (dual) :: c
		
		c%rp = abs(a%rp)
		c%ip = a%rp / abs(a%rp) * a%ip

	end function absDual

!---------------
! MAXVAL
!---------------
	function maxvalDual(a) result(c)
		type (dual), intent(in) :: a(:)
		type (dual) :: c
		integer :: ind(1)
		
		ind = maxloc(a%rp)

		c%rp = a(ind(1))%rp
		c%ip = a(ind(1))%ip

	end function maxvalDual

!---------------
! MAXVAL
!---------------
	! function matmulDual(a, b) result(c)
	! 	type (dual), intent(in) :: a(:,:), b(:,:)
	! 	type (dual) :: c(:,:)
				
	! 	! c%rp = matmul(a%rp,b%rp)
	! 	! c%ip = matmul(a%ip,b%rp)

	! end function matmulDual
			
end module dualClass