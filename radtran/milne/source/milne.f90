module milneMod
use vars
use atomic_functions
use math_functions, only : init_maths
use iso_c_binding, only: c_int, c_double
implicit none

contains

	subroutine c_setline(nLambda, lineData, waveOut) bind(c)
	integer(c_int), intent(in) :: nLambda
	real(c_double), intent(in) :: lineData(7)
	real(c_double), intent(out), dimension(nLambda) :: waveOut
	integer :: k
				
		line(1)%wave0 = lineData(1)
		line(1)%Jup = lineData(2)
		line(1)%Jlow = lineData(3)
		line(1)%gup = lineData(4)
		line(1)%glow = lineData(5)
		line(1)%lambdaInit = lineData(6)
		line(1)%lambdaStep = lineData(7)
		line(1)%nLambda = nLambda
		
		Stokes_Syn%nlambda = line(1)%nLambda
		
		if (associated(Stokes_Syn%lambda)) deallocate(Stokes_Syn%lambda)
		if (associated(Stokes_Syn%stokes)) deallocate(Stokes_Syn%stokes)
		allocate(Stokes_Syn%lambda(line(1)%nLambda))
		allocate(Stokes_Syn%stokes(4,line(1)%nLambda))
		
		Stokes_Syn%stokes = 0.d0
		
		do k = 1, line(1)%nLambda
			Stokes_Syn%lambda(k) = line(1)%lambdaInit + line(1)%lambdaStep * (k-1)						
		enddo
		
		waveOut = Stokes_Syn%lambda
		
		call init_maths
		
		if (allocated(zeeman_voigt)) deallocate(zeeman_voigt)
		if (allocated(zeeman_faraday)) deallocate(zeeman_faraday)
		
		if (allocated(ki)) deallocate(ki)
		if (allocated(kq)) deallocate(kq)
		if (allocated(ku)) deallocate(ku)
		if (allocated(kv)) deallocate(kv)
		if (allocated(fq)) deallocate(fq)
		if (allocated(fu)) deallocate(fu)
		if (allocated(fv)) deallocate(fv)
		
		if (allocated(delta)) deallocate(delta)
		if (allocated(stokes)) deallocate(stokes)
		
		if (allocated(profile)) deallocate(profile)
		if (allocated(v)) deallocate(v)
		
		allocate(profile(2,line(1)%nLambda))
		allocate(v(line(1)%nLambda))
				
		allocate(zeeman_voigt(3,line(1)%nLambda))
		allocate(zeeman_faraday(3,line(1)%nLambda))
		
		allocate(ki(line(1)%nLambda))
		allocate(kq(line(1)%nLambda))
		allocate(ku(line(1)%nLambda))
		allocate(kv(line(1)%nLambda))
		allocate(fq(line(1)%nLambda))
		allocate(fu(line(1)%nLambda))
		allocate(fv(line(1)%nLambda))
		
		allocate(delta(line(1)%nLambda))
		allocate(stokes(4,line(1)%nLambda))

		nLines = 1
		
	end subroutine c_setline

	subroutine c_addline(lineData) bind(c)
	real(c_double), intent(in) :: lineData(7)
	integer :: k

		nLines = nLines + 1
				
		line(nLines)%wave0 = lineData(1)
		line(nLines)%Jup = lineData(2)
		line(nLines)%Jlow = lineData(3)
		line(nLines)%gup = lineData(4)
		line(nLines)%glow = lineData(5)
		
	end subroutine c_addline
	
	subroutine c_milnesynth(nLambda, modelIn, muIn, stokesOut) bind(c)
	integer(c_int), intent(in) :: nLambda
	real(c_double), intent(in) :: modelIn(9), muIn
	real(c_double), intent(out), dimension(4,nLambda) :: stokesOut
	
 		model%Bfield = modelIn(1)
 		model%theta = modelIn(2)
 		model%chi = modelIn(3)
 		model%vmac = modelIn(4)
 		model%damping = modelIn(5)
 		model%B0 = modelIn(6)
 		model%B1 = modelIn(7)
 		model%doppler = modelIn(8)
 		model%kl = modelIn(9)
 		
 		model%mu = muIn
 				
 		call synthesize(model,line,Stokes_Syn)
 		 		
 		StokesOut = Stokes_Syn%stokes
 		
	end subroutine c_milnesynth

	subroutine c_milnesynthmany(nLambda, nModels, modelIn, muIn, stokesOut) bind(c)
	integer(c_int), intent(in) :: nLambda, nModels
	real(c_double) :: muIn
	real(c_double), intent(in), dimension(nModels,9) :: modelIn
	real(c_double), intent(out), dimension(4,nLambda,nModels) :: stokesOut
	integer(c_int) :: i
	
		do i = 1, nModels
			model%Bfield = modelIn(i,1)
			model%theta = modelIn(i,2)
			model%chi = modelIn(i,3)
			model%vmac = modelIn(i,4)
			model%damping = modelIn(i,5)
			model%B0 = modelIn(i,6)
			model%B1 = modelIn(i,7)
			model%doppler = modelIn(i,8)
			model%kl = modelIn(i,9)
			model%mu = muIn
					
			call synthesize(model,line,Stokes_Syn)
					
			StokesOut(:,:,i) = Stokes_Syn%stokes
		enddo
 		
	end subroutine c_milnesynthmany

end module milneMod