program test
use lteMod
implicit none
real(kind=8) :: lineListIn(11,3), lambdaAxis(100), atmosphereIn(7,64), stokesOut(5,100), dStokesOut(5,100,64)
integer :: i, j, nLambda, nDepths

	nDepths = 64

	call c_initatmosphere(nDepths)

	open(unit=12,file='lines.dat',action='read',status='old')
	do i = 1, 3
		read(12,*) (lineListIn(j,i),j=1,11)
	enddo
	close(12)

	open(unit=12,file='hsra_64.model',action='read',status='old')
	read(12,*)
	read(12,*)
	do i = 1, 64
		read(12,*) (atmosphereIn(j,i),j=1,7)
	enddo
	close(12)

	nLambda = 100
	do i = 1, 100
		lambdaAxis(i) = 10824 + 0.06*i
	enddo
	call c_initlines(3, lineListIn, nLambda, lambdaAxis)
	do i = 1, 5
		call c_synthlines(nDepths, atmosphereIn, nLambda, stokesOut, dStokesOut)
	enddo

end program test