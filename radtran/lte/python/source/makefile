COMPILER = gfortran
FLAGS = 
OPTIONS = -c -march=native -ffree-line-length-none -fPIC -O3
 
all: 
	make pylte
	
clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find ../ -maxdepth 2 -name "pylte.so" -delete

pylte: vars.o maths.o partition.o background.o lte.o synth.o hydros.o

vars.o: vars.f90
	$(COMPILER) $(OPTIONS)  vars.f90

maths.o: maths.f90 vars.o partition.o
	$(COMPILER) $(OPTIONS)  maths.f90
	
lte.o: lte.f90 vars.o maths.o hydros.o partition.o synth.o
	$(COMPILER) $(OPTIONS)  lte.f90
	
partition.o: partition.f90
	$(COMPILER) $(OPTIONS)  partition.f90
		
background.o: background.f90 maths.o
	$(COMPILER) $(OPTIONS)  background.f90
			
hydros.o: hydros.f90 vars.o partition.o maths.o background.o
	$(COMPILER) $(OPTIONS)  hydros.f90
	
synth.o: synth.f90 partition.o vars.o background.o hydros.o
	$(COMPILER) $(OPTIONS)  synth.f90
	
