all: clean
	gfortran -ffree-line-length-none -c standalone.f90 initialization.f90 loop.f90
	gfortran -ffree-line-length-none main.f90 initialization.f90 loop.f90 standalone.f90 -o surya

clean:
	rm -f surya
	rm -f *.dat
	rm -f *.o *.mod
	rm -f *~

