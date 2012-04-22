FC = gfortran -ffree-line-length-none

all: clean
	$(FC) -c initialization.f90 -g
	$(FC) -c staticfields.f90 loop.f90 -g
	$(FC) main.f90 initialization.f90 staticfields.f90 loop.f90 -o surya -g

clean:
	rm -f surya diffrot.dat psi.dat final.dat run.dat ber.dat butbot.dat rad.dat
	rm -f *~

refresh: clean
	rm -f *.o *.mod
