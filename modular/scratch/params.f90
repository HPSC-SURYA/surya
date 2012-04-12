! Overall parameter list for the code.

module params

	! problem size
	integer, parameter :: nmax=257,lmax=96
	integer, parameter :: n=nmax-1
	
	! output file names
	character(len=*) :: dr_filename='diffrot.dat'

	! constants
	double precision, parameter :: pi=4.0d0*atan(1.0d0)

end module params
