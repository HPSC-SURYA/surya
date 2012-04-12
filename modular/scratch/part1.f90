!
!	PART I. This is the ONLY part of the code in which Level I users
!	may want to make some changes.  Read Sect. 3 of the 'Guide' to
!	learn about the variables specified in this part.
!

module part1

	implicit none

	double precision ::	irelax=0  ! irelax=0 means code-init, irelax=1 means use init.dat
	double precision ::	tmax=1.0d0 ! max time of computation (1e8 s)
	double precision ::	v0=-29.0d0 ! amplitude of meridional flow (m/s)
	double precision ::	et0=2.6d0 ! poloidal diffusivity (1e12 cm2/s)
	double precision ::	et1=0.04d0 ! toroidal diffusivity (1e12 cm2/s)
	double precision ::	al0=25.0d0 ! alpha-effect amplitude (m/s)
	double precision ::	dt=0.0002d0 ! time step (1e8 s) ! NOTE: Must satisfy the Courant condition (explicit?)

end module part1
