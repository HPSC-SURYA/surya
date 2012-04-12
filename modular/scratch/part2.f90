!
!	PART II.  This is the part where the alpha coefficient, 
!	diffusivity, differential rotation and meridional circulation 
!	are specified. Level II Users wishing to make changes in 
!	this part should read Sect. 4 of the 'Guide'.
!

module part2

contains

	subroutine define_profiles(al, eta, etab, dom)
		use params
		use part1
	
		implicit none
		
		! general namelist declarations
		
		
		! part 1 namelist declarations
		double precision ::	irelax=0  ! irelax=0 means code-init, irelax=1 means use init.dat
		double precision ::	tmax=1.0d0 ! max time of computation (1e8 s)
		double precision ::	v0=-29.0d0 ! amplitude of meridional flow (m/s)
		double precision ::	et0=2.6d0 ! poloidal diffusivity (1e12 cm2/s)
		double precision ::	et1=0.04d0 ! toroidal diffusivity (1e12 cm2/s)
		double precision ::	al0=25.0d0 ! alpha-effect amplitude (m/s)
		double precision ::	dt=0.0002d0 ! time step (1e8 s) ! NOTE: Must satisfy the Courant condition (explicit?)
		
		! part2 namelist declarations		
		double precision :: 	pm ! outer radius (1e10 cm)
		double precision ::  lbf ! lower boundary fraction (of the radius)
		double precision :: 	pb ! inner radius (1e10 cm)
		double precision :: 	pw ! ???
		double precision :: 	qm ! range of theta
		double precision ::	dp=(pm-pb)/float(n) ! delta_radius
		double precision ::	dq=-qm/float(n) ! delta_theta
		
		double precision :: ra(nmax)
		
		double precision, intent(out) :: al(nmax,nmax), eta(nmax,nmax), etab(nmax,nmax), dom(nmax,nmax)
		
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		namelist /general/ nmax, lmax, n, dr_filename, pi
		namelist /part1/ irelax, tmax, v0, et0, et1, al0, dt
		namelist /part2/ pm, lbf, pw, qm
			pb=lbf*pm ! inner radius (1e10 cm)
			qm=pi ! range of theta
			dp=(pm-pb)/float(n) ! delta_radius
			dq=-qm/float(n) ! delta_theta
		
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		open(25,file=dr_filename,status='unknown',access='append')   ! dr_filename defined in module params

		do i = 1, nmax
		do j = 1, nmax

			ra(i)=pb+float(i-1)/float(n)*(pm-pb) 
			p=pb+float(i-1)/float(n)*(pm-pb)
			q=qm-float(j-1)/float(n)*qm
			co=dcos(q)

			! ALPHA PROFILE
				al(i,j) = 0.25d0*al0*co*(1.00d0+erf((p-0.95d0*pm)/(0.025d0*pm)))*(1.00d0-erf((p-pm)/(0.025d0*pm)))

			! DIFFUSIVITY PROFILES
		       eta(i,j) = 0.000220d0 + (et0/2.00d0)*(1.0d0 + erf((p-0.7d0*pm)/(0.025d0*pm)))
		       etab(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 + erf((p-0.72d0*pm)/(0.025d0*pm)))+ (et0/2.0d0)*(1.0d0+erf((p-0.975d0*pm)/(0.025d0*pm))) 

			! SOLAR DIFFERENTIAL ROTATION PROFILE
		       dom(i,j) = 271.9d0 + 0.50d0*(1.000d0 + erf((p-0.7d0*pm)/(0.025d0*pm)))*(289.5d0 - 39.4d0*co*co - 42.2d0*co*co*co*co - 271.9d0)         
		       write(25,'3(f13.5,1x)') q,p,dom(i,j)

		end do
		end do
	
		close(25)	
	
	end subroutine define_profiles
	
	
	
	subroutine define_profiles_midpoint()
		!	The do loops are closed.  We also need the diffusivity at
		!	mid-points in the grid for some calculations.  This is 
		!	obtained and stored now.
		
		

		do i=1,2*n+1
		do j=1,2*n+1
			 p=pb+float(i-1)*(pm-pb)/float(2*n)
			 etab2(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 + erf((p-0.72d0*pm)/(0.025d0*pm)))+ (et0/2.0d0)*(1.0d0+erf((p-0.975d0*pm)/(0.025d0*pm))) 
		end do
		end do
	
	end subroutine define_profiles_midpoint

end module part2
