program surya

   use initialization
   use loop
   use standalone
	
   implicit none

!      double precision :: mytimevalue

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialization
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

   ! read input namelists, allocate variables accordingly
   call read_input()
   call allocate_variables()

   mytimevalue=0.
   print*, 'first checkpoint: ', mytimevalue
  
   ! static field definitions
   call define_fields_grid()
   call define_fields_mid()
   call define_ddr()
   call define_mc_streamfunction()
   call define_mc_velocity()
  
   print*, 'second checkpoint: ', mytimevalue

   ! legendre polynomials
   call calculate_legendre()

   ! set initial conditions
   call initial_conditions()

   ! define matrix elements
   call define_matrix_elements()

   print*, 'initialization complete'
   

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Main Loop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   kend = tmax/dt
   mytimevalue=0.
   print*, 'mytimevalue set to zero:'
   print*, mytimevalue   ! this throws an error, why?
   stop 'checkpoint reached'
   do k=1,kend
      call advance_interior()
      call advance_boundaries()
      call magnetic_buoyancy()
      mytimevalue=mytimevalue+dt
      print*, 't is ', mytimevalue
      call record_data()
   end do


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Finalization
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   open(10,file=final_filename,status='unknown')
   do i=1,nmax
   do j=1,nmax
      p=pb+float(i-1)/float(n)*(pm-pb)
      q=qm-float(j-1)/float(n)*qm
      write(10,'(f13.7,1x,f13.7,1x,f13.7,1x,f13.7)') q, p, p*dsin(q)*u(i,j), ub(i,j)
   end do
   end do
   close(10)


end program surya
