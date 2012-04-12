! here are the routines contained in the main loop of the program

module loop

   use initialization
   use standalone

   implicit none

!      double precision :: mytimevalue

contains

   subroutine advance_interior()
   ! part V-A

      double precision :: br, bt

      ! solve explicit side of A evolution

         do i=2,n
            do j=2,n
               phi(i,j)=-d(i,j)*u(i-1,j)+(-e(i,j)+1.0d0)*u(i,j)-f(i,j)*u(i+1,j)+al(i,j)*ub(i,j)*dt/2.0d0
            end do
         end do

      ! do implicit matrix solve for A along constant r, for all r

         do i=2,n
            do j=2,n
               a1(j-1)=a(i,j)
               b1(j-1)=b(i,j)+1.0d0
               c1(j-1)=c(i,j)
               r(j-1)=phi(i,j)
            end do
            r(1)=phi(i,2)-a1(1)*u(i,1)
            r(n-1)=phi(i,n)-c1(n-1)*u(i,n+1)
            call tridag(a1,b1,c1,r,uu,n-1)
            do j=2,n
               u(i,j) = uu(j-1)
               phi(i,j)=-phi(i,j)+2.0d0*uu(j-1)
            end do
         end do

      ! solve explicit side of B evolution

         do i=2,n
            do j=2,n
               p=pb+float(i-1)/float(n)*(pm-pb)
               q=qm-float(j-1)/float(n)*qm
               br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt/(4.0d0*dq)
               bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/(4.0d0*p*dp)
               phib(i,j)=-db(i,j)*ub(i-1,j)+(-eb(i,j)+1.0d0)*ub(i,j)-fb(i,j)*ub(i+1,j)+br*dror(i,j)+bt*drot(i,j)
            end do
         end do

      ! do implicit solve for B along constant r, for all r

         do i=2,n
            do j=2,n
               a1(j-1)=ab(i,j)
               b1(j-1)=bb(i,j)+1.0d0
               c1(j-1)=cb(i,j)
               r(j-1)=phib(i,j)
            end do
            r(1)=phib(i,2)-a1(1)*ub(i,1)
            r(n-1)=phib(i,n)-c1(n-1)*ub(i,n+1)
            call tridag(a1,b1,c1,r,uub,n-1)
            do j=2,n
               ub(i,j) = uub(j-1)
               phib(i,j)=-phib(i,j)+2.0d0*uub(j-1)
               phi(i,j)=phi(i,j)+al(i,j)*ub(i,j)*dt/2.0d0
            end do
         end do
         


      !  NOW FOR OTHER DIRECTION

      ! do implicit solve for A along constant theta, for all theta

         do j=2,n
            do i=2,n
               a1(i-1)=d(i,j)
               b1(i-1)=e(i,j)+1.0d0
               c1(i-1)=f(i,j)
               r(i-1)=phi(i,j)
            end do
            r(1)=phi(2,j)-a1(1)*u(1,j)
            r(n-1)=phi(n,j)-c1(n-1)*u(n+1,j)
            call tridag(a1,b1,c1,r,uu,n-1)
            do i=2,n
               u(i,j)=uu(i-1)
               p=pb+float(i-1)/float(n)*(pm-pb)
               q=qm-float(j-1)/float(n)*qm
               br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt/(4.0d0*dq)
               bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/(4.0d0*p*dp)
               phib(i,j)=phib(i,j)+br*dror(i,j)+bt*drot(i,j)
            end do
         end do

      ! do implicit solve for B along constant theta, for all theta

         do j=2,n
            do i=2,n
               a1(i-1)=db(i,j)
               b1(i-1)=eb(i,j)+1.0d0
               c1(i-1)=fb(i,j)
               r(i-1)=phib(i,j)
            end do
            r(1)=phib(2,j)-a1(1)*ub(1,j)
            r(n-1)=phib(n,j)-c1(n-1)*ub(n+1,j)
            call tridag(a1,b1,c1,r,uub,n-1)
            do i=2,n
               ub(i,j)=uub(i-1)
            end do
         end do

   end subroutine advance_interior


   subroutine advance_boundaries()
   ! part V-B

      integer :: nup
      double precision :: cf

      ! Boundary condition at bottom

         do j=1,n+1
            u(1,j)=0.0d0
            ub(1,j)=0.0d0
         end do

      ! Boundary condition at top surface

         do nup=1,11
            do l=1,96
               ss1(l)=ss1_p(l)*ss(n,dq)
            end do
            do j=1,n+1
               call coeff(n,ss1,dq,j,cf)
               u(n+1,j)=u(n,j)+cf*dp
            end do
         end do

         do j=1,n+1
            ub(n+1,j)=0.0d0
         end do

      ! Boundary conditions at the poles

         do i=2,n+1
            u(i,n+1)=0.0d0
            ub(i,n+1)=0.0d0
            u(i,1)=0.0d0
            ub(i,1)=0.0d0
         end do

   end subroutine advance_boundaries


   subroutine magnetic_buoyancy()
   ! part V-C

      double precision :: ber, bold, tau, qjer
      integer :: ier, jer, nmb, kmb, kch, n2, nt

      tau = .0088d0
      nmb = tau/dt
      kmb = k/nmb
      kch = nmb*kmb	! this is like a ghetto modulus?

      if (k.eq.kch) then

         n2 = n/2.
         nt = n - 6.
	
         do i=1,n2
         do j=2,n
            if (ra(i).gt.0.71*pm) then
               if (ABS(ub(i,j)).GT.0.8d0) then
                jer=j
                ier=i
                qjer=qm-float(jer-1)*qm/float(n)	
                ber=ub(ier,jer)
                bold=ub(nt,jer)
                ub(nt,jer)=((ra(ier)/ra(nt))*0.5d0*ber) + bold
                ub(ier,jer)=0.5d0*ber
                 open(25,file='ber.dat',status='unknown',access='append')
                 write(25,'(d15.8,1x,d15.8,1x,d15.8)') mytimevalue, qjer, ber
                 close(25)
               end if
            end if
         end do
         end do

      end if

   end subroutine magnetic_buoyancy


   subroutine record_data()
   ! part V-D

      double precision :: br, p20, tdiff
      integer :: istep

      ! WRITING IN THE FILES 'rad.dat' and 'butbot.dat' after certain intervals of time

      if (mytimevalue.gt.0.0d0) then

         p20 = pi/20.
         istep = mytimevalue/p20
         tdiff = mytimevalue - float(istep)*p20

         if(.not.(tdiff.ge.dt)) then

            open(17,file=radial_filename,status='unknown',access='append') 
            open(18,file=bottom_filename,status='unknown',access='append') 
            do j = 2, n
               q=qm-float(j-1)/float(n)*qm
            br = -(u(n-1,j-1)*dsin(q-dq) - u(n-1,j+1)*dsin(q+dq))/(2.0d0*dq*dsin(q))/pm
            write(17,'(f13.7,1x,f13.7,1x,f13.7)') mytimevalue, q, br
            write(18,'(f13.7,1x,f13.7,1x,f13.7)') mytimevalue, q, ub(45,j)
            end do
            close(17)
            close(18)
            close(19)

         end if

      end if	

      ! WRITING IN THE FILE 'run.dat' after every 40 time steps

      if (k/40*40.eq.k) then

         open(95,file=run_filename,status='unknown',access='append')
!         write(95,'(d15.9,1x,d15.9,1x,d15.9,1x,d15.9,1x,d15.9,1x,d15.9)') t, ub(nmax/2,nmax/2), ub(nmax/2+1,nmax/2+1), ub(nmax/2+2,nmax/2+2), u(nmax/2,nmax/2), u(nmax/2+1,nmax/2+1)
         write(95,'(d15.9,1x,d15.9,1x,d15.9,1x,d15.9,1x,d15.9,1x,d15.9)') mytimevalue, ub(2,3), ub(3,2), ub(4,5), u(2,3), u(3,2)
         close(95)

!         print*, 'nmax = ', nmax
!         print*, 'dimensions = ', shape(u)
!         print*, 'dimensions = ', shape(ub)
!         print*, 't is ', t
 !        print*, 'u(2,3) is ', u(2,3)
 !        print*, 'ub(2,3) is ', ub(2,3)

         stop 'checkpoint reached'

      end if

   end subroutine record_data

end module loop
