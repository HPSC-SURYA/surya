! definitions of arrays and fields that don't change in the main loop (they're static)

module staticfields

   use initialization

   implicit none
   save

contains

   ! i think all these "field definitions" can probably be merged into one nested do-loop (not that it would speed things up much)

   subroutine define_fields_grid()
   ! define those fields with grid point resolution

      double precision :: co

	   open(25,file=dr_filename,status='unknown',access='append')   
	   do i = 1, nmax

         ! i-dependent stuff
         ra(i)=pb+float(i-1)/float(n)*(pm-pb)
         p=pb+float(i-1)/float(n)*(pm-pb)

	   do j = 1, nmax

         ! alpha, diffusivities, and differential rotation. DR is written out to dr_filename
         q=qm-float(j-1)/float(n)*qm
         co = dcos(q)
         al(i,j) = 0.25d0*al0*co*(1.00d0+erf((p-0.95d0*pm)/(0.025d0*pm)))*(1.00d0-erf((p-pm)/(0.025d0*pm)))                                     ! ALPHA PROFILE
         eta(i,j) = 0.000220d0 + (et0/2.00d0)*(1.0d0 + erf((p-0.7d0*pm)/(0.025d0*pm)))                                                          ! DIFFUSIVITY PROFILE
         etab(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 + erf((p-0.72d0*pm)/(0.025d0*pm)))+ (et0/2.0d0)*(1.0d0+erf((p-0.975d0*pm)/(0.025d0*pm)))   ! DIFFUSIVITY PROFILE
         dom(i,j) = 271.9d0 + 0.50d0*(1.000d0 + erf((p-0.7d0*pm)/(0.025d0*pm)))*(289.5d0 - 39.4d0*co*co - 42.2d0*co*co*co*co - 271.9d0)         ! SOLAR DIFFERENTIAL ROTATION PROFILE
         write(25,'(f13.5,1x,f13.5,1x,f13.5)') q, p, dom(i,j)

		end do
    	end do
    	close(25)

   end subroutine define_fields_grid


   subroutine define_fields_mid()
   ! define those fields with midpoint precision

      do i=1,2*n+1
      do j=1,2*n+1

         p=pb+float(i-1)*(pm-pb)/float(2*n)
         etab2(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 + erf((p-0.72d0*pm)/(0.025d0*pm)))+ (et0/2.0d0)*(1.0d0+erf((p-0.975d0*pm)/(0.025d0*pm)))

      end do
      end do

   end subroutine define_fields_mid

   subroutine define_ddr()
   ! define the derivatives of the DR profile

      do i=2,n
      do j=2,n

         dror(i,j)=(dom(i+1,j)-dom(i-1,j))/(2.0d0*dp)  
         drot(i,j)=(dom(i,j+1)-dom(i,j-1))/(2.0d0*dq)

      end do
      end do

   end subroutine define_ddr

   subroutine define_mc_streamfunction()
   ! define the streamfunction (psi) for the meridional circulation
   ! mid-point resolution

      double precision :: beta1, beta2, pp, del, gm, p0, exq0, exqm, gau

		beta1=1.5d0
		beta2=1.8d0
		pp=0.61d0*pm
		del=2.0000001d0
		gm=3.47d0
		p0=(pm-pb)/4.00d0     
   
      open(82,file=mc_filename,status='unknown',access='append')

      do i=1,2*n+1
      do j=1,2*n+1

         p=pb+float(i-1)/float(2*n)*(pm-pb)
         q=qm-float(j-1)/float(2*n)*qm

         if(q.le.pi/2.0d0) then
         exq0=dexp(-beta1*q**del)
         exqm=dexp(beta2*(q-pi/2.0d0))
         gau=dexp(-((p-p0)/gm)**2)
         psi(i,j)=(p-pp)*dsin(pi*(p-pp)/(pm-pp))*(1.00d0-exq0)*(1.00d0-exqm)*gau
         end if 

         if(q.gt.pi/2.0d0) then
         exq0=dexp(-beta1*(pi-q)**del)
         exqm=dexp(beta2*(pi/2.0d0-q))
         gau=dexp(-((p-p0)/gm)**2)
         psi(i,j)=-(p-pp)*dsin(pi*(p-pp)/(pm-pp))*(1.00d0-exq0)*(1.00d0-exqm)*gau
         end if

         if(p.le.pp)then
         psi(i,j)=0.0d0
         end if
         write(82,'(f13.5,1x,f13.5,1x,f13.5)') q, p, psi(i,j)

      end do
      end do

      close(82)

   end subroutine define_mc_streamfunction


   subroutine define_mc_velocity()
   ! calculate the MC velocities from the streamfunction
   ! mid-point resolution

      do i=2,2*n
      do j=2,2*n

         p=pb+float(i-1)*(pm-pb)/float(2*n)
         q=qm-float(j-1)*qm/float(2*n)
         vp1(i,j)=0.95d0*v0*(psi(i,j+1)-psi(i,j-1))/(p**2*dsin(q)*dq*(pm/p-0.95d0)**1.5)
         vq1(i,j)=-0.95d0*v0*(psi(i+1,j)-psi(i-1,j))/(p*dsin(q)*dp*(pm/p-0.95d0)**1.5)
         deta(i,j)=(etab2(i+1,j)-etab2(i-1,j))/(dp)
         vpb(i,j)=(vp1(i,j)-deta(i,j))*dt/(2.d0*p*dp)
         vp(i,j)=vp1(i,j)*dt/(2.d0*p*dp)
         vq(i,j)=vq1(i,j)*dt/(2.d0*p*dsin(q)*dq)

      end do
      end do

   end subroutine define_mc_velocity


   subroutine calculate_legendre()

      double precision :: cl, x

      do l=1,lmax

         cl=1.0d0/(1.0d0+float(l)/float(l+1)*(pm/pw)**(2*l+1))
         ss1_p(l)=float(2*l+1)/(float(l)*pm)*(cl-float(l)/float(l+1)*(1.0d0-cl))

      do j=1,n+1

         qm=4.0d0*atan(1.0d0)
         q=qm+float(j-1)*dq
         x=dcos(q)
         pl(j,l) = plgndr(l,1,x)
         sn(j) = dsin(q)

      end do
      end do

   end subroutine calculate_legendre


   subroutine initial_conditions()

      double precision :: ursintheta

      if (irelax.eq.0) then

         do i = 1,nmax
         do j = 1,nmax

            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
            u(i,j) = 0.0d0
            ub(i,j) = 1.0*dsin(1.0d0*q)*dsin(pi*((p-pb)/(pm-pb)))

         end do
         end do

      else

         open(12,file=init_filename,status='unknown')

         do i=1,n+1
         do j=1,n+1

            read(12,'(f13.7,1x,f13.7,1x,f13.7,1x,f13.7)') q, p, ursintheta, ub(i,j)
            if(j.eq.n+1) then
            ! boundary condition
               u(i,j)=0.0d0
            else
               u(i,j)=ursintheta/(p*dsin(q))
            end if

         end do
         end do

         close(12)

      end if

   end subroutine initial_conditions


   subroutine define_matrix_elements()
   ! pre-compute coefficients for solvers
   ! includes mesh geometry and velocity, since these never change

      do i=2,2*n
      do j=2,2*n

         dvp(i,j)=(vp1(i+1,j)-vp1(i-1,j))/(dp)

      end do
      end do

      do i=2,n
      do j=2,n

         p=pb+float(i-1)/float(n)*(pm-pb)
         q=qm-float(j-1)/float(n)*qm
         a(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)+vq(2*i-1,2*j-1)*dsin(q-dq/2.0d0)*(1.0d0+vq(2*i-1,2*j-2)*dsin(q-dq))/2.0d0)
         b(i,j)=-(-eta(i,j)*dt/((p*dq)**2)-eta(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-eta(i,j)*dt/(4.0d0*(p*dsin(q))**2)-vq(2*i-1,2*j-1)*(dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j)*dsin(q))-dsin(q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-2)*dsin(q)))/2.0d0)
         c(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)+eta(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-vq(2*i-1,2*j-1)*dsin(q+dq/2.0d0)*(1.0d0-vq(2*i-1,2*j)*dsin(q+dq))/2.0d0)
         d(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)+vp(2*i-1,2*j-1)*(p-dp/2.0d0)*(1.0d0+vp(2*i-2,2*j-1)*(p-dp))/2.0d0-0.0d0*eta(i,j)*dt/(p*dp))
         e(i,j)=-(-eta(i,j)*dt/(dp**2)-1.0d0*eta(i,j)*dt/(p*dp)-eta(i,j)*dt/(4.0d0*(p*dsin(q))**2)-vp(2*i-1,2*j-1)*(dp+p*((p+dp/2.0d0)*vp(2*i,2*j-1)+(p-dp/2.0d0)*vp(2*i-2,2*j-1)))/2.0d0)
         f(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)+eta(i,j)*dt/(p*dp)-vp(2*i-1,2*j-1)*(p+dp/2.0d0)*(1.0d0-vp(2*i,2*j-1)*(p+dp))/2.0d0)

      end do
      end do

      do i=2,n
      do j=2,n

         p=pb+float(i-1)/float(n)*(pm-pb)
         q=qm-float(j-1)/float(n)*qm
         ab(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)+vq(2*i-1,2*j-2)*dsin(q-dq/2.0d0)*(1.0d0+vq(2*i-1,2*j-3)*dsin(q-dq))/2.0d0)
         bb(i,j)=-(-etab(i,j)*dt/((p*dq)**2)-etab(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-etab(i,j)*dt/(4.0d0*(p*dsin(q))**2)-(vq(2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j-1)*dsin(q))-vq(2*i-1,2*j-2)*dsin(q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-1)*dsin(q)))/2.0d0)
         cb(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)+etab(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-vq(2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0-vq(2*i-1,2*j+1)*dsin(q+dq))/2.0d0)
         db(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)-0.0d0*etab(i,j)*dt/(p*dp)+vpb(2*i-1,2*j-1)*(p-dp/2.0d0)*(1.0d0 +vpb(2*i-2,2*j-1)*(p-dp))/(2.0d0))
         eb(i,j)=-(-etab(i,j)*dt/(dp**2)-1.0d0*etab(i,j)*dt/(p*dp)-etab(i,j)*dt/(4.0d0*(p*dsin(q))**2)-dvp(2*i-1,2*j-1)*dt/2.0d0-vpb(2*i-1,2*j-1)*(dp+vpb(2*i,2*j-1)*p*(p+dp/2.0d0)+(p-dp/2.0d0)*p*vpb(2*i-2,2*j-1))/2.0d0)
         fb(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)+etab(i,j)*dt/(p*dp)-vpb(2*i-1,2*j-1)*(p+dp/2.0d0)*(1.0d0 -vpb(2*i,2*j-1)*(p+dp))/2.0d0)

      end do
      end do

   end subroutine define_matrix_elements

end module staticfields
