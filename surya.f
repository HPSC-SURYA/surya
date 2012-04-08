! Formatted and annotated version


!                  The solar dynamo code SURYA
!
!
!   INTRODUCTION :
!     This is a 2-dimensional code for solving the kinematic dynamo
!     equation in the solar convection zone with a meridional
!     circulation.  The potential user of this code should first
!     read the accompanying 'Guide', which explains the code and
!     describes how it can be run. This code is divided in several
!     parts for the benifit of the users. These different parts
!     are explained in the 'Guide'.  Somebody familiar with the
!     solar dynamo proble should be able to run this code after
!     reading the 'Guide'.  This 'Guide' divides the potential
!     users in three levels and provides appropriate instructions
!     for users at these three different levels.
!              
!     This code was developed over the years at Indian Institute
!     of Science, Bangalore, by Arnab Rai Choudhuri and his
!     successive PhD students - Mausumi Dikpati, Dibyendu Nandy,
!     Piyali Chatterjee.
!
        implicit real*8(a-h,o-z)
		! nmax sets mesh points in both directions
		! there are 3 other places where this number is set
        parameter (nmax=257,lmax=96)
        common u(nmax,nmax),t,it
        common/lmx/l
        common /plsincos/pl(nmax,lmax),sn(nmax)
        double precision a(nmax,nmax),b(nmax,nmax),c(nmax,nmax),
     & d(nmax,nmax),
     &  e(nmax,nmax),f(nmax,nmax),vp(2*nmax,2*nmax),vq(2*nmax,2*nmax)
     &  ,a1(nmax),b1(nmax),c1(nmax),r(nmax),phi(nmax,nmax),
     &  phib(nmax,nmax),ss1(nmax),uub(nmax),
     &  uu(nmax),al(nmax,nmax),dom(nmax,nmax),ub(nmax,nmax),ra(nmax)
        double precision ab(nmax,nmax),bb(nmax,nmax),cb(nmax,nmax),
     &  db(nmax,nmax),eb(nmax,nmax),fb(nmax,nmax),eta(nmax,nmax),
     &  dror(nmax,nmax),drot(nmax,nmax),vp1(2*nmax,2*nmax),
     &  vq1(2*nmax,2*nmax),psi(2*nmax,2*nmax),etab(nmax,nmax),
     &  ss1_p(nmax),etab2(2*nmax,2*nmax),dvp(2*nmax,2*nmax),
     &  vpb(2*nmax,2*nmax),deta(2*nmax,2*nmax)
        external ss,erf
        
        n=nmax-1
	  pi=4.0d0*atan(1.0d0)
!
!----------------------------------------------------------------------
!	
!	PART I. This is the ONLY part of the code in which Level I users
!	may want to make some changes.  Read Sect. 3 of the 'Guide' to
!	learn about the variables specified in this part.  
!

! irelax=0 means code-init, irelax=1 means use init.dat
	irelax=0  
! max time of computation (1e8 s)
	tmax=1.0d0 
! amplitude of meridional flow (m/s)
	v0=-29.0d0
! poloidal diffusivity (1e12 cm2/s)
	et0=2.6d0
! toroidal diffusivity (1e12 cm2/s)
	et1=0.04d0
! alpha-effect amplitude (m/s)
	al0=25.0d0
! time step (1e8 s)
! NOTE: Must satisfy the Courant condition (explicit?)
	dt=0.0002d0

!
!----------------------------------------------------------------------
!
!	PART II.  This is the part where the alpha coefficient, 
!	diffusivity, differential rotation and meridional circulation 
!	are specified. Level II Users wishing to make changes in 
!	this part should read Sect. 4 of the 'Guide'.
!

! outer radius (1e10 cm)
	pm=6.96d0
! inner radius (1e10 cm)
	pb=0.55d0*pm
! ???
	pw=500.0d0
! range of theta
	qm=pi
! delta_radius
	dp=(pm-pb)/float(n)
! delta_theta
	dq=-qm/float(n)

!	
!	Here begin the do loops to calculate profiles of alpha, 
!	diffusivity and differential rotation at all the grid 
!	points. The differential rotation is written in the file 
!	'diffrot.dat'.
!
	open(25,file='diffrot.dat',status='unknown',access='append')   
	do i = 1, nmax
	do j = 1, nmax
	    ra(i)=pb+float(i-1)/float(n)*(pm-pb) 
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
	    co = dcos(q)
! ALPHA PROFILE
	   al(i,j) = 0.25d0*al0*co*(1.00d0+
     &	     erf((p-0.95d0*pm)/(0.025d0*pm)))*(1.00d0-erf((p
     &       -pm)/(0.025d0*pm)))
! DIFFUSIVITY PROFILES
          eta(i,j) = 0.000220d0 + (et0/2.00d0)*(1.0d0 +
     &       erf((p-0.7d0*pm)/(0.025d0*pm)))
          etab(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 +
     &       erf((p-0.72d0*pm)/(0.025d0*pm)))+ (et0/2.0d0)
     &       *(1.0d0+erf((p-0.975d0*pm)/(0.025d0*pm))) 
! SOLAR DIFFERENTIAL ROTATION PROFILE
          dom(i,j) = 271.9d0 + 0.50d0*(1.000d0 +
     &       erf((p-0.7d0*pm)/(0.025d0*pm)))*(289.5d0 -
     &       39.4d0*co*co - 42.2d0*co*co*co*co -
     &       271.9d0)         
          write(25,27)q,p,dom(i,j)
  27      format(3(f13.5,1x))
		end do
    	end do
    	close(25)
! 
!	The do loops are closed.  We also need the diffusivity at
!	mid-points in the grid for some calculations.  This is 
!	obtained and stored now.
!
	do i=1,2*n+1
	do j=1,2*n+1
	    p=pb+float(i-1)*(pm-pb)/float(2*n)
	    etab2(i,j) = 0.00022d0 + (et1/2.00d0)*(1.0d0 +
     &       erf((p-0.72d0*pm)/(0.025d0*pm)))+ (et0/2.0d0)
     &       *(1.0d0+erf((p-0.975d0*pm)/(0.025d0*pm))) 
	end do
	end do
!	
!	Derivatives of differential rotation are obtained and stored
!	now for future use.
!
	    do i=2,n
	do j=2,n
	    dror(i,j)=(dom(i+1,j)-dom(i-1,j))/(2.0d0*dp)  
	    drot(i,j)=(dom(i,j+1)-dom(i,j-1))/(2.0d0*dq)
		end do
   		end do
!
! MERIDIONAL CIRCULATION. Note that this is calculated both at the 
!   grid-points and mid-points. First the stream function psi(i,j)is 
!   calculated. It is written in the file 'psi.dat'.	
!
       beta1=1.5d0
       beta2=1.8d0
       pp=0.61d0*pm
       del=2.0000001d0
       gm=3.47d0
       p0=(pm-pb)/4.00d0        
	 open(82,file='psi.dat',status='unknown',access='append')     
	 do i=1,2*n+1
	 do j=1,2*n+1
           p=pb+float(i-1)/float(2*n)*(pm-pb)
           q=qm-float(j-1)/float(2*n)*qm
!
           if(q.le.pi/2.0d0) then
           exq0=dexp(-beta1*q**del)
           exqm=dexp(beta2*(q-pi/2.0d0))
           gau=dexp(-((p-p0)/gm)**2)
           psi(i,j)=(p-pp)*dsin(pi*(p-pp)/(pm-pp))*(1.00d0-exq0)*
     &         (1.00d0-exqm)*gau
           end if 
!
           if(q.gt.pi/2.0d0) then
           exq0=dexp(-beta1*(pi-q)**del)
           exqm=dexp(beta2*(pi/2.0d0-q))
           gau=dexp(-((p-p0)/gm)**2)
           psi(i,j)=-(p-pp)*dsin(pi*(p-pp)/(pm-pp))*(1.00d0-exq0)*
     &         (1.00d0-exqm)*gau
           end if
! 
           if(p.le.pp)then
           psi(i,j)=0.0d0
           end if
           write(82,34)q,p,psi(i,j)
 34        format(3(f13.5,1x))
         end do
         end do
           close(82)
!
!	Now components of velocity are calculated at grid-points and 
!	mid-points from the stream function psi(i,j).
!
	do i=2,2*n
	do j=2,2*n
           p=pb+float(i-1)*(pm-pb)/float(2*n)
           q=qm-float(j-1)*qm/float(2*n)
           vp1(i,j)=0.95d0*v0*(psi(i,j+1)-psi(i,j-1))/(p**2*dsin(q)*dq
     &          *(pm/p-0.95d0)**1.5)
           vq1(i,j)=-0.95d0*v0*(psi(i+1,j)-psi(i-1,j))/(p*dsin(q)*dp
     &          *(pm/p-0.95d0)**1.5)
           deta(i,j)=(etab2(i+1,j)-etab2(i-1,j))/(dp)
           vpb(i,j)=(vp1(i,j)-deta(i,j))*dt/(2.d0*p*dp)
           vp(i,j)=vp1(i,j)*dt/(2.d0*p*dp)
           vq(i,j)=vq1(i,j)*dt/(2.d0*p*dsin(q)*dq)
        end do
        end do
!
!------------------------------------------------------------
!
!	PART II-A.  Legendre polynomials are calculated here with the
!	help of a subroutine and stored for future use.
!
           do l=1,96
            cl=1.0d0/(1.0d0+float(l)/float(l+1)*(pm/pw)**(2*l+1))
            ss1_p(l)=float(2*l+1)/(float(l)*pm)*(cl-
     &        float(l)/float(l+1)*(1.0d0-cl))
        do j=1,n+1
        qm=4.0d0*atan(1.0d0)
        q=qm+float(j-1)*dq
        x=dcos(q)
        pl(j,l) = plgndr(l,1,x)
        sn(j) = dsin(q)
        end do
           end do
!
!------------------------------------------------------------------
!
!	PART III. This is the part where initial values of the 
!	variables u(i,j) and ub(i,j) are chosen, before beginning 
!	the time advancement.

 42    format(4(f13.7,1x))
	if(irelax.eq.0) then
	  do i = 1,n+1
	  do j = 1,n+1
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
	  u(i,j) = 0.0d0
	  ub(i,j) = 1.0*dsin(1.0d0*q)*dsin(pi
     &      *((p-pb)/(pm-pb)))
	  end do
	  end do
	  go to 10
	end if

        open(12,file='init.dat',status='unknown')
        do i=1,n+1
         do j=1,n+1
          read(12,42) q,p,u(i,j),ub(i,j)
           if(j.eq.n+1) then
            u(i,j)=0.0d0
           else
            u(i,j)=u(i,j)/(p*dsin(q))
           end if
         end do
        end do
        close(12)
 44     format(4(f13.5,1x))
!
!----------------------------------------------------------------
!
!	PART IV. For advancing in time, we need some matrix 
!	elements arising out of the difference schemes used. We 
!	calculate and store the matrix elements now.  Look at 
!	Sect. 5 of the 'Guide' for a discussion about these 
!	matrix elements.  Only Level III Users need to bother
!	about the matrix elements.
!

! pre-compute coefficients for solvers
! includes mesh geometry and velocity, since these never change
10	do i=2,2*n
	do j=2,2*n	          
	  	dvp(i,j)=(vp1(i+1,j)-vp1(i-1,j))/(dp)
	end do
	end do  
       do i=2,n
         do j=2,n
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
        a(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)+
     & vq(2*i-1,2*j-1)*dsin(q-dq/
     & 2.0d0)*(1.0d0+vq(2*i-1,2*j-2)*dsin(q-dq))/2.0d0)
        b(i,j)=-(-eta(i,j)*dt/((p*dq)**2)-
     & eta(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-
     & eta(i,j)*dt/(4.0d0*(p*dsin(q))**2)-vq(
     & 2*i-1,2*j-1)*(dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j)*dsin(q))-
     & dsin(q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-2)*dsin(q)))/2.0d0)
        c(i,j)=-(eta(i,j)*dt/(2.0d0*(p*dq)**2)+
     & eta(i,j)*dt/(2.0d0*dtan(q)*p*p*
     & dq)-vq(2*i-1,2*j-1)*dsin(q+dq/2.0d0)*(1.0d0-vq(2*i-1,2*j)*dsin(q+
     & dq))/2.0d0)
        d(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)+vp(2*i-1,2*j-1)*(p-dp/2.0d0)*
     & (1.0d0+vp(2*i-2,2*j-1)*(p-dp))/2.0d0-0.0d0*eta(i,j)*dt/(p*dp))
        e(i,j)=-(-eta(i,j)*dt/(dp**2)-
     & 1.0d0*eta(i,j)*dt/(p*dp)-eta(i,j)*dt/(4.0d0*
     & (p*dsin(q))**2)-vp(2*i-1,2*j-1)*(dp+p*((p+dp/2.0d0)*vp(2*i,2*j-
     & 1)+(p-dp/2.0d0)*vp(2*i-2,2*j-1)))/2.0d0)
        f(i,j)=-(eta(i,j)*dt/(2.0d0*dp**2)+
     & eta(i,j)*dt/(p*dp)-vp(2*i-1,2*j-1)*
     & (p+dp/2.0d0)*(1.0d0-vp(2*i,2*j-1)*(p+dp))/2.0d0)
        end do
        end do
!
         do i=2,n
         do j=2,n
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
        ab(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)+
     & vq(2*i-1,2*j-2)*dsin(q-dq/
     & 2.0d0)*(1.0d0+vq(2*i-1,2*j-3)*dsin(q-dq))/2.0d0)
        bb(i,j)=-(-etab(i,j)*dt/((p*dq)**2)-
     & etab(i,j)*dt/(2.0d0*dtan(q)*p*p*dq)-
     & etab(i,j)*dt/(4.0d0*(p*dsin(q))**2)-(vq(
     & 2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0+vq(2*i-1,2*j-1)*dsin(q))
     &-vq(2*i-1,2*j-2)*dsin(
     & q-dq/2.0d0)*(1.0d0-vq(2*i-1,2*j-1)*dsin(q)))/2.0d0)
        cb(i,j)=-(etab(i,j)*dt/(2.0d0*(p*dq)**2)+
     & etab(i,j)*dt/(2.0d0*dtan(q)*p*p*
     & dq)-vq(2*i-1,2*j)*dsin(q+dq/2.0d0)*(1.0d0-vq(2*i-1,2*j+1)*dsin(q+
     & dq))/2.0d0)
        db(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)-0.0d0*etab(i,j)*dt/(p*dp)+
     & vpb(2*i-1,2*j-1)*(p-dp/2.0d0)*(
     & 1.0d0 +vpb(2*i-2,2*j-1)*(p-dp))/(2.0d0))
        eb(i,j)=-(-etab(i,j)*dt/(dp**2)-
     & 1.0d0*etab(i,j)*dt/(p*dp)-etab(i,j)*dt/(4.0d0*
     & (p*dsin(q))**2)-dvp(2*i-1,2*j-1)*dt/2.0d0-vpb(2*i-1,2*j-1)*(dp+
     & vpb(2*i,2*j-1)*p*(p+dp/2.0d0)+(p-dp/2.0d0)*p*
     & vpb(2*i-2,2*j-1))/2.0d0)
        fb(i,j)=-(etab(i,j)*dt/(2.0d0*dp**2)+
     & etab(i,j)*dt/(p*dp)-vpb(2*i-1,2*j-1)*
     & (p+dp/2.0d0)*(1.0d0 -vpb(2*i,2*j-1)*(p+dp))/2.0d0)
          end do
          end do
!
!-------------------------------------------------------------------------
!
!	PART V. Now we come to the central part of the programme where the
!	time advancement takes place. 'k' is the counter to keep tab on
!	the time.  In each step of the do loop 'do k=1,kend', the magnetic
!	fields are advanced through one time step 'dt'.
!
        kend = tmax/dt
        t=0.0d0
        do k=1,kend
!
! PART V-A. CALCULATING TIME ADVANCED MAGNETIC FIELDS AT INTERIOR GRID 
! POINTS. Sect. 5 of the 'Guide' explains how magnetic fields are 
! advanced by solving tridiagonal matrices with the subroutine 'tridag'
!
! solve explicit side of A evolution
          do i=2,n
           do j=2,n
            phi(i,j)=-d(i,j)*u(i-1,j)+(-e(i,j)+1.0d0)*u(i,j)
     &        -f(i,j)*u(i+1,j)+al(i,j)*ub(i,j)*dt/2.0d0
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
! solve explicit sde of B evolution
          do i=2,n
           do j=2,n
            p=pb+float(i-1)/float(n)*(pm-pb)
            q=qm-float(j-1)/float(n)*qm
            br = (u(i,j+1)*dsin(q+dq)-u(i,j-1)*dsin(q-dq))*dt/(4.0d0*dq)
	    bt=(u(i-1,j)*(p-dp)-u(i+1,j)*(p+dp))*dsin(q)*dt/(4.0d0*p*dp)
            phib(i,j)=-db(i,j)*ub(i-1,j)+(-eb(i,j)+1.0d0)*ub(i,j)
     &        -fb(i,j)*ub(i+1,j)+br*dror(i,j)+bt*drot(i,j)
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
            phib(i,j)=phib(i,j)
     &        +br*dror(i,j)+bt*drot(i,j)
           end do
          end do
!
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
!
! PART V-B. BOUNDARY CONDITIONS. These are discussed in Sect. 2.1 of .
! the 'Guide'. 
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
!
! PART V-C.  INCORPORATING THE MAGNETIC BUOYANCY by radial transport
! of toroidal field exceeding a specified critical value:
	   tau = .0088d0
	   nmb = tau/dt
           kmb = k/nmb
           kch = nmb*kmb
           if(k.eq.kch)then
	     n2 = n/2.
	     nt = n - 6.		
           do i=1,n2
           do j=2,n
            if (ra(i).gt.0.71*pm)then
            if (ABS(ub(i,j)).GT.0.8d0)then
             jer=j
             ier=i
             qjer=qm-float(jer-1)*qm/float(n)	
             ber=ub(ier,jer)
             bold=ub(nt,jer)
             ub(nt,jer)=((ra(ier)/ra(nt))*0.5d0*ber) + bold
             ub(ier,jer)=0.5d0*ber
              open(25,file='ber.dat',status='unknown',access='append')
              write(25,47) t,qjer,ber
  47          format(3(d15.8,1x))
              close(25)
            end if
            end if
           end do
           end do
           end if
!
          t=t+dt
!
! PART V-D. WRITING IN THE FILES 'rad.dat' and 'butbot.dat' after 
! certain intervals of time
!
        if(t.gt.0.0d0) then
	  p20 = pi/20.
	  istep = t/p20
	  tdiff = t - float(istep)*p20
 	    if(tdiff.ge.dt)go to 55
  	  open(17,
     &file='rad.dat',status='unknown',access='append') 
  	  open(18,
     &file='butbot.dat',status='unknown',access='append') 
 	 do j = 2, n
            q=qm-float(j-1)/float(n)*qm
 	    br = -(u(n-1,j-1)*dsin(q-dq) - u(n-1,j+1)
     & *dsin(q+dq))/(2.0d0*dq*dsin(q))/pm
 	    write(17,37) t, q, br
 	    write(18,37) t, q, ub(45,j)
 	  end do
  37      format(3(f13.7,1x))
 	  close(17)
 	  close(18)
 	  close(19)
  55    end if	
!
! WRITING IN THE FILE 'run.dat' after every 40 time steps
!
	if(k/40*40.eq.k)then
 	  open(95,
     &file='run.dat',status='unknown',access='append')
! TODO: Generalize for grid
	  write(95,90)t,ub(45,74),ub(44,74),ub(45,54),
     &   u(120,122),u(120,6)
  90      format(6(d15.9,1x))
	  close(95)
	end if

         end do
! The mammoth do loop which began at the beginning of PART V ends
! here. 
!
!---------------------------------------------------------------
!	PART VI. Now the final state is written in the file 
!	'final.dat'
!
         open(10,
     &file='final.dat',status='unknown')
         do i=1,n+1
          do j=1,n+1
           p=pb+float(i-1)/float(n)*(pm-pb)
           q=qm-float(j-1)/float(n)*qm
           write(10,42) q,p,p*dsin(q)*u(i,j),ub(i,j)
          end do
         end do
         close(10)
	
        stop
        end

!
!     MAIN PROGRAM ENDS HERE.
! -----------------------------------------------------
!
!     SUBROUTINE FOR INVERTING TRIDIAGONAL MATRIX :
        subroutine tridag(a,b,c,r,u,n)
        implicit real*8 (a-h,o-z)
        parameter (nmax=257)
        dimension gam(nmax),a(nmax),b(nmax),c(nmax),r(nmax),u(nmax)
        if(b(1).eq.0.d0) pause
        bet=b(1)
        u(1)=r(1)/bet
        do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.d0) pause
          u(j)=(r(j)-a(j)*u(j-1))/bet
  11    continue
        do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
  12    continue
        return
        end


!
!     THREE SUBROUTINES NEEDED FOR UPPER BOUNDARY CONDITION :
!       The subroutine plgndr calculates Legendre polynomials.
        subroutine coeff(n,ss1,dq,j,cf)
        implicit real*8(a-h,o-z)
        external plgndr
        parameter (nmax=257,lmax=96)
        dimension ss1(nmax)
        common u(nmax,nmax),t,it
        common/lmx/l
        common /plsincos/pl(nmax,lmax),sn(nmax)
!       qm=4.0d0*atan(1.0d0)
!       q=qm+float(j-1)*dq
!       x=dcos(q)
        cf=0.0d0
        do l=1,96
!       cf=cf+ss1(l)*plgndr(l,1,x)
        cf=cf+ss1(l)*pl(j,l)
        end do
        return
        end



        function ss(n,dq)
        implicit real*8(a-h,o-z)
        external plgndr
        parameter (nmax=257,lmax=96)
        dimension y(nmax)
        common u(nmax,nmax),t,it
        common/lmx/l
        common /plsincos/pl(nmax,lmax),sn(nmax)
!       pb=3.52d0
        do j=1,n+1
!       qm=4.d0*atan(1.0d0)
!       q=qm+float(j-1)*dq
!       x=dcos(q)
!       y(j)=plgndr(l,1,x)*u(n+1,j)*dsin(q)
        y(j)=pl(j,l)*u(n+1,j)*sn(j)
        end do
        term1=0.0d0
        do k1=2,n,2
        term1=term1+y(k1)
        end do
        term2=0.0d0
        do k2=3,n-1,2
        term2=term2+y(k2)
        end do
        ss=dq/3.0d0*(y(1)+4.0d0*term1+2.0d0*term2+y(n+1))
        return
        end


        function plgndr(l,m,x)
        implicit real*8(a-h,o-z)
        if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.0d0)pause 'bad arguments'
        pmm=1.0d0
        if(m.gt.0) then
          somx2=dsqrt((1.0d0-x)*(1.0d0+x))
          fact=1.0d0
          do 11 i=1,m
             pmm=-pmm*fact*somx2
             fact=fact+2.0d0
  11      continue
        end if
        if(l.eq.m) then
          plgndr=pmm
        else
          pmmp1=x*(2*m+1)*pmm
          if(l.eq.m+1) then
             plgndr=pmmp1
          else
             do 12 ll=m+2,l
                 pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                 pmm=pmmp1
                 pmmp1=pll
  12         continue
           plgndr=pll
          end if
        end if
        return
        end

!
      function erf (x)
      implicit real*8 (a-h,o-z)
      data p/0.3275911d0/,a1/0.254829592d0/,a2/-0.284496736d0/,
     &     a3/1.421413741d0/,a4/-1.453152027d0/,a5/1.061405429d0/
      sig=sign(1.0d0,x)
      xabs=dabs(x)
      if (xabs.gt.20) then
!     
!     *** vermeidung von underflow ***
!     ***  avoidance  of underflow ***
!
       erf=sig
      else
       t=1.0d0/(1.0d0+p*xabs)
       erf=1.0d0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-xabs*xabs)
       erf=erf*sig
      endif
      return
      end

