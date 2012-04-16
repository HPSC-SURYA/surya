! read input, allocate variables, and define some functions that are commonly use throughout the code

module initialization

      implicit none
      save  ! might need this to "preserve data values"

      ! general namelist declarations
      integer :: nmax, lmax, n
      character(len=32) :: dr_filename, mc_filename, init_filename, radial_filename, bottom_filename, run_filename, final_filename
		
      ! part 1 namelist declarations
      integer :: irelax
      double precision :: tmax, v0, et0, et1, al0, dt
		
      ! part2 namelist declarations		
      double precision :: pm, lbf,  pw, qmf
      double precision :: pb, qm, dp, dq

      ! some constants
      double precision :: pi=3.1415926535897932384626433832795028841972d+0

      ! main variables, size to be allocated once nmax is known
      double precision, allocatable :: a(:,:),b(:,:),c(:,:),d(:,:),e(:,:),f(:,:)
      double precision, allocatable :: vp(:,:),vq(:,:),a1(:),b1(:),c1(:),r(:),phi(:,:),phib(:,:),ss1(:),uub(:),uu(:),al(:,:),dom(:,:),ub(:,:),ra(:)
      double precision, allocatable :: ab(:,:),bb(:,:),cb(:,:),db(:,:),eb(:,:),fb(:,:)
      double precision, allocatable :: eta(:,:),dror(:,:),drot(:,:),vp1(:,:),vq1(:,:),psi(:,:),etab(:,:),ss1_p(:),etab2(:,:),dvp(:,:),vpb(:,:),deta(:,:)
		
      ! these particular variables used to be in common blocks:
      double precision, allocatable :: u(:,:), pl(:,:), sn(:)
      double precision :: t
      integer :: it, l

      ! loop counters and scalar coordinates
      integer :: i, j, k, kend
      double precision :: p, q
		

contains

   subroutine read_input()
   ! read in the input parameters from inputs.txt

      namelist /general/ nmax, lmax, dr_filename, mc_filename, init_filename, radial_filename, bottom_filename, run_filename, final_filename
      namelist /part1/ irelax, tmax, v0, et0, et1, al0, dt
      namelist /part2/ pm, lbf, pw, qmf

      open(99, file='input.txt', status='old')

      read(99, nml=general)
         n=nmax-1

      read(99, nml=part1)

      read(99, nml=part2)
         pb=lbf*pm                  ! bottom of the domain, lbf*R
         qm=qmf*4.0d0*atan(1.0d0)   ! range of theta, qmf*PI

   end subroutine read_input


   subroutine allocate_variables()
   ! allocate variables based on problem size

      allocate( a(nmax,nmax),b(nmax,nmax),c(nmax,nmax),d(nmax,nmax),e(nmax,nmax),f(nmax,nmax) )
      allocate( vp(2*nmax,2*nmax),vq(2*nmax,2*nmax),a1(nmax),b1(nmax),c1(nmax),r(nmax),phi(nmax,nmax),phib(nmax,nmax),ss1(nmax),uub(nmax),uu(nmax),al(nmax,nmax),dom(nmax,nmax),u(nmax,nmax),ub(nmax,nmax),ra(nmax) )
      allocate( ab(nmax,nmax),bb(nmax,nmax),cb(nmax,nmax),db(nmax,nmax),eb(nmax,nmax),fb(nmax,nmax) )
      allocate( eta(nmax,nmax),dror(nmax,nmax),drot(nmax,nmax),vp1(2*nmax,2*nmax),vq1(2*nmax,2*nmax),psi(2*nmax,2*nmax),etab(nmax,nmax),ss1_p(nmax),etab2(2*nmax,2*nmax),dvp(2*nmax,2*nmax),vpb(2*nmax,2*nmax),deta(2*nmax,2*nmax) )
      allocate( pl(nmax,lmax),sn(nmax) )

   end subroutine allocate_variables

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! everything from here on was previously in standalone.f90
   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   subroutine tridag(a,b,c,d,x,n)
   ! SUBROUTINE FOR INVERTING TRIDIAGONAL MATRIX
   ! this is the original algorithm, translated into F90 compliant format

      integer, intent(in) :: n
      double precision, dimension(nmax), intent(in) :: a, b, c, d
      double precision, dimension(nmax), intent(inout) :: x

      integer :: j
      double precision :: gam(nmax), bet

      if (.not.(b(1).eq.0.d0)) then
         bet=b(1)
         x(1)=d(1)/bet
         do j=2,n
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j)*gam(j)
            if(.not.(bet.eq.0.d0)) x(j)=(d(j)-a(j)*x(j-1))/bet
         end do
         do j=n-1,1,-1
          x(j)=x(j)-gam(j+1)*x(j+1)
         end do
      end if

   end subroutine tridag
	
	subroutine tridag_wiki(a,b,c,v,x,n)
      ! From Wikipedia (Tridiagonal Matrix Algorithm)

      implicit none
		!        a - sub-diagonal (means it is the diagonal below the main diagonal)
		!        b - the main diagonal
		!        c - sup-diagonal (means it is the diagonal above the main diagonal)
		!        v - right part
		!        x - the answer
		!        n - number of equations
 
      integer,intent(in) :: n
      real(8),dimension(n),intent(in) :: a,b,c,v
      real(8),dimension(n),intent(out) :: x
      real(8),dimension(n) :: bp,vp
      real(8) :: m
      integer i
 
		! Make copies of the b and v variables so that they are unaltered by this sub
        bp(1) = b(1)
        vp(1) = v(1)
 
        !The first pass (setting coefficients):
		do i = 2,n
         m = a(i)/bp(i-1)
         bp(i) = b(i) - m*c(i-1)
         vp(i) = v(i) - m*vp(i-1)
		end do
 
         x(n) = vp(n)/bp(n)
        !The second pass (back-substition)
		do i = n-1, 1, -1
          x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
		end do
 
	end subroutine tridag_wiki


	subroutine coeff(n,ss1,dq,j,cf)
   !     THREE SUBROUTINES NEEDED FOR UPPER BOUNDARY CONDITION :
   !       The subroutine plgndr calculates Legendre polynomials.

      ! double precision :: cf

      integer, intent(in) :: n, j
      double precision, intent(in) :: ss1(nmax)
      double precision, intent(out) :: cf

      ! this input does nothing
      double precision, intent(in) :: dq

      cf=0.0d0
      do l=1,96
         cf=cf+ss1(l)*pl(j,l)
      end do

	end subroutine coeff



	function ss(n,dq)

      integer, intent(in) :: n
      double precision, intent(in) :: dq

      double precision :: y(nmax)
      double precision :: term1, term2
      integer :: i, k1, k2
      double precision :: ss

      do j=1,n+1
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

	end function ss


   function plgndr(l,m,x)
      
      double precision :: plgndr
      integer, intent(in) :: l, m
      double precision, intent(in) :: x

      integer :: ll
      double precision ::  pmm, pmmp1, somx2, fact, pll

      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.0d0) stop 'function plgndr: bad arguments'

      pmm=1.0d0

      if(m.gt.0) then
         somx2=dsqrt((1.0d0-x)*(1.0d0+x))
         fact=1.0d0
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.0d0
         end do
      end if

      if(l.eq.m) then
         plgndr=pmm
      else
         pmmp1=x*(2*m+1)*pmm
         if(l.eq.m+1) then
            plgndr=pmmp1
         else
            do ll=m+2,l
               pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm=pmmp1
               pmmp1=pll
            end do
            plgndr=pll
         end if
      end if

   end function plgndr


   ! Error Function
   function erf(x)

      double precision :: erf
      double precision, intent(in) :: x

      double precision :: p, a1, a2, a3, a4, a5
      double precision :: sig, xabs, tmp

      data p/0.3275911d0/, a1/0.254829592d0/, a2/-0.284496736d0/, a3/1.421413741d0/, a4/-1.453152027d0/, a5/1.061405429d0/

      sig=sign(1.0d0,x)
      xabs=dabs(x)
      if (xabs.gt.20) then
         erf=sig
      else
         tmp=1.0d0/(1.0d0+p*xabs)
         erf=1.0d0-((((a5*tmp+a4)*tmp+a3)*tmp+a2)*tmp+a1)*tmp*exp(-xabs*xabs)
         erf=erf*sig
      endif

   end function erf


end module initialization














