! this module contains the standalone subroutines and functions used throughout SURYA

module standalone

   use initialization

contains	

   subroutine tridag(a,b,c,d,x,n)
   ! SUBROUTINE FOR INVERTING TRIDIAGONAL MATRIX
   ! this is the original algorithm, translated into F90 compliant format

      integer, intent(in) :: n
      double precision, dimension(nmax), intent(in) :: a, b, c, d
      double precision, dimension(nmax), intent(inout) :: x

      double precision :: gam(nmax) ! local

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
      
      integer, intent(in) :: l, m
      double precision, intent(in) :: x

      double precision ::  pmm=1.0d0, somx2, fact

      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.0d0) stop 'function plgndr: bad arguments'

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

      double precision :: sig, xabs

      data p/0.3275911d0/,a1/0.254829592d0/,a2/-0.284496736d0/,a3/1.421413741d0/,a4/-1.453152027d0/,a5/1.061405429d0/

      sig=sign(1.0d0,x)
      xabs=dabs(x)
      if (xabs.gt.20) then
         erf=sig
      else
         t=1.0d0/(1.0d0+p*xabs)
         erf=1.0d0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-xabs*xabs)
         erf=erf*sig
      endif

   end function erf

end module standalone
