! this module contains the standalone subroutines and functions used throughout SURYA

module standalone

contains	
	
	subroutine tridag(a,b,c,v,x,n)
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
 
	end subroutine tridag


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
	end subroutine coeff



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
	end function ss


	function plgndr(l,m,x)
        implicit real*8(a-h,o-z)
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
        return
	end function plgndr


   ! Error Function
	function erf(x)
      implicit real*8 (a-h,o-z)
      data p/0.3275911d0/,a1/0.254829592d0/,a2/-0.284496736d0/,a3/1.421413741d0/,a4/-1.453152027d0/,a5/1.061405429d0/
      sig=sign(1.0d0,x)
      xabs=dabs(x)
      if (xabs.gt.20) then
!     
!     ***  avoidance  of underflow ***
!
      erf=sig
      else
      t=1.0d0/(1.0d0+p*xabs)
      erf=1.0d0-((((a5*t+a4)*t+a3)*t+a2)*t+a1)*t*exp(-xabs*xabs)
      erf=erf*sig
      endif
      return
	end function erf

end module standalone
