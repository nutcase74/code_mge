subroutine progon3(a,b,c,d,x,n)
!
!c***********************************************************
!c        Program Progon3 for solving linear system         *
!c        of the form                                       *
!c                                                          *
!c        a(1)*x(1) + b(1)*x(2)      + c(1)*x(n) = d(1)     *
!c        ............................................      *
!c        c(i)*x(i-1) + a(i)*x(i)  + b(i)*x(i+1) = d(i)     *
!c        .............................................     *
!c        b(n)*x(1)     + a(n)*x(n-1) + b(n)*x(n)= d(n)     *
!c        with tridiagonal matrix by sweep method           *
!c    n - number of equations                               *
!c  a,b,                                                    *
!c  c,d - arrays (of size n) containing elements of the     *
!c        matrix' diagonals and of the right parts          *
!
!
!c  w,s,t,                                                  *
!c  u,v - work arrays (of size n+1)                         *
!c        Output:                                           *
!c    x - array (of size n) containing solution of the      *
!c        system                                            *
!c***********************************************************
!
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  integer, intent(in) :: n
  real(dp), dimension(n), intent(in)    :: a,b,c,d
  real(dp), dimension(n), intent(inout) :: x
!
! LOCAL VARIABLES
!
  real(dp), dimension(n+1) :: w,s,t,u,v
  integer :: i, i1
  real(dp) :: z 
!  print *,"=============================="
!  print *," c(1) ", c(1)
!  print *," b(n) ", b(n)
!  print *,"=============================="
  u(1) = 0.0d0
  v(1) = 0.0d0
  w(1) = 1.0d0
  do i = 1,n
    i1 = i + 1
    z  = 1.0d0/(a(i) + c(i)*v(i))
    v(i1) = -b(i)*z
    u(i1) = (-c(i)*u(i) + d(i))*z
    w(i1) = -c(i)*w(i)*z
  enddo
  s(n) = 1.0d0
  t(n) = 0.0d0
  do i = n-1,1, -1
    s(i) = v(i+1)*s(i+1) + w(i+1)
    t(i) = v(i+1)*t(i+1) + u(i+1)
  enddo
  x(n) = (d(n) - b(n)*t(1) - c(n)*t(n-1))/(a(n) + b(n)*s(1) + c(n)*s(n-1))
  do i = 1, n-1
    x(i) = s(i)*x(n) + t(i)
  enddo
  return
end subroutine progon3 
