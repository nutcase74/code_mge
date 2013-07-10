interface

subroutine trilinear_approx(vs, coord, xn)
  implicit none
  integer, parameter  :: dp=kind(1.0d0)
  real(dp), dimension(8,3), intent(in) :: vs
  real(dp), dimension(3),   intent(in) :: coord
  real(dp), dimension(3), intent(inout):: xn
end subroutine trilinear_approx

subroutine trilinear_factor(a,b,c,coef8)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), intent(in) :: a,b,c
  real(dp), dimension(8), intent(inout) :: coef8
end subroutine trilinear_factor

subroutine func_trilinear(vs, coef)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), dimension(8,3), intent(in) :: vs
  real(dp), dimension(7,3), intent(inout) :: coef
end subroutine func_trilinear

subroutine func_trilinear_jocabian(vs, xc, coef, b)
!
! compute Jocabian for trilinear interpolation
!
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), dimension(8,3), intent(in) :: vs
  real(dp), dimension(3),   intent(in) :: xc
  real(dp), dimension(7,3), intent(inout) :: coef
  real(dp), dimension(3,3), intent(inout) :: b
end subroutine func_trilinear_jocabian

end interface
