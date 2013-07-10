subroutine ceflux(w,p,ss,fs)
!
! CALCULATE INVISCID FLUX CONTRIBUTION ON FACE SS OF A CELL
!
! w  : conserv. variables
! p  : pressure
! ss : face-vector
! 
! 2012-08-02 :: CHAMGE ARRAY TO DEFER TYPE
!
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), dimension(:), intent(in) :: w
  real(dp), dimension(:), intent(in) :: ss
  real(dp),               intent(in) :: p
  real(dp), dimension(:), intent(out) :: fs

  fs(1) = sum(w(2:4)*ss(1:3))
  fs(2) = fs(1)*w(2)/w(1) + p*ss(1)
  fs(3) = fs(1)*w(3)/w(1) + p*ss(2)
  fs(4) = fs(1)*w(4)/w(1) + p*ss(3)
  fs(5) = fs(1)*(w(5) + p)/w(1)

  return
end subroutine ceflux
 
