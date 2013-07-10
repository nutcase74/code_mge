subroutine depvars_all(i0,j0,k0,im,jm,km,nconv,cv,dv)
!
!  EVAL DEPENDENT VARAIBLES, (PRESENTLY,  ONLY PRESSURE IS CALCULATED)
!  
!  use typelib4
  use typepara
  implicit none
  integer, parameter  :: dp=kind(1.0d0)
  integer, intent(in) :: i0,j0,k0,im,jm,km,nconv
! input/output variables
  real(dp), dimension(i0:im,j0:jm,k0:km,nconv), intent(in)    :: cv
  real(dp), dimension(i0:im,j0:jm,k0:km),       intent(inout) :: dv

  integer :: i,j,k
  real(dp) :: qq

  do k=k0,km
    do j=j0,jm
      do i=i0,im
        qq        = 0.5d0*sum(cv(i,j,k,2:4)**2.0d0)/cv(i,j,k,1)
        dv(i,j,k) = (gamma -1.0d0)*dim(cv(i,j,k,5),qq)
      end do
    end do
  end do

  return
end subroutine depvars_all


subroutine depvars_one(nconv,cv,dv)
!
!  EVAL DEPENDENT VARAIBLES, (PRESENTLY,  ONLY PRESSURE IS CALCULATED)
!  
  use typepara
  implicit none
  integer, parameter  :: dp=kind(1.0d0)
  integer, intent(in) :: nconv
! input/output variables
  real(dp), dimension(nconv), intent(in) :: cv
  real(dp), intent(inout) :: dv

  real(dp) :: qq

  qq = 0.5d0*sum(cv(2:4)*cv(2:4) )/cv(1)
  dv = (gamma -1.0d0)*dim(cv(5),qq)

  return
end subroutine depvars_one
