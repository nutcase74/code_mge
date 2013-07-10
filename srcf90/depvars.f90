subroutine depvars_all(cv,dv)
!
!  EVAL DEPENDENT VARAIBLES, (PRESENTLY,  ONLY PRESSURE IS CALCULATED)
!  
  use typepara
  implicit none
  integer, parameter  :: dp=kind(1.0d0)

!  real(dp), dimension(:,:,:,:), intent(in)  :: cv
!  real(dp), dimension(:,:,:), intent(inout) :: dv

  real(dp), dimension(:,:,:,:), pointer  :: cv
  real(dp), dimension(:,:,:),   pointer  :: dv

  integer :: i,j,k
  real(dp) :: qq
  integer, dimension(4) :: low,upp

  low = (/ (lbound(cv,i),i=1,4) /)
  upp = (/ (ubound(cv,i),i=1,4) /)
  !print *, "DEPVAR LOW/UPP: ",low,upp

  do k=low(3), upp(3)
    do j=low(2), upp(2)
      do i=low(1), upp(1)

        qq        = 0.5d0*sum(cv(i,j,k,2:4)*cv(i,j,k,2:4) )/cv(i,j,k,1)
        dv(i,j,k) = (gamma -1.0d0)*dim(cv(i,j,k,5),qq)

      end do
    end do
  end do
  return
end subroutine depvars_all


subroutine depvars_one(cv,dv)
!
!  EVAL DEPENDENT VARAIBLES, (PRESENTLY,  ONLY PRESSURE IS CALCULATED)
!  
  use typepara
  implicit none
  integer, parameter  :: dp=kind(1.0d0)
  real(dp), dimension(:), intent(in) :: cv
  real(dp), intent(inout) :: dv
  real(dp) :: qq

   qq = 0.5d0*sum(cv(2:4)*cv(2:4) )/cv(1)
   dv = (gamma -1.0d0)*dim(cv(5),qq)

  return
end subroutine depvars_one
