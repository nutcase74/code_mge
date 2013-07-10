subroutine movco(blk,mgl,nmesh,mblks)
! restrict flow variables and residual from finer grid to coarse
! grid by subroutine collc
! =========================================================================
  use typelib4
  use typepara
  implicit none
  integer, parameter :: dp=kind(1.0d0)

! input/output variables
  type(mgrid), pointer :: blk
  integer, intent(in) :: mgl
  integer, intent(in) :: nmesh
  integer, dimension(nmesh), intent(in) :: mblks

! local variables
  integer :: n,i,j,k
  type(mgrid), pointer :: btr


  do n=1,mblks(mgl)

    call collc(blk,mgl,n)

    btr=>get_mgrid(blk, mgl, n)

    call depvars(btr%w, btr%p)

    do i=2,btr%il
      do j=2,btr%jl
        do k=2,btr%kl
          btr%wn(i,j,k,:) = btr%w(i,j,k,:)
        end do
      end do
    end do

  end do

!!**********************************************************
!! removed EMBED & OVERLAPP IN EULER AS THEY ARE OBSOLETE
!!    call bc_embed(blk,mgl)   ! 
    call bc_overlap(blk,mgl)
!************************************************************
  do n=1,mblks(mgl)
    btr=>get_mgrid(blk, mgl, n)
    call bcond( btr )
  end do

  return
end subroutine movco


