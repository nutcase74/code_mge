subroutine movfin(blk,mgl,mlvl,nmesh,mblks)
! prolong from coarse grid to finer grid
! ===============================================================
  use typelib4
  use typepara
  implicit none

  integer, parameter :: dp=kind(1.0d0)

  type(mgrid), pointer :: blk
  integer, intent(in)  :: mlvl,mgl
  integer, intent(in)  :: nmesh
  integer, dimension(nmesh), intent(in) :: mblks
! local variables
  integer :: n
  type(mgrid), pointer :: btr

!  if(mgl-1.gt.mlvl) return

  do n=1,mblks(mgl)
  !  if(mgl-1.le. mlvl) then
       call addw(blk, mgl, mlvl, n)
  !  end if
  end do

  do n=1,mblks(mgl)
    btr=>get_mgrid(blk,mgl,n)
    call depvars(btr%w, btr%p)
  end do

  call bc_embed(blk,mgl)
!!**********************************************************
!! removed OVERLAPP IN EULER AS THEY ARE OBSOLETE
  call bc_overlap(blk,mgl)
!************************************************************
  do n=1,mblks(mgl)
    btr=>get_mgrid(blk,mgl,n)
    call bcond( btr )
  end do

  return
end subroutine movfin
