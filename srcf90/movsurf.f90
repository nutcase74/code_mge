subroutine movsurf(mgrida) !, level, popul)
!
! GET MOVING BOUDNARY CONDITIONS AT RAY-SURFACE INTERSECTION POINTS, VALUES HAVE BEEN
! EVALUATED AND STORED IN {MGRIDA}
!
  use typelib4
  implicit none
  type(mgrid), pointer :: mgrida, mgtr
  integer :: m, b, k, level
  integer, dimension(:), pointer :: popul

  call eval_popul(mgrida,popul)
  level = size(popul)

  do m = 1, level
    do b = 1, popul(m)

      mgtr    =>get_mgrid(mgrida, m, b)

      !one_blk =>get_mgrid(blk,    m, b)
      !print *," ASSOCIATED mgtr   ", associated(mgtr)
      !print *," ASSOCIATED one_blk", associated(one_blk)
      !write(unit=io,fmt=*) "LEVEL ",m,"BLOCK ",b
      !print *," ibid ", ibid, " mgtr%id ", mgtr%id, " nmir ", iblk(b,m)%znmir

      print *,mgtr%code, " MOVSURF ", mgtr%nmir

      do k = 1, mgtr%nmir
        mgtr%z(k)%pnew = mgtr%z(k)%psurf + mgtr%z(k)%disp

        !one_blk%z(k)%norv  = mgtr%z(k)%norv 
        !one_blk%z(k)%vel   = mgtr%z(k)%vel
        !write(unit=io,fmt=*) " ppt ",iblk(b,m)%ixnode(k)%new_Ppt, " rn0 ",iblk(b,m)%ixnode(k)%rn,&
        !" rn ",iblk(b,m)%ixnode(k)%new_rn, " qb ",iblk(b,m)%ixnode(k)%qb
      end do
    end do
  end do

  return
end subroutine movsurf
