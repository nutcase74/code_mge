subroutine movsurf(blk, mgrida) !, level, popul)
!
! GET MOVING BOUDNARY CONDITIONS AT RAY-SURFACE INTERSECTION POINTS, VALUES HAVE BEEN
! EVALUATED AND STORED IN {MGRIDA}
!
  use typelib4
  implicit none
  type(mgrid), pointer :: blk, one_blk
  type(mgrid), pointer :: mgrida, mgtr
  integer :: m, b, k, level  !, io, ioerror
  integer, dimension(:), pointer :: popul

  character(len=64) :: fname
  print *,"  MVSURF_COUPLED, ASSO BLK    " , associated(blk)
  print *,"  MVSURF_COUPLED, ASSO MGRIDA " , associated(mgrida)

  !call channel(io)
  !fname="00mvsurf.dat"
  !open (unit=io, file=fname, status="unknown", action="write", form="formatted", &
  !        position="rewind", iostat=ioerror)
  call eval_popul(blk,popul)
  level = size(popul)

  do m = 1, level
    do b = 1, popul(m)
      !ibid = sum(iind(1:m-1))+b

      print *," LEVEL, SERIAL ", m, b

      mgtr    =>get_mgrid(mgrida, m, b)
      one_blk =>get_mgrid(blk,    m, b)

      print *," ASSOCIATED mgtr   ", associated(mgtr)
      print *," ASSOCIATED one_blk", associated(one_blk)

      !write(unit=io,fmt=*) "LEVEL ",m,"BLOCK ",b
      !print *," ibid ", ibid, " mgtr%id ", mgtr%id, " nmir ", iblk(b,m)%znmir

      do k = 1, mgtr%nmir
        one_blk%z(k)%pnew  = mgtr%z(k)%psurf + mgtr%z(k)%disp
        one_blk%z(k)%norv  = mgtr%z(k)%norv 
        one_blk%z(k)%vel   = mgtr%z(k)%vel
        !write(unit=io,fmt=*) " ppt ",one_blk%z(k)%pnew, " rn0 ",one_blk%z(k)%nvec,&
        !" rn ",one_blk%z(k)%norv, " qb ",one_blk%z(k)%vel
        !write(unit=io,fmt=*) " disp ",mgtr%z(k)%disp,  " qb ",one_blk%z(k)%vel
      end do
    end do
  end do
  !close(io)
  return
end subroutine movsurf
