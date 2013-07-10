subroutine write_solut_one(blk, job, error)
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)

  type(mgrid), pointer         :: blk
  character(len=*), intent(in) :: job
  integer, intent(inout)       :: error

  integer :: m, n, io !, nconv
  integer :: i, j, k, d, low(4), upp(4)
  character(len=64) :: fname, tcode
  type(mgrid), pointer :: b
!
! N.B. NCONV has to be defined in mgrid, define here for now
!
!  nconv = 5
  error = 0
!  fname = trim(job)//"_restart.dat"
!  call channel(io)
!  open (unit=io, file=fname, status="unknown", action="write", form="formatted")

  b => blk
  do while (associated(b))

    call i2char(b%code,tcode,3,0)
    fname = trim(job)//"_"//trim(tcode)//".dat"
    call channel(io)
    open (unit=io, file=fname, status="unknown", action="write", form="formatted")

    low = (/ (lbound(b%w,d),d=1,4) /)
    upp = (/ (ubound(b%w,d),d=1,4) /)

    write(unit=io,fmt="(5I5)") b%code, low, upp
    write(unit=io,fmt="(5f16.7)") ((((b%w(i,j,k,d),i=low(1),upp(1)),j=low(2),upp(2)),&
                                                   k=low(3),upp(3)),d=low(4),upp(4))

    close(io)

    b=>b%nxt
  end do

  return
end subroutine write_solut_one


subroutine write_solut(blk, job, error)
!
! WRITE FLOW BLOCK AS RESTART 
!
! NOTE: ENTIRE RANGE OF BLOCK IS WRITTEN 
!
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)

  type(mgrid), pointer         :: blk
  character(len=*), intent(in) :: job
  integer, intent(inout)       :: error

  integer :: m, n, io, low(4), upp(4) !nconv
  integer :: i, j, k, d
  character(len=64) :: fname
  type(mgrid), pointer :: b
!
! N.B. NCONV has to be defined in mgrid, define here for now
!
!  nconv = 5
  error = 0
  fname = trim(job)//"_restart.dat"
  call channel(io)
  open (unit=io, file=fname, status="unknown", action="write", form="formatted")
  print *,"....WRITING RESTART FILE "//trim(fname)

  b => blk
  do while (associated(b))

    low = (/ (lbound(b%w,d),d=1,4) /)
    upp = (/ (ubound(b%w,d),d=1,4) /)

    write(unit=io,fmt="(5I5)") b%code,low,upp
    write(unit=io,fmt="(5f16.7)") ((((b%w(i,j,k,d),i=low(1),upp(1)),j=low(2),upp(2)),&
                                                   k=low(3),upp(3)),d=low(4),upp(4))

!    write(unit=io,fmt="(5I5)") b%code,b%ib,b%jb,b%kb,nconv
!    write(unit=io,fmt="(5f16.7)") ((((b%w(i,j,k,d),i=1,b%ib),j=1,b%jb),&
!                                                   k=1,b%kb),d=1,nconv)

    b=>b%nxt
  end do
  close(io)
  return
end subroutine write_solut
