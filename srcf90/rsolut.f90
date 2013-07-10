subroutine read_solut(blk, imodel, error)
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), pointer :: blk
  character(len=*), intent(in) :: imodel
  integer, intent(inout) :: error

  type(mgrid), pointer :: b
  integer :: io, code, nconv, nconv_param, low(4), upp(4)
  integer :: i, j, k, d, ib, jb, kb
  character(len=64) :: fname

  error = 0
  nconv_param = 5

  fname = trim(imodel)//"_restart.dat"
  call channel(io)
  open (unit=io, file=fname, status="unknown", action="read", form="formatted")

  print *,"...READING RESTART FILE "//trim(fname)
  b=>blk
  do while (associated(b))

    read(unit=io,fmt="(5I5)") code, low, upp !ib,jb,kb,nconv
!    write(*,fmt="(a,9I5)") "CODE LOW UPP ", code,low,upp

    if (code/=b%code .or. upp(1)/=b%ib .or. upp(2)/=b%jb.or. upp(3)/=b%kb &
                     .or. upp(4)/=nconv_param) then
      print *," INCONSISTENT FILE DOF & GRID DIMENSION "
      stop
    end if

    read(unit=io,fmt="(5f16.7)") ((((b%w(i,j,k,d),i=low(1),upp(1)),j=low(2),upp(2)),&
                                                  k=low(3),upp(3)),d=low(4),upp(4))

    b=>b%nxt
  end do
  close(io)
  return
end subroutine read_solut
