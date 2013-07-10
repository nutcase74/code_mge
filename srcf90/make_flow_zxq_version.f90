subroutine make_flow(blk, imcyc)
!
! flow field initilization with uniform flow
! =================================================================
!
!
  use typelib4
  use typepara
  implicit none

  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), pointer :: blk
  integer, intent(in)  :: imcyc

  type(mgrid), pointer :: b
  integer :: ib,jb,kb

  b=>blk
  do while (associated(b))

    ib = b%ib
    jb = b%jb
    kb = b%kb

!    print *," ****************************************************** "
!    print *," INIT_XF    CODE ", b%code
!    print *," ib, jb, kb : ",ib,jb,kb

    ! allocate necessary arrays as public variables saved through computation

    allocate(b%w (0:ib,0:jb,0:kb,5))
    allocate(b%p (0:ib,0:jb,0:kb))
    allocate(b%w1(0:ib,0:jb,0:kb,5))
    allocate(b%w2(0:ib,0:jb,0:kb,5))
    allocate(b%ws(0:ib,0:jb,0:kb,5))
    allocate(b%dw(0:ib,0:jb,0:kb,5))
    allocate(b%fw(0:ib,0:jb,0:kb,5))
    allocate(b%wn(ib-1,jb-1,kb-1,5))
    allocate(b%wr(ib-1,jb-1,kb-1,5))
    allocate(b%qf(ib-1,jb-1,kb-1,5))
    allocate(b%radi(0:ib,0:jb,0:kb),b%radj(0:ib,0:jb,0:kb))
    allocate(b%radk(0:ib,0:jb,0:kb),b%vdtim(0:ib,0:jb,0:kb))
    allocate(b%epsx(0:ib,0:jb,0:kb),b%epsy(0:ib,0:jb,0:kb), b%epsz(0:ib,0:jb,0:kb))
    allocate(b%rtrms( imcyc ))
    if(eqntype == 'n' .or. eqntype == 'N') &
         allocate(b%vt(0:ib,0:jb,0:kb))

    b%w    = 0.0d0
    b%p    = 0.0d0
    b%w1   = 0.0d0
    b%w2   = 0.0d0
    b%ws   = 0.0d0
    b%dw   = 0.0d0
    b%fw   = 0.0d0
    b%wn   = 0.0d0
    b%wr   = 0.0d0
    b%qf   = 0.0d0
    b%radi = 0.0d0
    b%radj = 0.0d0
    b%radk = 0.0d0
    b%vdtim= 0.0d0
    b%epsx = 0.0d0
    b%epsy = 0.0d0
    b%epsz = 0.0d0
    if(eqntype == 'n' .or. eqntype == 'N') &
         b%vt   = viscl
!
! ----- ALLOCATE RTRMS FOR CONVERGENCE HISTORY
!

    b%rtrms = 0.0d0

    b=>b%nxt
  end do

  return
end subroutine make_flow
