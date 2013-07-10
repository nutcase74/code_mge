subroutine init_flow_mp(blk, lrflow, imodel)
!
! flow field initilization with uniform flow
!
!
  use typelib4
  use mpilib

  implicit none

  integer, parameter :: dp=kind(1.0d0)
  type(mgrid),      pointer    :: blk
  character(len=*), intent(in) :: imodel
  logical,          intent(in) :: lrflow

  type(mgrid), pointer :: btr
  integer :: i, j, k, n, error
  integer :: ib,jb,kb

!  integer :: mgl
!  real(dp) :: qq

  print *," STARTING INIT FLOW_MP, lrflow ", lrflow

  if (lrflow) then
    
    if (myrank==imaster) call read_solut(blk, imodel, error)
    call bcast_grid(blk, cvar_dat, error)
    call bcast_movbc(blk)

  else

    btr=>blk
    do while (associated(btr))

      ib = btr%ib
      jb = btr%jb
      kb = btr%kb

      print *," INIT_FLOW MP::ASSOCiATED b%vt ", associated(btr%vt), " VISCL ", viscl 

      if (eqntype == 'n' .or. eqntype == 'N') then
        btr%vt = viscl
      end if

!      print *," DONE READING IN FILE, IM ", myrank, " b%VT =",btr%vt 

      do k=0,kb
        do j=0,jb
          do i=0,ib
            btr%w(i,j,k,1) = rho0
            btr%w(i,j,k,2) = rho0*u0
            btr%w(i,j,k,3) = rho0*v0
            btr%w(i,j,k,4) = rho0*w0
            btr%w(i,j,k,5) = rho0*h0 - p0
          end do
        end do
      end do
      btr=>btr%nxt
    end do
  end if



  btr=>blk
  do while (associated(btr))

    call depvars( btr%w, btr%p)

    if (execmode==1) then
      !! 
      !! for unsteady calculation, give the initial value of flow variables at 
      !! time step "n" and "n-1" in the implicit backward time difference
      !!
      btr%w1 = btr%w
      btr%w2 = btr%w
    end if

    !if (myrank/=imaster) then
    call bc_cutcell(btr)
    call bcond(btr)
    !endif
    btr=>btr%nxt
  end do
  return
end subroutine init_flow_mp
