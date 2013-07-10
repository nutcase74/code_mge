subroutine init_flow(blk, imodel)
!
! flow field initilization with uniform flow
! =================================================================
!
!
  use typelib4
  use typepara
  use ctrlib
  implicit none

  integer, parameter :: dp=kind(1.0d0)

  type(mgrid), pointer :: blk
  character(len=*), intent(in) :: imodel

  type(mgrid), pointer :: btr

  integer :: i, j, k, n, error
  integer :: ib,jb,kb
  integer :: mgl
  real(dp) :: qq, ptmp, Ttmp, rat

  if (iflow==1) then
    call read_solut(blk, imodel, error)
  else

    btr=>blk
    do while (associated(btr))

      ib = btr%ib
      jb = btr%jb
      kb = btr%kb

      if(eqntype == 'n' .or. eqntype == 'N') then
        btr%vt = viscl
      end if

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

    call depvars(0,0,0,btr%ib, btr%jb, btr%kb, 5, btr%w, btr%p)

    if(runtype==1) then
      !! 
      !! for unsteady calculation, give the initial value of flow variables at 
      !! time step "n" and "n-1" in the implicit backward time difference
      !!
      btr%w1 = btr%w
      btr%w2 = btr%w
    end if

    call bc_cutcell(btr)
    call bcond(btr)

    btr=>btr%nxt
  end do
  return
end subroutine init_flow
