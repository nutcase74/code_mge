subroutine init_flow(blk, lrflow, imodel)
!
! flow field initilization with uniform flow
! =================================================================
!
!
  use typelib4
  implicit none

  integer, parameter :: dp=kind(1.0d0)

  type(mgrid), pointer :: blk
  character(len=*), intent(in) :: imodel
  logical, intent(in) :: lrflow

  type(mgrid), pointer :: btr

  integer :: i, j, k, error
  integer :: ib,jb,kb
!  integer :: mgl
!  real(dp) :: qq

  print *," INIT FLOW lrflow ", lrflow

  if ( lrflow ) then
    !
    !  READ FLOW FROM RESTART FILE
    !
    print *," READ SOLUT "
    call read_solut(blk, imodel, error)
    print *," DONE READ SOLUT "

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
          !qq = 0.5d0*(b%w(i,j,k,2)**2.0d0+b%w(i,j,k,3)**2.0d0+b%w(i,j,k,4)**2.0d0)/b%w(i,j,k,1)
          !b%p(i,j,k) = (gamma -1.0d0)*(b%w(i,j,k,5) -qq)

            if(ibm_on) then
              if(btr%blank(i,j,k) == 0) then
                btr%w(i,j,k,2:4) = 0.0d0
                btr%w(i,j,k,5)   = rho0*gamma*ei0-p0
              end if
            end if

          end do
        end do
      end do
      btr=>btr%nxt

    end do
  end if

  btr=>blk
  do while (associated(btr))

    !call depvars(0,0,0,btr%ib, btr%jb, btr%kb, 5, btr%w, btr%p)

    call depvars( btr%w, btr%p)

    if (execmode==1) then
      !! 
      !! for unsteady calculation, give the initial value of flow variables at 
      !! time step "n" and "n-1" in the implicit backward time difference
      !!
      btr%w1 = btr%w
      btr%w2 = btr%w
    end if

    if(ibm_on) then
!       call bc_IBM1(btr)
!       call bc_pressure(btr)
!       call bc_velocity(btr) 
    else
       call bc_cutcell(btr)
    end if

    call bcond(btr)

    btr=>btr%nxt
  end do
  return
end subroutine init_flow
