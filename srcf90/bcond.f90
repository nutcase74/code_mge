subroutine bcond(b)
! subroutine treating boundary conditions for each non-zero bctype
! six faces are treated independently, use 1-D Riemann invariables.
! bctype = -20: far field
!          -50: symmetric plane
!          -60: periodic boundary condition
! ===================================================================
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b
! local variables
  integer  :: i
! *******************************************************************

  do i = 1,6

    if(b%bctype(i) /= 0) then

      if(b%bctype(i) == -20) then
        call bc_far(b,i)
      else if(b%bctype(i) == -50) then
        call bc_symm(b,i)
      else if(b%bctype(i) == -60) then
        call bc_periodic(b,i)
        if(i .le. 3) then
          if(b%bctype(i+3) /= -60) then
            print *, "Periodic boundary must appear in pairs !!!"
            stop
          end if
        end if
      else if(b%bctype(i) == -80) then
        call bc_wall(b,i)
      else
        print *, "bctype :",b%bctype(i),"not implemented !!!" 
      end if
    end if

  end do

! set values on block corner cells
  call bc_corner(b)

  return
end subroutine bcond
