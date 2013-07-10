subroutine blk2mgrid(blk, mgr, n, popul, p_tot, atcellcenter)
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), pointer :: blk
  type(mgrid), pointer :: mgr
  integer, intent(in) ::  n
  integer, dimension(n), intent(in) :: popul
  real(dp), intent(in) :: p_tot
  integer, intent(inout) :: atcellcenter

  integer :: m, ilevel
  type(mgrid), pointer :: ptr, one_blk

  do ilevel = 1, n
    do m = 1, popul(ilevel)
      ptr    =>get_mgrid(mgr, ilevel, m)
      one_blk=>get_mgrid(blk, ilevel, m)
      if ( (.not. associated(ptr)).or. (.not. associated(one_blk)) ) then
        print *," FATAL ERROR: FAIL TO RETRIEVE MGRID FOR LEVEL/SERIAL ", n,m
        stop
      end if
      call upload_solut(one_blk, ptr, p_tot, atcellcenter)
    end do
  end do
  return
end subroutine blk2mgrid


