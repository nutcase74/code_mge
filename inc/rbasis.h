interface

  subroutine rbfmq_basis(x, xs, cc, ns, xdim, rbasis)
  implicit none
  integer, parameter       :: dp=kind(1.0d0)
  integer, intent(in)      :: ns, xdim
  real(dp), dimension(xdim), intent(in)       :: x
  real(dp), dimension(ns, xdim), intent(in)   :: xs
  real(dp), dimension(ns), intent(inout)      :: rbasis
  real(dp), intent(in) :: cc
  end subroutine rbfmq_basis

  subroutine rbfinvmq_basis(x, xs, cc, ns, xdim, rbasis)
  implicit none
  integer, parameter       :: dp=kind(1.0d0)
  integer, intent(in)      :: ns, xdim
  real(dp), dimension(xdim), intent(in)       :: x
  real(dp), dimension(ns, xdim), intent(in)   :: xs
  real(dp), dimension(ns), intent(inout)      :: rbasis
  real(dp), intent(in) :: cc
  end subroutine rbfinvmq_basis

  subroutine rbfgauss_basis(x, xs, cc, ns, xdim, rbasis)
  implicit none
  integer, parameter       :: dp=kind(1.0d0)
  integer, intent(in)      :: ns, xdim
  real(dp), dimension(xdim), intent(in)       :: x
  real(dp), dimension(ns, xdim), intent(in)   :: xs
  real(dp), dimension(ns), intent(inout)      :: rbasis
  real(dp), intent(in) :: cc
  end subroutine rbfgauss_basis

  subroutine rbfplate_basis(x, xs, cc, ns, xdim, rbasis)
  implicit none
  integer, parameter       :: dp=kind(1.0d0)
  integer, intent(in)      :: ns, xdim
  real(dp), dimension(xdim), intent(in)       :: x
  real(dp), dimension(ns, xdim), intent(in)   :: xs
  real(dp), dimension(ns), intent(inout)      :: rbasis
  real(dp), intent(in) :: cc
  end subroutine rbfplate_basis

  subroutine mqfconst(x, n, dof, idx, ns, refcen, lref, csq)
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), parameter :: epsil = 1.0E-2
  integer, intent(in) :: n, dof
  real(dp), dimension(n,dof), intent(in) :: x
  integer, dimension(:), pointer :: idx
  integer, intent(inout)  :: ns 
  real(dp), dimension(dof), intent(inout) :: refcen
  real(dp), intent(inout) :: lref, csq
  end subroutine mqfconst

end interface
