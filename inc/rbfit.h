  subroutine rbf_coef(x,f,n,ndim, ns,xs, coef, cc, kk, err, funcs)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in)                      :: n, ndim, ns
  real(dp), dimension(n,ndim), intent(in)  :: x
  real(dp), dimension(n), intent(in)       :: f
  real(dp), dimension(ns,ndim), intent(in) :: xs
  real(dp), dimension(ns), intent(inout)   :: coef
  real(dp), intent(in)                     :: cc
  real(dp), intent(inout)                  :: kk
  real(dp), intent(inout)                  :: err

  interface
    subroutine funcs(x,xs,cc,ns,xdim,rbasis)
    implicit none
    integer, parameter       :: dp=kind(1.0d0)
    integer, intent(in)      :: ns, xdim
    real(dp), dimension(xdim), intent(in)       :: x
    real(dp), dimension(ns, xdim), intent(in)   :: xs
    real(dp), dimension(ns), intent(inout)      :: rbasis
    real(dp), intent(in) :: cc
    end subroutine funcs
  end interface
  end subroutine rbf_coef

  subroutine rbf_nodal(x,n,ndim, ns,xs, xfield, coef, cc, funcs)
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, intent(in)                      :: n, ndim, ns
  real(dp), dimension(n,ndim), intent(in)  :: x
  real(dp), dimension(ndim), intent(in)    :: xfield
  real(dp), dimension(ns,ndim), intent(in) :: xs
  real(dp), dimension(n), intent(inout)    :: coef
  real(dp), intent(in)                     :: cc

  interface
    subroutine funcs(x,xs,cc,ns,xdim,rbasis)
    implicit none
    integer, parameter       :: dp=kind(1.0d0)
    integer, intent(in)      :: ns, xdim
    real(dp), dimension(xdim), intent(in)       :: x
    real(dp), dimension(ns, xdim), intent(in)   :: xs
    real(dp), dimension(ns), intent(inout)      :: rbasis
    real(dp), intent(in) :: cc
    end subroutine funcs
  end interface
  end subroutine rbf_nodal
