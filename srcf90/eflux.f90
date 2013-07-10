subroutine eflux(b)
!
!
! calculate inviscid flux thru grid boundary, it's easy for orthogonal
! uniform Cartesian grids, cell volume calculated in a simple fashion.
!
! ======================================================================

  use typelib4
  use typepara

  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b

  integer :: i, j, k, n, m, zz
  integer :: il,jl,kl,ie,je,ke,ib,jb,kb
  integer, dimension(3) :: icut

  real(dp), dimension(5) :: fs1, fs2
  real(dp), dimension(:,:,:,:), pointer :: efsx,efsy,efsz
! **********************************************************************

  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb
  allocate( efsx(0:ib,0:jb,0:kb,5),&
            efsy(0:ib,0:jb,0:kb,5),&
            efsz(0:ib,0:jb,0:kb,5) )

  efsx = 0.0d0
  efsy = 0.0d0
  efsz = 0.0d0

! euler flux in i-direction
  do k=2,kl
    do j=2,jl
      do i=2,ie
        call ceflux(b%w(i,  j,k,:),b%p(i,  j,k),b%si(i,  j,k,:),fs1)
        call ceflux(b%w(i-1,j,k,:),b%p(i-1,j,k),b%si(i-1,j,k,:),fs2)
        efsx(i,j,k,:) = 0.5d0*(fs1(:)+fs2(:)) 
      end do
    end do
  end do

  do k=2,kl
    do j=2,jl
      do i=2,il
        b%dw(i,j,k,:) = efsx(i,j,k,:) -efsx(i+1,j,k,:)
      end do
    end do
  end do

! euler flux in j-direction
  if(i2d==0) then
    do k=2,kl
      do j=2,je
        do i=2,il
          call ceflux(b%w(i,j,  k,:),b%p(i,j,  k),b%sj(i,j,  k,:),fs1)
          call ceflux(b%w(i,j-1,k,:),b%p(i,j-1,k),b%sj(i,j-1,k,:),fs2)
          efsy(i,j,k,:) = 0.5d0*(fs1(:)+fs2(:)) 
        end do
      end do
    end do

    do k=2,kl
      do j=2,jl
        do i=2,il
          b%dw(i,j,k,:) = b%dw(i,j,k,:) + efsy(i,j,k,:) -efsy(i,j+1,k,:)
        end do
      end do
    end do
  end if

! euler flux in k-direction
  do k=2,ke
    do j=2,jl
      do i=2,il
        call ceflux(b%w(i,j,k,  :),b%p(i,j,k  ),b%sk(i,j,k,  :),fs1)
        call ceflux(b%w(i,j,k-1,:),b%p(i,j,k-1),b%sk(i,j,k-1,:),fs2)
        efsz(i,j,k,:) = 0.5d0*(fs1(:)+fs2(:)) 
      end do
    end do
  end do

  do k=2,kl
    do j=2,jl
      do i=2,il
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + efsz(i,j,k,:) -efsz(i,j,k+1,:)
      end do
    end do
  end do

   !debug
   !do k = 2,kl
   !  do j = 2,jl
   !    do i = 2,il
   !      if (b%dw(i,j,k,1) .ne. b%dw(i,j,k,1) ) then
   !        print *,"eflux1 ijk",i,j,k,"dw",b%dw(i,j,k,1),b%w(i,j,k,1)
   !        stop
   !      end if
   !    end do
   !  end do
   !end do

! --- INVISCID FLUX THRU WALL, i.e. CUT-CELLS,  
  if(.not. ibm_on)  call efluxb(b,efsx,efsy,efsz)

  deallocate(efsx,efsy,efsz)

   !debug
   !do k = 2,kl
   !  do j = 2,jl
   !    do i = 2,il
   !      if (b%dw(i,j,k,1) .ne. b%dw(i,j,k,1) ) then
   !        print *,"eflux2 ijk",i,j,k,"dw",b%dw(i,j,k,1),b%w(i,j,k,1)
   !        stop
   !      end if
   !    end do
   !  end do
   !end do

  return

contains

subroutine ceflux(w,p,ss,fs)
!
! CALCULATE INVISCID FLUX CONTRIBUTION ON FACE SS OF A CELL
!
! w  : conserv. variables
! p  : pressure
! ss : face-vector
! 
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), dimension(5), intent(in) :: w
  real(dp), dimension(3), intent(in) :: ss
  real(dp), intent(in) :: p
  real(dp), dimension(5), intent(out) :: fs

  fs(1) = sum(w(2:4)*ss(1:3))
  fs(2) = fs(1)*w(2)/w(1) + p*ss(1)
  fs(3) = fs(1)*w(3)/w(1) + p*ss(2)
  fs(4) = fs(1)*w(4)/w(1) + p*ss(3)
  fs(5) = fs(1)*(w(5) + p)/w(1)

  return
end subroutine ceflux
 


end subroutine eflux
