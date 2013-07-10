subroutine efluxb(b,efsx,efsy,efsz)
!***********************************************
!
! EULER FLUX AT WALL BOUNDARY
!
!***********************************************
  use typepara
  use typelib4

  implicit none
  integer, parameter :: dp=kind(1.0d0)
! input/output variables
  type(mgrid), intent(inout) :: b
  real(dp), dimension(0:b%ib,0:b%jb,0:b%kb,5), intent(inout) :: efsx,efsy,efsz
!
! local variables
!
  type(txcut), pointer :: tmp
  integer :: i, j, k, id, ncl, zz, m
  integer :: is,js,ks
  integer, dimension(3) :: icut, isol
  real(dp) :: psol
  real(dp) :: dis2, dis4, ep
  real(dp) :: lm, lp, lam
  real(dp), dimension(5) :: wsol
  real(dp), dimension(3) :: ss, sc
  real(dp), dimension(5) :: fs1, fs2, fs3, fs, fa
  real(dp), dimension(3,5) :: diw

  do m = 1,b%nmir

    tmp=>b%z(m)

    id   = tmp%faceid
    zz   = tmp%zz
    isol = tmp%scid
    icut = tmp%ccid
    psol = tmp%psol
    wsol = tmp%wsol

! compute the flux at the boundary

    i  = icut(1) 
    j  = icut(2) 
    k  = icut(3) 
    is = isol(1) 
    js = isol(2) 
    ks = isol(3) 

    if(zz==1) then
      sc(:) = b%si(i,j,k,:)
      ss(:) = b%si(is,js,ks,:)
    else if(zz==2) then
      sc(:) = b%sj(i,j,k,:)
      ss(:) = b%sj(is,js,ks,:)
    else if(zz==3) then
      sc(:) = b%sk(i,j,k,:)
      ss(:) = b%sk(is,js,ks,:)
    end if

    if(abs(b%w(i,j,k,1)) < 1.0e-8) then
       print *, "Zero density in bcwall_eflux : ",i,j,k,b%w(i,j,k,1)
       stop
    end if

    call ceflux(b%w(i,j,k,:),b%p(i,j,k),sc,fs1)
    call ceflux(wsol(:),psol,ss,fs2)
    fs = 0.5d0*(fs1+fs2) 

    if(zz == 1) then

      b%dw(i,j,k,:) = b%dw(i,j,k,:) - (efsx(i,j,k,:) -efsx(i+1,j,k,:))
      ! x-direction
      if(id == 4) then
        call ceflux(b%w(i-1,j,k,:),b%p(i-1,j,k),b%si(i-1,j,k,:),fs3)
        fa = 0.5d0*(fs1+fs3)
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + fa(:) - fs(:)
      else if(id == 1) then
        call ceflux(b%w(i+1,j,k,:),b%p(i+1,j,k),b%si(i+1,j,k,:),fs3)
        fa = 0.5d0*(fs1+fs3)
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + fs(:) - fa(:)
      end if

    else if(zz == 2) then

      b%dw(i,j,k,:) = b%dw(i,j,k,:) - (efsy(i,j,k,:) -efsy(i,j+1,k,:))
      ! y-direction
      if(id == 5) then
        call ceflux(b%w(i,j-1,k,:),b%p(i,j-1,k),b%sj(i,j-1,k,:),fs3)
        fa = 0.5d0*(fs1+fs3)
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + fa(:)-fs(:)
      else if(id == 2) then
        call ceflux(b%w(i,j+1,k,:),b%p(i,j+1,k),b%sj(i,j+1,k,:),fs3)
        fa = 0.5d0*(fs1+fs3)
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + fs(:) - fa(:)
      end if

    else if(zz == 3) then

      b%dw(i,j,k,:) = b%dw(i,j,k,:) - (efsz(i,j,k,:) -efsz(i,j,k+1,:))
      ! z-direction
      if(id == 6) then
        call ceflux(b%w(i,j,k-1,:),b%p(i,j,k-1),b%sk(i,j,k-1,:),fs3)
        fa = 0.5d0*(fs1+fs3)
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + fa(:) - fs(:)
      else if(id == 3) then
        call ceflux(b%w(i,j,k+1,:),b%p(i,j,k+1),b%sk(i,j,k+1,:),fs3)
        fa = 0.5d0*(fs1+fs3)
        b%dw(i,j,k,:) = b%dw(i,j,k,:) + fs(:) - fa(:)
      end if

    end if

  end do

  return
end subroutine efluxb

