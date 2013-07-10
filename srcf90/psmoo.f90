subroutine psmoo(b)
!
! subroutine for implicit residual smoothing
! smooth coefficients are explicitly given and constant thru computation.
! solve tridiagonal matrix with progon3()
! ===================================================================

  use typepara
  use typelib4

  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b

! local variables

  integer :: i, j, k, n
  integer :: il,jl,kl
  integer :: ib,jb,kb,mx
  real(dp), dimension(:), pointer :: aa,bb,cc

! *******************************************************************

  if(smoo .le. 0.0d0) return  ! no smoothing

  il = b%il
  jl = b%jl
  kl = b%kl

  ib = b%ib
  jb = b%jb
  kb = b%kb
  
  mx = max(max(ib,jb),kb)
  allocate(aa(mx),bb(mx),cc(mx))
  aa = 0.0d0
  bb = 0.0d0
  cc = 0.0d0

! cells with 0 residual like solid cells would not smooth.
! smooth in i-direction
  do k=2,kl
    do j=2,jl
      do i=2,il
        aa(i) = 1.0d0 + 2.0d0*b%epsx(i,j,k)
        bb(i) = -b%epsx(i,j,k)
        cc(i) = -b%epsx(i,j,k)
      end do
      bb(il)=0.0d0
      cc(2) =0.0d0
      do n=1,5
        call progon3(aa(2:il),bb(2:il),cc(2:il),b%dw(2:il,j,k,n),b%dw(2:il,j,k,n),il-1)
      end do
    end do
  end do

! smooth in y-direction
  do k=2,kl
    do i=2,il
      do j=2,jl
        aa(j) = 1.0d0 + 2.0d0*b%epsy(i,j,k)
        bb(j) = -b%epsy(i,j,k)
        cc(j) = -b%epsy(i,j,k)
      end do
      bb(jl)=0.0d0
      cc(2) =0.0d0
      do n=1,5
        call progon3(aa(2:jl),bb(2:jl),cc(2:jl),b%dw(i,2:jl,k,n),b%dw(i,2:jl,k,n),jl-1)
      end do
    end do
  end do

! smooth in z-direction
  do i=2,il
    do j=2,jl
      do k=2,kl
        aa(k) = 1.0d0 + 2.0d0*b%epsz(i,j,k)
        bb(k) = -b%epsz(i,j,k)
        cc(k) = -b%epsz(i,j,k)
      end do
      bb(kl)=0.0d0
      cc(2) =0.0d0
      do n=1,5
        call progon3(aa(2:kl),bb(2:kl),cc(2:kl),b%dw(i,j,2:kl,n),b%dw(i,j,2:kl,n),kl-1)
      end do
    end do
  end do


  if (associated(aa)) deallocate(aa)
  if (associated(bb)) deallocate(bb)
  if (associated(cc)) deallocate(cc)

  return
end subroutine psmoo

subroutine psmoo_old(b)
!
! subroutine for implicit residual smoothing
! smooth coefficients are explicitly given and constant thru computation.
! ===================================================================
!  use typeblock
!  use typepara

  use typepara
  use typelib4

  implicit none
  integer, parameter :: dp=kind(1.0d0)
! input variables
!  type(block_type), pointer :: b
  type(mgrid), intent(inout) :: b
! local variables
  integer :: i, j, k, n
  integer :: il,jl,kl
  integer :: im, jm, km, ilm, jlm, klm
  integer :: ir, jr, kr, irp, jrp, krp
  real(dp) :: aa, bb
  real(dp), dimension(:,:), pointer :: epsx0, epsy0, epsz0
!  real(dp) :: smoo
! *******************************************************************

  if(smoo .le. 0.0d0) return  ! no smoothing

  ilm = b%ncell(1) !nx =il-1
  jlm = b%ncell(2) !y  =jl-1
  klm = b%ncell(3) !z  =kl-1

  il = b%il
  jl = b%jl
  kl = b%kl
  allocate(epsx0(2:il,2),epsy0(2:jl,2),epsz0(2:kl,2))

! coefficients in x-direction
  
  aa = 1.0d0 +2.0d0*smoo
  bb = -smoo

  epsx0(2,1) = aa
  epsx0(2,2) = bb/aa

  do i=3,ilm
    epsx0(i,1) = aa -bb*epsx0(i-1,2)
    epsx0(i,2) = bb/epsx0(i,1)
  end do
  epsx0(il,1) = aa -bb*epsx0(ilm,2)

! coefficients in y-direction
  
  epsy0(2,1) = aa
  epsy0(2,2) = bb/aa
  do j=3,jlm
    epsy0(j,1) = aa -bb*epsy0(j-1,2)
    epsy0(j,2) = bb/epsy0(j,1)
  end do
  epsy0(jl,1) = aa -bb*epsy0(jlm,2)

! coefficients in z-direction

  epsz0(2,1) = aa
  epsz0(2,2) = bb/aa
  do k=3,klm
    epsz0(k,1) = aa -bb*epsz0(k-1,2)
    epsz0(k,2) = bb/epsz0(k,1)
  end do
  epsz0(kl,1) = aa -bb*epsz0(klm,2)

! cells with 0 residual like solid cells would not smooth.
! smooth in i-direction

  do n=1,5
    do k=2,kl
      do j=2,jl
        if(b%blank(2,j,k)/=0)&
        b%dw(2,j,k,n) = b%dw(2,j,k,n)/epsx0(2,1)
        do i=3,il
          im = i -1
          if(b%blank(i,j,k)/=0)&
          b%dw(i,j,k,n) = (b%dw(i,j,k,n) -bb*b%dw(im,j,k,n))/epsx0(i,1)
        end do
        
        do i=2,ilm
          ir  = il +1 -i
          irp = ir +1
          if(b%blank(ir,j,k)/=0)&
          b%dw(ir,j,k,n) = b%dw(ir,j,k,n) -epsx0(ir,2)*b%dw(irp,j,k,n)
        end do
      end do       

! smooth in y-direction
      do i=2,il
        if(b%blank(i,2,k)/=0)&
        b%dw(i,2,k,n) = b%dw(i,2,k,n)/epsy0(2,1)
        do j=3,jl
          jm = j -1
          if(b%blank(i,j,k)/=0)&
          b%dw(i,j,k,n) = (b%dw(i,j,k,n) -bb*b%dw(i,jm,k,n))/epsy0(j,1)
        end do

        do j=2,jlm
          jr  = jl +1 -j
          jrp = jr +1
          if(b%blank(i,jr,k)/=0)&
          b%dw(i,jr,k,n) = b%dw(i,jr,k,n) -epsy0(jr,2)*b%dw(i,jrp,k,n)
        end do
      end do
    end do

! smooth in z-direction

    do j=2,jl
      do i=2,il
        if(b%blank(i,j,2)/=0)&
        b%dw(i,j,2,n) = b%dw(i,j,2,n)/epsz0(2,1)
      end do
!   end do

    do k=3,kl
      km = k -1
     !do j=2,jl
        do i=2,il
          if(b%blank(i,j,k)/=0)&
          b%dw(i,j,k,n) = (b%dw(i,j,k,n) -bb*b%dw(i,j,km,n))/epsz0(k,1)
        end do
     !end do
    end do

    do k=2,klm
      kr  = kl +1 -k
      krp = kr +1
     !do j=2,jl
        do i=2,il
          if(b%blank(i,j,kr)/=0)&
          b%dw(i,j,kr,n) = b%dw(i,j,kr,n) -epsz0(kr,2)*b%dw(i,j,krp,n)
        end do
      end do
    end do
  end do

  deallocate(epsx0,epsy0,epsz0)
  return
end subroutine psmoo_old
