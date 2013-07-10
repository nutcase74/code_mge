subroutine dfluxb(b,istage,dfsx,dfsy,dfsz)
!*******************************************************************************
!
!  Artificial Dissipation at wall boundary
!
!*******************************************************************************
  use typepara
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
! input/output variables
  type(mgrid), intent(inout) :: b
  integer, intent(in) :: istage
  real(dp), dimension(0:b%ib,0:b%jb,0:b%kb,5), intent(in) :: dfsx,dfsy,dfsz
!
! local variables
!
  type(txcut), pointer :: tmp
  integer :: i, j, k, id, ncl, zz, m
  integer, dimension(3) :: icut  !, isol
  real(dp) :: psol
  real(dp) :: dis2, dis4, ep
  real(dp) :: lm, lp, lam
  real(dp), dimension(5) :: wsol, var_refl
  real(dp), dimension(3) :: fp
  real(dp), dimension(5) :: fs1, fs2
  real(dp), dimension(3,5) :: diw
!---------------------------------------------------------------------------------
  do m = 1,b%nmir

    tmp=>b%z(m)

    id   = tmp%faceid
    zz   = tmp%zz
    icut = tmp%ccid
    psol = tmp%psol
    wsol = tmp%wsol
    wsol(5) = wsol(5) + psol  ! replace energy by enthalpy

! compute the flux at the boundary
    i  = icut(1)
    j  = icut(2)
    k  = icut(3)

    if(zz==1) then
      if(id == 4) then
        b%fw(i  ,j,k,:) = b%fw(i  ,j,k,:) -(dfsx(i,  j,k,:) -dfsx(i+1,j,k,:))
        b%fw(i-1,j,k,:) = b%fw(i-1,j,k,:) -(dfsx(i-1,j,k,:) -dfsx(i,j,k,:))
        fp(1) = abs(b%p(i-2,j,k) -2.0d0*b%p(i-1,j,k) +b%p(i,j,k))/&
                abs(b%p(i-2,j,k) +2.0d0*b%p(i-1,j,k) +b%p(i,j,k))
        fp(2) = abs(b%p(i-1,j,k) -2.0d0*b%p(i,j,k) +psol)/&
                abs(b%p(i-1,j,k) +2.0d0*b%p(i,j,k) +psol)
        ep    = max(fp(1),fp(2))
        lm    = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam   = lm
        dis2     = vis2*ep
        dis4     = vis4
        dis4     = dim(dis4,dis2)
        dis2     = dis2*cdis(istage)
        dis4     = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k,:) -b%w(i-1,j,k,:)
        diw(2,:) = wsol(:) -b%w(i,j,k,:)
        diw(3,:) = diw(2,:)
        fs2(:)   = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs2(:)   = dis2*diw(2,:) -dis4*fs2(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs2(:) = fs2(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k,:),wsol(:),fs2(:),b%si(i+1,j,k,:),fs2(:))
        end if

!! modified by tsllc 01-12-10
!! set fp(1) to 0 if (i-3) is negative
        if ((i-3)>=0) then
          fp(3) = abs(b%p(i-1,j,k) -2.0d0*b%p(i-2,j,k) +b%p(i-3,j,k))/&
                  abs(b%p(i-1,j,k) +2.0d0*b%p(i-2,j,k) +b%p(i-3,j,k))
        else
          fp(3) = 0.0d0
        end if
        ep       = max(fp(1),fp(2),fp(3))
        lm       = b%radi(i-1,j,k) +b%radj(i-1,j,k) +b%radk(i-1,j,k)
        lp       = b%radi(i,  j,k) +b%radj(i,  j,k) +b%radk(i,  j,k)
        lam      = 0.5d0*(lm +lp)
        dis2     = vis2*ep
        dis4     = vis4
        dis4     = dim(dis4,dis2)
        dis2     = dis2*cdis(istage)
        dis4     = dis4*cdis(istage)
        diw(1,:) = b%w(i-1,j,k,:) -b%w(i-2,j,k,:)
        diw(2,:) = b%w(i,j,k,:) -b%w(i-1,j,k,:)
        diw(3,:) = wsol(:) -b%w(i,j,k,:)
        fs1(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs1(:) = dis2*diw(2,:) -dis4*fs1(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs1(:) = fs1(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i-1,j,k,:),b%w(i,j,k,:),fs1(:),b%si(i,j,k,:),fs1(:))
        end if

        b%fw(i  ,j,k,:) = b%fw(i  ,j,k,:) +(fs1(:) -fs2(:))
        b%fw(i-1,j,k,:) = b%fw(i-1,j,k,:) +(dfsx(i-1,j,k,:) -fs1(:))

      else if(id==1) then

        b%fw(i  ,j,k,:) = b%fw(i  ,j,k,:) -(dfsx(i,j,k,:) -dfsx(i+1,j,k,:))
        b%fw(i+1,j,k,:) = b%fw(i+1,j,k,:) -(dfsx(i+1,j,k,:) -dfsx(i+2,j,k,:))

        fp(1) = abs(b%p(i+1,j,k) -2.0d0*b%p(i,j,k) +psol)/&
                abs(b%p(i+1,j,k) +2.0d0*b%p(i,j,k) +psol)

        fp(2) = abs(b%p(i+2,j,k) -2.0d0*b%p(i+1,j,k) +b%p(i,j,k))/&
                abs(b%p(i+2,j,k) +2.0d0*b%p(i+1,j,k) +b%p(i,j,k))

!! modified by tsllc 01-12-10
!! set fp(4) to 0 if (i+3) > b%ib
        if ((i+3)<=b%ib) then
          fp(3) = abs(b%p(i+3,j,k) -2.0d0*b%p(i+2,j,k) +b%p(i+1,j,k))/&
                  abs(b%p(i+3,j,k) +2.0d0*b%p(i+2,j,k) +b%p(i+1,j,k))
        else
          fp(3) = 0.0d0
        end if

        ep    = max(fp(1),fp(2),fp(3))
        lm    = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lp    = b%radi(i+1,j,k) +b%radj(i+1,j,k) +b%radk(i+1,j,k)
        lam   = 0.5d0*(lm +lp)
        dis2  = vis2*ep
        dis4  = vis4
        dis4  = dim(dis4,dis2)
        dis2  = dis2*cdis(istage)
        dis4  = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k,:) -wsol(:)
        diw(2,:) = b%w(i+1,j,k,:) -b%w(i,j,k,:)
        diw(3,:) = b%w(i+2,j,k,:) -b%w(i+1,j,k,:)
        fs2(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs2(:) = dis2*diw(2,:) -dis4*fs2(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs2(:) = fs2(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k,:),b%w(i+1,j,k,:),fs2(:),b%si(i+1,j,k,:),fs2(:))
        end if

        ep    = max(fp(1),fp(2))
        lp    = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam   = lp
        dis2  = vis2*ep
        dis4  = vis4
        dis4  = dim(dis4,dis2)
        dis2  = dis2*cdis(istage)
        dis4  = dis4*cdis(istage)
        diw(2,:) = b%w(i,j,k,:) -wsol(:)
        diw(1,:) = diw(2,:)
        diw(3,:) = b%w(i+1,j,k,:) -b%w(i,j,k,:)
        fs1(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs1(:) = dis2*diw(2,:) -dis4*fs1(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs1(:) = fs1(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(wsol(:),b%w(i,j,k,:),fs1(:),b%si(i,j,k,:),fs1(:))
        end if

        b%fw(i  ,j,k,:) = b%fw(i  ,j,k,:) +(fs1(:) -fs2(:))
        b%fw(i+1,j,k,:) = b%fw(i+1,j,k,:) +(fs2(:) -dfsx(i+2,j,k,:))
      end if
    else if(zz==2) then
      if(id == 5) then
        b%fw(i,j  ,k,:) = b%fw(i,j  ,k,:) -(dfsy(i,j,k,:) -dfsy(i,j+1,k,:))
        b%fw(i,j-1,k,:) = b%fw(i,j-1,k,:) -(dfsy(i,j-1,k,:) -dfsy(i,j,k,:))
        fp(1) = abs(b%p(i,j,k) -2.0d0*b%p(i,j-1,k) +b%p(i,j-2,k))/&
                abs(b%p(i,j,k) +2.0d0*b%p(i,j-1,k) +b%p(i,j-2,k))
        fp(2) = abs(psol -2.0d0*b%p(i,j,k) +b%p(i,j-1,k))/&
                abs(psol +2.0d0*b%p(i,j,k) +b%p(i,j-1,k))
        ep    = max(fp(1),fp(2))
        lm    = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam   = lm
        dis2     = vis2*ep
        dis4     = vis4
        dis4     = dim(dis4,dis2)
        dis2     = dis2*cdis(istage)
        dis4     = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k,:) -b%w(i,j-1,k,:)
        diw(2,:) = wsol(:) -b%w(i,j,k,:)
        diw(3,:) = diw(2,:)
        fs2(:)   = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs2(:)   = dis2*diw(2,:) -dis4*fs2(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs2(:) = fs2(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k,:),wsol(:),fs2(:),b%sj(i,j+1,k,:),fs2(:))
        end if

!! modified by tsllc 01-12-10
!! set fp(1) to 0 if (j-3) is negative
        if ((j-3)>=0) then
          fp(3) = abs(b%p(i,j-1,k) -2.0d0*b%p(i,j-2,k) +b%p(i,j-3,k))/&
                  abs(b%p(i,j-1,k) +2.0d0*b%p(i,j-2,k) +b%p(i,j-3,k))
        else
          fp(3) = 0.0d0
        end if
        ep = max(fp(1),fp(2),fp(3))
        lm = b%radi(i,j-1,k) +b%radj(i,j-1,k) +b%radk(i,j-1,k)
        lp = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam      = 0.5d0*(lm +lp)
        dis2     = vis2*ep
        dis4     = vis4
        dis4     = dim(dis4,dis2)
        dis2     = dis2*cdis(istage)
        dis4     = dis4*cdis(istage)
        diw(1,:) = b%w(i,j-1,k,:) -b%w(i,j-2,k,:)
        diw(2,:) = b%w(i,j,k,:) -b%w(i,j-1,k,:)
        diw(3,:) = wsol(:) -b%w(i,j,k,:)
        fs1(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs1(:) = dis2*diw(2,:) -dis4*fs1(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs1(:) = fs1(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j-1,k,:),b%w(i,j,k,:),fs1(:),b%sj(i,j,k,:),fs1(:))
        end if

        b%fw(i,j  ,k,:) = b%fw(i,j  ,k,:) +(fs1(:) -fs2(:))
        b%fw(i,j-1,k,:) = b%fw(i,j-1,k,:) +(dfsy(i,j-1,k,:) -fs1(:))
      else if(id==2) then
        b%fw(i,j  ,k,:) = b%fw(i,j  ,k,:) -(dfsy(i,j,k,:) -dfsy(i,j+1,k,:))
        b%fw(i,j+1,k,:) = b%fw(i,j+1,k,:) -(dfsy(i,j+1,k,:) -dfsy(i,j+2,k,:))
        fp(1) = abs(b%p(i,j+1,k) -2.0d0*b%p(i,j,k) +psol)/&
                abs(b%p(i,j+1,k) +2.0d0*b%p(i,j,k) +psol)
        fp(2) = abs(b%p(i,j+2,k) -2.0d0*b%p(i,j+1,k) +b%p(i,j,k))/&
                abs(b%p(i,j+2,k) +2.0d0*b%p(i,j+1,k) +b%p(i,j,k))

!! modified by tsllc 01-12-10
!! set fp(4) to 0 if (j+3) > b%jb
        if ((j+3)<=b%jb) then
          fp(3) = abs(b%p(i,j+3,k) -2.0d0*b%p(i,j+2,k) +b%p(i,j+1,k))/&
                  abs(b%p(i,j+3,k) +2.0d0*b%p(i,j+2,k) +b%p(i,j+1,k)) 
        else
          fp(3) = 0.0d0
        end if

        ep = max(fp(1),fp(2),fp(3))
        lm = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lp = b%radi(i,j+1,k) +b%radj(i,j+1,k) +b%radk(i,j+1,k)
        lam   = 0.5d0*(lm +lp)
        dis2  = vis2*ep
        dis4  = vis4
        dis4  = dim(dis4,dis2)
        dis2  = dis2*cdis(istage)
        dis4  = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k,:) -wsol(:)
        diw(2,:) = b%w(i,j+1,k,:) -b%w(i,j,k,:)
        diw(3,:) = b%w(i,j+2,k,:) -b%w(i,j+1,k,:)
        fs2(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs2(:) = dis2*diw(2,:) -dis4*fs2(:)
        
        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs2(:) = fs2(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k,:),b%w(i,j+1,k,:),fs2(:),b%sj(i,j+1,k,:),fs2(:))
        end if

        ep  = max(fp(1),fp(2))
        lp  = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam = lp
        dis2  = vis2*ep
        dis4  = vis4
        dis4  = dim(dis4,dis2)
        dis2  = dis2*cdis(istage)
        dis4  = dis4*cdis(istage)
        diw(2,:) = b%w(i,j,k,:) -wsol(:)
        diw(1,:) = diw(2,:)
        diw(3,:) = b%w(i,j+1,k,:) -b%w(i,j,k,:)
        fs1(:)   = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs1(:)   = dis2*diw(2,:) -dis4*fs1(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs1(:) = fs1(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(wsol(:),b%w(i,j,k,:),fs1(:),b%sj(i,j,k,:),fs1(:))
        end if

        b%fw(i,j  ,k,:) = b%fw(i,j  ,k,:) +(fs1(:) -fs2(:))
        b%fw(i,j+1,k,:) = b%fw(i,j+1,k,:) +(fs2(:) -dfsy(i,j+2,k,:))
      end if
    else if(zz==3) then
      if(id==6) then
        b%fw(i,j,k  ,:) = b%fw(i,j,k  ,:) -(dfsz(i,j,k,:) -dfsz(i,j,k+1,:))
        b%fw(i,j,k-1,:) = b%fw(i,j,k-1,:) -(dfsz(i,j,k-1,:) -dfsz(i,j,k,:))
        fp(1) = abs(b%p(i,j,k) -2.0d0*b%p(i,j,k-1) +b%p(i,j,k-2))/&
                abs(b%p(i,j,k) +2.0d0*b%p(i,j,k-1) +b%p(i,j,k-2))
        fp(2) = abs(psol -2.0d0*b%p(i,j,k) +b%p(i,j,k-1))/&
                abs(psol +2.0d0*b%p(i,j,k) +b%p(i,j,k-1))
        ep  = max(fp(1),fp(2))
        lm  = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam = lm
        dis2     = vis2*ep
        dis4     = vis4
        dis4     = dim(dis4,dis2)
        dis2     = dis2*cdis(istage)
        dis4     = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k,:) -b%w(i,j,k-1,:)
        diw(2,:) = wsol(:) -b%w(i,j,k,:)
        diw(3,:) = diw(2,:)
        fs2(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs2(:) = dis2*diw(2,:) -dis4*fs2(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs2(:) = fs2(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k,:),wsol(:),fs2(:),b%sk(i,j,k+1,:),fs2(:))
        end if

!! modified by tsllc 01-12-10
!! set fp(1) to 0 if (k-3) is negative
        if ((k-3)>=0) then
          fp(3) = abs(b%p(i,j,k-1) -2.0d0*b%p(i,j,k-2) +b%p(i,j,k-3))/&
                  abs(b%p(i,j,k-1) +2.0d0*b%p(i,j,k-2) +b%p(i,j,k-3))
        else
          fp(3) = 0.0d0
        end if

        ep = max(fp(1),fp(2),fp(3))
        lm = b%radi(i,j,k-1) +b%radj(i,j,k-1) +b%radk(i,j,k-1)
        lp = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam      = 0.5d0*(lm +lp)
        dis2     = vis2*ep
        dis4     = vis4
        dis4     = dim(dis4,dis2)
        dis2     = dis2*cdis(istage)
        dis4     = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k-1,:) -b%w(i,j,k-2,:)
        diw(2,:) = b%w(i,j,k,:) -b%w(i,j,k-1,:)
        diw(3,:) = wsol(:) -b%w(i,j,k,:)
        fs1(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs1(:) = dis2*diw(2,:) -dis4*fs1(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs1(:) = fs1(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k-1,:),b%w(i,j,k,:),fs1(:),b%sk(i,j,k,:),fs1(:))
        end if

        b%fw(i,j,k,  :) = b%fw(i,j,k  ,:) +(fs1(:) -fs2(:))
        b%fw(i,j,k-1,:) = b%fw(i,j,k-1,:) +(dfsz(i,j,k-1,:) -fs1(:))
      else if(id==3) then
        b%fw(i,j,k,  :) = b%fw(i,j,k,  :) -(dfsz(i,j,k,  :) -dfsz(i,j,k+1,:))
        b%fw(i,j,k+1,:) = b%fw(i,j,k+1,:) -(dfsz(i,j,k+1,:) -dfsz(i,j,k+2,:))
        fp(1) = abs(b%p(i,j,k+1) -2.0d0*b%p(i,j,k) +psol)/&
                abs(b%p(i,j,k+1) +2.0d0*b%p(i,j,k) +psol)
        fp(2) = abs(b%p(i,j,k+2) -2.0d0*b%p(i,j,k+1) +b%p(i,j,k))/&
                abs(b%p(i,j,k+2) +2.0d0*b%p(i,j,k+1) +b%p(i,j,k))

!! modified by tsllc 01-12-10
!! set fp(4) to 0 if (k+3) > b%kb
        if ((k+3)<=b%kb) then
          fp(3) = abs(b%p(i,j,k+3) -2.0d0*b%p(i,j,k+2) +b%p(i,j,k+1))/&
                  abs(b%p(i,j,k+3) +2.0d0*b%p(i,j,k+2) +b%p(i,j,k+1)) 
        else
          fp(3) = 0.0d0
        end if

        ep = max(fp(1),fp(2),fp(3))
        lm = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lp = b%radi(i,j,k+1) +b%radj(i,j,k+1) +b%radk(i,j,k+1)
        lam   = 0.5d0*(lm +lp)
        dis2  = vis2*ep
        dis4  = vis4
        dis4  = dim(dis4,dis2)
        dis2  = dis2*cdis(istage)
        dis4  = dis4*cdis(istage)
        diw(1,:) = b%w(i,j,k,:) -wsol(:)
        diw(2,:) = b%w(i,j,k+1,:) -b%w(i,j,k,:)
        diw(3,:) = b%w(i,j,k+2,:) -b%w(i,j,k+1,:)
        fs2(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs2(:) = dis2*diw(2,:) -dis4*fs2(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs2(:) = fs2(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k,:),b%w(i,j,k+1,:),fs2(:),b%sk(i,j,k+1,:),fs2(:))
        end if

        ep  = max(fp(1),fp(2))
        lp  = b%radi(i,j,k) +b%radj(i,j,k) +b%radk(i,j,k)
        lam = lp
        dis2  = vis2*ep
        dis4  = vis4
        dis4  = dim(dis4,dis2)
        dis2  = dis2*cdis(istage)
        dis4  = dis4*cdis(istage)
        diw(2,:) = b%w(i,j,k,:) -wsol(:)
        diw(1,:) = diw(2,:)
        diw(3,:) = b%w(i,j,k+1,:) -b%w(i,j,k,:)
        fs1(:) = diw(3,:) -2.0d0*diw(2,:) +diw(1,:)
        fs1(:) = dis2*diw(2,:) -dis4*fs1(:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          fs1(:) = fs1(:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(wsol(:),b%w(i,j,k,:),fs1(:),b%sk(i,j,k,:),fs1(:))
        end if

        b%fw(i,j,k  ,:) = b%fw(i,j,k  ,:) +(fs1(:) -fs2(:))
        b%fw(i,j,k+1,:) = b%fw(i,j,k+1,:) +(fs2(:) -dfsz(i,j,k+2,:))
      end if
    end if
  end do

  return
end subroutine dfluxb

