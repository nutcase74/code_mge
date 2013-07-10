subroutine dflux(b,istage)
! calculate artificial dissipation
! use Jameson's mixed first order and third order difference formula
! ==================================================================
!  use typeblock
  use typelib4
  use typepara
  implicit none
  integer, parameter :: dp=kind(1.0d0)
! input/output variables
  type(mgrid), intent(inout) :: b
  integer, intent(in) :: istage
! variables
  integer  :: i, j, k, n, mx
  integer :: il,jl,kl,ie,je,ke,ib,jb,kb
  real(dp) :: lm, lp, lam, dis2, dis4, ep, pb
  real(dp), dimension(:), pointer :: fp
  real(dp), dimension(:,:), pointer :: diw
  real(dp), dimension(:,:,:,:), pointer :: dfsx,dfsy,dfsz

  real(dp) :: rr, ee, qq, c2, hh, cs, qn
  real(dp) :: e0, e1, e2, e3, s1, s2, s3
  real(dp) :: vd, qd, ak1, ak2, ak3, gm1
  real(dp) :: epsn, epsl
  real(dp), dimension(3) :: us, ss
  real(dp), dimension(5) :: fs
  real(dp), dimension(:,:,:), allocatable :: we
! ******************************************************************
  if(cdis(istage) .le. 0.0) return

  gm1 = gamma - 1.0d0

  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb

  mx = max(ib,max(jb,kb))
  allocate(fp(0:mx), diw(0:mx,5))
  allocate(dfsx(0:ib,0:jb,0:kb,5))
  allocate(dfsy(0:ib,0:jb,0:kb,5),dfsz(0:ib,0:jb,0:kb,5))
  allocate(we(0:ib,0:jb,0:kb))

  dfsx = 0.0d0; dfsy = 0.0d0; dfsz = 0.0d0

  fp   = 0.0d0
  ep   = 0.0d0
  lam  = 0.0d0
  dis2 = 0.0d0
  dis4 = 0.0d0
  diw  = 0.0d0
  !!!   REPLACE  energy by enthalpy
  do k=0, kb
    do j=0, jb
      do i=0, ib
        we( i,j,k)   = b%w(i,j,k,5)
        b%w(i,j,k,5) = b%w(i,j,k,5) + b%p(i,j,k)
      end do
    end do
  end do

 ! dissipation in i-direction
  do k=2,kl
    do j=2,jl
      do i=1,ie
        pb = abs(b%p(i-1,j,k) +2.0d0*b%p(i,j,k) +b%p(i+1,j,k))
        if(abs(pb) .gt. 0.0d0) then
          fp(i) = abs(b%p(i-1,j,k) -2.0d0*b%p(i,j,k) +b%p(i+1,j,k))/pb
        else
          fp(i) = 0.0d0
        end if
      end do
      fp(0) = 0.0d0
      !fp(1) = 0.0d0
      !fp(ie) = 0.0d0
      fp(ib) = 0.0d0
      do i=2,ie
        diw(i,:) = b%w(i,j,k,:) -b%w(i-1,j,k,:)
      end do
      diw(1, :) = diw(2,:)
      diw(ib,:) = diw(ie,:)
      do i=2,ie
        ep   = max(fp(i-2),fp(i-1),fp(i),fp(i+1))
        lp   = b%radi(i,j,k)+b%radj(i,j,k)+b%radk(i,j,k)
        lm   = b%radi(i-1,j,k) +b%radj(i-1,j,k) +b%radk(i-1,j,k)
        lam  = 0.5d0*(lm +lp)
 	dis2 = vis2*ep
        dis4 = vis4
        dis4 = dim(dis4,dis2)
        dis2 = dis2*cdis(istage)
        dis4 = dis4*cdis(istage)
        dfsx(i,j,k,:) = diw(i+1,:) -2.0d0*diw(i,:) +diw(i-1,:)
        dfsx(i,j,k,:) = dis2*diw(i,:) -dis4*dfsx(i,j,k,:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          dfsx(i,j,k,:) = dfsx(i,j,k,:) * lam
          ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i-1,j,k,:),b%w(i,j,k,:), &
                            dfsx(i,j,k,:),b%si(i,j,k,:),fs(:))
          dfsx(i,j,k,:)= fs(:)
        else
          print *, "dissip must be 'scalar' or 'matrix', dissip=", dissip
          stop
        end if

      end do
    end do
  end do

  do n=1,5
    do k=2,kl
      do j=2,jl
        do i=2,il
          b%fw(i,j,k,n) = b%fw(i,j,k,n) + (dfsx(i,j,k,n) -dfsx(i+1,j,k,n))
        end do
      end do
    end do
  end do

! dissipation in j-direction
  if(i2d==0) then
  do k=2,kl
    do i=2,il
      do j=1,je
        pb = abs(b%p(i,j+1,k) +2.0d0*b%p(i,j,k) +b%p(i,j-1,k))
        if(abs(pb) .gt. 0.0d0) then
          fp(j) = abs(b%p(i,j+1,k) -2.0d0*b%p(i,j,k) +b%p(i,j-1,k))/pb
        else 
          fp(j) = 0.0d0
        end if
      end do
      fp(0) = 0.0d0
      !fp(1) = 0.0d0
      !fp(je) = 0.0d0
      fp(jb) = 0.0d0
      do j=2,je
        diw(j,:) = b%w(i,j,k,:) -b%w(i,j-1,k,:)
      end do
      diw(1,:)  = diw(2,:)
      diw(jb,:) = diw(je,:)
      do j=2,je
        ep   = max(fp(j-2),fp(j-1),fp(j),fp(j+1))
        lp   = b%radi(i,j,  k) +b%radj(i,j,  k) +b%radk(i,j,  k)
        lm   = b%radi(i,j-1,k) +b%radj(i,j-1,k) +b%radk(i,j-1,k)
        lam  = 0.5d0*(lm +lp)
        dis2 = vis2*ep
        dis4 = vis4
        dis4 = dim(dis4,dis2)
        dis2 = dis2*cdis(istage)
        dis4 = dis4*cdis(istage)
        dfsy(i,j,k,:) = diw(j+1,:) -2.0d0*diw(j,:) +diw(j-1,:)
        dfsy(i,j,k,:) = dis2*diw(j,:) -dis4*dfsy(i,j,k,:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          dfsy(i,j,k,:)   = dfsy(i,j,k,:) * lam
         ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j-1,k,:),b%w(i,j,k,:), &
                            dfsy(i,j,k,:),b%sj(i,j,k,:),fs(:))
          dfsy(i,j,k,:)= fs(:)
        end if

      end do
    end do
  end do
  end if

  do n=1,5
    do k=2,kl
      do j=2,jl
        do i=2,il
          b%fw(i,j,k,n) = b%fw(i,j,k,n) +(dfsy(i,j,k,n) -dfsy(i,j+1,k,n))
        end do
      end do
    end do
  end do
 
! dissipation in k-direction
  do i=2,il
    do j=2,jl
      do k=1,ke
        pb = abs(b%p(i,j,k-1) +2.0d0*b%p(i,j,k) +b%p(i,j,k+1))
        if(abs(pb) .gt. 0.0d0) then
          fp(k) = abs(b%p(i,j,k-1) -2.0d0*b%p(i,j,k) +b%p(i,j,k+1))/pb
        else 
          fp(k) = 0.0d0
        end if
      end do
      fp(0) = 0.0d0
      !fp(1) = 0.0d0
      !fp(ke) = 0.0d0
      fp(kb) = 0.0d0
      do k=2,ke
        diw(k,:) = b%w(i,j,k,:) -b%w(i,j,k-1,:)
      end do
      diw(1, :) = diw(2,:)
      diw(kb,:) = diw(ke,:)
      do k=2,ke
        ep   = max(fp(k-2),fp(k-1),fp(k), fp(k+1))
        lp   = b%radi(i,j,k  ) +b%radj(i,j,k  ) +b%radk(i,j,k  )
        lm   = b%radi(i,j,k-1) +b%radj(i,j,k-1) +b%radk(i,j,k-1)
        lam  = 0.5d0*(lm +lp)
        dis2 = vis2*ep
        dis4 = vis4
        dis4 = dim(dis4,dis2)
        dis2 = dis2*cdis(istage)
        dis4 = dis4*cdis(istage)
        dfsz(i,j,k,:) = diw(k+1,:) -2.0d0*diw(k,:) +diw(k-1,:)
        dfsz(i,j,k,:) = dis2*diw(k,:) -dis4*dfsz(i,j,k,:)

        if(dissip == 'scalar' .or. dissip == 'SCALAR') then
          dfsz(i,j,k,:)   = dfsz(i,j,k,:) * lam
        ! Matrix dissipation
        else if(dissip == 'matrix' .or. dissip == 'MATRIX') then
          call dflux_matrix(b%w(i,j,k-1,:),b%w(i,j,k,:), &
                            dfsz(i,j,k,:),b%sk(i,j,k,:),fs(:))
          dfsz(i,j,k,:)= fs(:)
        end if

      end do
    end do
  end do

  do n=1,5
    do k=2,kl
      do j=2,jl
        do i=2,il
          b%fw(i,j,k,n) = b%fw(i,j,k,n) +(dfsz(i,j,k,n) -dfsz(i,j,k+1,n))
        end do
      end do
    end do
  end do

   !debug
   !do k = 2,kl
   !  do j = 2,jl
   !    do i = 2,il
   !      if (b%fw(i,j,k,1) .ne. b%fw(i,j,k,1) ) then
   !        print *,"dflux1 ijk",i,j,k,"fw",b%fw(i,j,k,1),b%w(i,j,k,1)
   !        stop
   !      end if
   !    end do
   !  end do
   !end do

! --- ARTIFICIAL DISSIPATION AT SOLID WALL, I.E. CUT-CELLS
  if(.not. ibm_on) call dfluxb(b,istage,dfsx,dfsy,dfsz)
!
!     replace the enthalpy by the energy
!
  do k=0,kb
    do j=0,jb
      do i=0,ib
        b%w(i,j,k,5)  = we(i,j,k)
      enddo
    enddo
  enddo

  deallocate(fp,diw,we)
  deallocate(dfsx,dfsy,dfsz)

   !debug
   !do k = 2,kl
   !  do j = 2,jl
   !    do i = 2,il
   !      if (b%fw(i,j,k,1) .ne. b%fw(i,j,k,1) ) then
   !        print *,"dflux2 ijk",i,j,k,"fw",b%fw(i,j,k,1),b%w(i,j,k,1)
   !        stop
   !      end if
   !    end do
   !  end do
   !end do

  return
end subroutine dflux

subroutine dflux_matrix(u1,u2,dw,ss,ff)
!
! CALCULATE THE MATRIX IN MATRIX DISSIPATION
!
  use typepara
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  real(dp), dimension(5), intent(in) :: u1,u2,dw
  real(dp), dimension(3), intent(in) :: ss
  real(dp), dimension(5), intent(out) :: ff
  
  real(dp) :: us(3),gam1,epsn,epsl
  real(dp) :: rr,ee,qq,c2,hh,cs,qn
  real(dp) :: e0,e1,e2,e3,s1,s2,s3
  real(dp) :: ak1,ak2,ak3,vd,qd

  gam1 = gamma - 1.0d0
  epsn = 1.0d0
  epsl = 1.0d0
  !epsn = 0.25d0
  !epsl = 0.025d0

  rr      = 0.5d0*(u1(1)  +u2(1))
  us(1:3) = 0.5d0*(u1(2:4)+u2(2:4)) / rr
  ee      = 0.5d0*(u1(5)  +u2(5))/rr
  
  qq = 0.5d0*(sum(us*us))
  c2 = gamma*gam1*(ee-qq)
  hh = c2/gam1 + qq
  
  cs = sqrt(c2*sum(ss*ss))
  qn = sum(ss*us)
  
  e0 = abs(qn) + cs
  e1 = max(epsn*e0, abs(qn+cs))
  e2 = max(epsn*e0, abs(qn-cs))
  e3 = max(epsl*e0, abs(qn   ))
  
  s1 = 0.5d0*(e1+e2)
  s2 = 0.5d0*(e1-e2)
  s3 = e3
  
  ak1   = (s1-s3)/c2*gam1
  ak2   = s2/cs
  ak3   = c2/gam1/sum(ss*ss)
  
  vd = qq*dw(1) - us(1)*dw(2) - us(2)*dw(3) - us(3)*dw(4) + dw(5)
  qd = qn*dw(1) - ss(1)*dw(2) - ss(2)*dw(3) - ss(3)*dw(4)

  ff(1)= ak1*       vd              -ak2*qd                      +s3*dw(1)
  ff(2)= ak1*(us(1)*vd-ss(1)*ak3*qd)+ak2*(gam1*ss(1)*vd-us(1)*qd)+s3*dw(2)
  ff(3)= ak1*(us(2)*vd-ss(2)*ak3*qd)+ak2*(gam1*ss(2)*vd-us(2)*qd)+s3*dw(3)
  ff(4)= ak1*(us(3)*vd-ss(3)*ak3*qd)+ak2*(gam1*ss(3)*vd-us(3)*qd)+s3*dw(4)
  ff(5)= ak1*(   hh*vd-   qn*ak3*qd)+ak2*(gam1*qn*vd-hh*qd)      +s3*dw(5)

  return
end subroutine dflux_matrix
