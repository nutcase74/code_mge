subroutine bc_cutcell(b)

!****************************************************************************
!
! INTERPOLATE STATE VARIABLES AT SOLID CELLS
! USING CURVATURE-CORRECTED SYMMETRY TECHNIQUE
!
!****************************************************************************


  use typelib
  use typepara
  use ctrlib

  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b
!
! local variables
!
  integer :: i , j, k, m, ncl, n, zz, id
  integer :: isol(3), isol0(3), mm, ctr, jx, icut(3)
  integer, dimension(3) :: jxa, di
  integer, dimension(:,:), pointer :: cloud
  real(dp) :: nn,cc,dpdn ,dpdn1,dpdn2, qn,dx
  real(dp) :: psol, p_refl, dis, velnor, veltan, ptmp
  real(dp), dimension(:), pointer :: vphi
  real(dp), dimension(2) :: curv
  real(dp), dimension(3) :: nvec, vel_tanp, tan_plan, ppt, qb, qt, nnew
  real(dp), dimension(3) :: xs,xc,us,uc
  real(dp), dimension(5) :: var_refl, wsol, wtmp
  real(dp), dimension(:,:), pointer :: var 
  real(dp), dimension(:), pointer :: varp
! *****************************************************************************

  !print *," BCUT SAY size(b%z) ", size(b%z), " NMIR=", b%nmir

  dx = b%cellvolume**(1.0/3.0d0)

  do m = 1,b%nmir
    ncl  = b%z(m)%ipol%npoint
    zz   = b%z(m)%zz
    isol = b%z(m)%scid
    icut = b%z(m)%ccid
    id   = b%z(m)%faceid

    if (approxbc) then
       ppt  = b%z(m)%pnew
       nvec = b%z(m)%norv
       qb   = b%z(m)%vel/vel_fs   ! normalize the boundary velocity by freestream velocity
    else
       ppt  = b%z(m)%psurf
       nvec = b%z(m)%nvec
       qb   = 0.0d0 
    end if

    nn      = sqrt(sum(nvec*nvec))
    nvec    = nvec/nn
    curv    = b%z(m)%radcur
    dis     = 2.0d0*b%z(m)%rr

    psol      = 0.0d0
    wsol(1:5) = 0.0d0

    if(ncl == 0) then
       print *, "Warning: ZERO cloud points found !!!"
       b%z(m)%psol = psol
       b%z(m)%wsol = wsol
       return
    end if

    allocate(cloud(ncl,3),vphi(1:ncl),var(ncl,5),varp(ncl))
    cloud = 0
    vphi = 0.0d0
    var  = 0.0d0
    varp = 0.0d0
    
    cloud = b%z(m)%ipol%ijkdex
    vphi  = b%z(m)%ipol%coef

    do i=1,ncl
       var(i,:) =b%w(cloud(i,1),cloud(i,2),cloud(i,3),:) 
       varp(i)  =b%p(cloud(i,1),cloud(i,2),cloud(i,3))
    end do

! compute the values on the reflected node
    var_refl = 0.0d0
    do j=1,5
       do i=1,ncl
          var_refl(j)  = var_refl(j) +vphi(i)*var(i,j)
       end do
    end do
    
    p_refl = 0.0d0
    do i=1,ncl
       p_refl = p_refl +vphi(i)*varp(i)
    end do

    if (eqntype == 'e' .or. eqntype == 'E') then
!------------------------------------------------------------------------
!
! Compute Solid cell value Using curvature-corrected symmetry technique
!
!------------------------------------------------------------------------
! values at solid cell center

      psol = p_refl
      wsol = var_refl

      if (.false.) then
        ! correction of the pressure
        vel_tanp(1:3) = (var_refl(2:4)-sum(var_refl(2:4)*nvec(1:3))*nvec(1:3))/var_refl(1)
        tan_plan(1:3) = vel_tanp(1:3)/sqrt(sum(vel_tanp(1:3)*vel_tanp(1:3)))
        dpdn = sum(vel_tanp(1:3)**2.0d0)*sqrt(curv(1)**2+curv(2)**2)*wsol(1)
        !dpdn = min(psol,dpdn*dis)
        psol = psol - dpdn*dis

        !correction of density
        !if(psol/p_refl < 0.0d0) then
        !   print *,m,ppt,psol
        !   wsol(1)   = var_refl(1) 
        !else 
        !   wsol(1)   = var_refl(1) *(psol/p_refl)**(1.0d0/gamma) 
        !end if
        !if(psol/p_refl < 0.d0)  write(*, '(I4,12F12.7)') m,ppt(1:3),var_refl(2:4),wsol(2:4),psol,p_refl
        
        !correction of velocity
        if(.false.) then
          cc = 2.0d0*gamma/(gamma-1.0d0)
          cc = cc *(p_refl/var_refl(1) - psol/wsol(1))
          veltan = sqrt(sum(vel_tanp(1:3)**2.0d0)+cc)
          velnor = -sum(var_refl(2:4)*nvec(1:3))/var_refl(1)
          wsol(2:4) = (veltan*tan_plan(1:3) + velnor*nvec(1:3))*wsol(1)
        end if
      end if 
    
      !psol = max(psol, 0.0d0)

      qn   = 2.0d0*(nvec(1)*(wsol(2)-wsol(1)*qb(1)) &
                  + nvec(2)*(wsol(3)-wsol(1)*qb(2)) &
                  + nvec(3)*(wsol(4)-wsol(1)*qb(3)) )

      wsol(2:4) = wsol(2:4)-qn*nvec(1:3)
      wsol(5)   = psol/(gamma-1.0d0)+0.5d0*sum(wsol(2:4)**2.0d0)/wsol(1)
      !write(*, '(I4,13F12.6)') m,ppt(1:3),var_refl(2:4),qb(1:3),wsol(2:4),qb(1)/wsol(2)*100.d0

      do n=1,5
        if(wsol(n) .ne. wsol(n)) then
          print *, "wsol is NaN : ", wsol
          print *, m,psol,dpdn*dis,ppt,curv,var_refl(1)
          print *, "number of cloud points :", ncl
          do i=1,ncl
            write(*,'("cloud :",3I4,7F12.5)') cloud(i,1:3),vphi(i), var(i,1:5), varp(i)
          end do
          stop
        end if
      end do

!--------------------------------------------------------------
! no slip boundary condition for Navier-Stokes equation
! use linear interpolation to get the flow variable
! from the reflected node and wall values 
!--------------------------------------------------------------
    else if(eqntype == 'n' .or. eqntype == 'N') then
      psol = p_refl
      wsol = var_refl

      !------------------------------------------------------------
      !option A: mirror
      !------------------------------------------------------------
      !qn   = 2.0d0*(nvec(1)*(wsol(2)-wsol(1)*qb(1)) &
      !            + nvec(2)*(wsol(3)-wsol(1)*qb(2)) &
      !            + nvec(3)*(wsol(4)-wsol(1)*qb(3)) )
      !wsol(2:4) = wsol(2:4)-qn*nvec(1:3)
      ! keep the normal value unchanged, pointing inward (-nvec)
      ! reverse the tangential direction.
      !nnew(1:3) = -nvec(1:3)
      !qn        = sum(wsol(2:4)*nnew(1:3))
      !qt(1:3)   = wsol(2:4) - qn*nnew(1:3)
      !wsol(2:4) = qn*nnew(1:3)-qt(1:3) 
      !print *, real(wsol(2:4))
      !------------------------------------------------------------
      !option B: linear interpolation 
      !(reflected nodes and surface point at the intersection)
      wsol(2:4) = 2.0d0*qb(1:3)*var_refl(1) - var_refl(2:4)
      !print *, real(wsol(2:4))
      !------------------------------------------------------------
      !option C: linear interpolation 
      !(cut cell and surface point at the intersection ppt)
      !------------------------------------------------------------
      !di = isol-icut
      !xc = b%x(icut(1)-di(1),icut(2)-di(2),icut(3)-di(3),:)+dx/2.0d0
      !xs = b%x(isol(1),isol(2),isol(3),:)+dx/2.0d0
      
      !do i=1,3
      !  if(icut(i)-di(i) /= isol(i)) exit
      !end do   
   !   print *,i,dx,real(xs(i)-xc(i)),real(ppt(i)-xc(i))

      !uc = b%w(icut(1)-di(1),icut(2)-di(2),icut(3)-di(3),2:4)/b%w(icut(1),icut(2),icut(3),1)
      !us = uc+(xs(i)-xc(i))/(ppt(i)-xc(i))*(qb-uc) 
   !   print *, icut,isol
   !   print *, real(xs),real(ppt),real(xc)

      !wsol(2:4) = us*var_refl(1)
      !------------------------------------------------------------
   !   print *, real(wsol(2:4)),real(qb),real(uc)
   !   print *
   !   vel_tanp(1:3) = (var_refl(2:4)-sum(var_refl(2:4)*nvec(1:3))*nvec(1:3))/var_refl(1)
   !   dpdn = sum(vel_tanp(1:3)**2.0d0)*sqrt(curv(1)**2+curv(2)**2)*wsol(1)
   !   print *, real(curv), real(dis), real(psol),real(dpdn*dis)
   !   psol = psol - dpdn*dis

      wsol(5)   = psol/(gamma-1.0d0)+0.5d0*sum(wsol(2:4)**2.0d0)/wsol(1)

    end if

    b%z(m)%psol = psol
    b%z(m)%wsol = wsol

      ! if solid cell has more than one cut cell neighbor
      ! take average of wsol computed from all the cut cell neighbors

    if(.false.) then
      ctr      = 1
      jxa(ctr) = m

        wtmp = wsol
        ptmp = psol
        do jx = 1,m-1
          isol0 = b%z(jx)%scid
          if(all(isol==isol0)) then
            ctr = ctr + 1
            wtmp = wtmp + b%z(jx)%wsol
            ptmp = ptmp + b%z(jx)%psol
            jxa(ctr) = jx
          end if
        end do
        if(ctr > 1) then
          wsol = wtmp/ctr
          psol = ptmp/ctr
          wsol(5) = psol/(gamma-1.0d0)+0.5d0*sum(wsol(2:4)**2.0d0)/wsol(1)
        end if

      do i=1,ctr
        b%z(jxa(i))%psol = psol
        b%z(jxa(i))%wsol = wsol
      end do
    end if

    !b%z(m)%psol = 0.0d0
    !b%z(m)%wsol(2:4) = 0.0d0
    !b%p(isol(1),isol(2),isol(3))   = psol
    !if(b%blank(isol(1),isol(2),isol(3)) == 0) then
    !   b%w(isol(1),isol(2),isol(3),1)   = rho0
    !   b%w(isol(1),isol(2),isol(3),2:4) = 0.0d0
    !   b%w(isol(1),isol(2),isol(3),5)   = rho0*gamma*ei0-p0
    !   call depvars(b%w(isol(1),isol(2),isol(3),:), b%p(isol(1),isol(2),isol(3)))
    !end if

    !b%z(m)%wsol(1)   = rho0
    !b%z(m)%wsol(2:4) = 0.0
    !b%z(m)%wsol(5)   = rho0*gamma*ei0-p0
    !call depvars(b%z(m)%wsol(:),b%z(m)%psol)

    deallocate(cloud, vphi, var, varp)

  end do


  return
end subroutine bc_cutcell

