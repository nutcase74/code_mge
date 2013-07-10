subroutine upload_solut_xcross(b, ptr, p_tot, atcellcenter)
  use typepara
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(in)      :: b
  type(mgrid), intent(inout)   :: ptr
  real(dp), intent(in)         :: p_tot 
  integer, intent(inout)       :: atcellcenter

  integer :: i,j,k, nval
  real(dp) :: rr,uu,vv,ww,re,qq,pp,mm,ptloss
  integer, dimension(3) :: lowc, uppc, ie
!!
!! SET FLAG FOR VALUE AT CELL CENTER 
!!
  atcellcenter = 1

  lowc = ptr%lowcell
  uppc = ptr%uppcell

  nval = 7

  ptr%up2date  = 1

  if (ptr%isolut/=7) then
    ptr%isolut = 7
    nval = 7

    if (associated(ptr%varname)) deallocate(ptr%varname)
    if (associated(ptr%solut))   deallocate(ptr%solut)

    allocate( ptr%solut(lowc(1):uppc(1),lowc(2):uppc(2),lowc(3):uppc(3), nval ),&
              ptr%varname(nval) )
  end if 

  if (.false.) then
  print *," ******************************************************************************"
  print *," ptr%LEVEL, SERIAL, NCEL ",ptr%level, ptr%serial, ptr%ncell
  print *," ISOLUT ", ptr%isolut
  print *," ASSOCIATED VARNAME ", associated(ptr%varname)
  print *," ASSOCIATED SOLUT   ", associated(ptr%solut)
  print *," LBOUND SOLUT       ", lbound(ptr%solut,1), lbound(ptr%solut,2), lbound(ptr%solut,3)
  print *," UBOUND SOLUT P     ", ubound(ptr%solut,1), ubound(ptr%solut,2), ubound(ptr%solut,3)
  print *," ******************************************************************************"
  end if

  ptr%varname(1) =  'rho'
  ptr%varname(2) =  'u'
  ptr%varname(3) =  'v'
  ptr%varname(4) =  'w'
  ptr%varname(5) =  'cp'
  ptr%varname(6) =  'mach'
  ptr%varname(7) =  'pt-loss'

  do k = 2, b%kl
    do j = 2, b%jl
      do i = 2, b%il
        ! "variables=x,y,z,r,u,v,w,cp,mach,ptloss"
        rr = b%w(i,j,k,1)
        uu = b%w(i,j,k,2)/rr
        vv = b%w(i,j,k,3)/rr
        ww = b%w(i,j,k,4)/rr
        re = b%w(i,j,k,5)
        qq = 0.5d0*(uu**2 + vv**2 + ww**2)*rr
        pp = (gamma -1.0d0)*dim(re,qq)
        mm = sqrt((uu*uu+vv*vv+ww*ww)/gamma/pp*rr)
        ptloss = pp*(1.0d0+0.2d0*mm*mm)**3.5d0/p_tot - 1.0d0

        pp = (pp/p0 - 1.0d0)*2.0d0/(gamma*rm*rm)

        !! CELL INDEX IN COMPUTATONAL MESH 
        ie = (/i,j,k/)
        ptr%solut(ie(1),ie(2),ie(3),1) = rr
        ptr%solut(ie(1),ie(2),ie(3),2) = uu
        ptr%solut(ie(1),ie(2),ie(3),3) = vv
        ptr%solut(ie(1),ie(2),ie(3),4) = ww
        ptr%solut(ie(1),ie(2),ie(3),5) = pp
        ptr%solut(ie(1),ie(2),ie(3),6) = mm
        ptr%solut(ie(1),ie(2),ie(3),7) = ptloss

      end do
    end do
  end do
  return
end subroutine upload_solut_xcross


subroutine upload_solut(blk, p_tot, atcellcenter)
  use typepara
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), pointer         :: blk
  real(dp), intent(in)         :: p_tot 
  integer, intent(inout)       :: atcellcenter

  integer :: i,j,k, nval
  real(dp) :: rr,uu,vv,ww,re,qq,pp,mm,ptloss
  integer, dimension(3) :: lowc, uppc, ie
  type(mgrid), pointer  :: ptr

  integer, dimension(3) :: icut,isol,isol0
  integer :: ctr,jx,id,id0,m
  real(dp) :: ptmp
  real(dp), dimension(5) :: wtmp
!!
!! SET FLAG FOR VALUE AT CELL CENTER 
!!
  atcellcenter = 1

  ptr=>blk
  do while (associated(ptr))

    lowc = ptr%lowcell
    uppc = ptr%uppcell

    ptr%up2date  = 1

    if(eqntype == 'e' .or. eqntype == 'E') then
      nval = 7
      if (ptr%isolut/=7) then
        ptr%isolut = 7
        nval = 7

        if (associated(ptr%varname)) deallocate(ptr%varname)
        if (associated(ptr%solut))   deallocate(ptr%solut)

        allocate( ptr%solut(lowc(1):uppc(1),lowc(2):uppc(2),lowc(3):uppc(3), nval ),&
                  ptr%varname(nval) )
      end if
    else
      nval = 11
      if (ptr%isolut/=11) then
        ptr%isolut = 11
        nval = 11

        if (associated(ptr%varname)) deallocate(ptr%varname)
        if (associated(ptr%solut))   deallocate(ptr%solut)

        allocate( ptr%solut(lowc(1):uppc(1),lowc(2):uppc(2),lowc(3):uppc(3), nval ),&
                  ptr%varname(nval) )
      end if
 
    end if 

    if (.false.) then
    print *," ******************************************************************************"
    print *," ptr%LEVEL, SERIAL, NCEL ",ptr%level, ptr%serial, ptr%ncell
    print *," ISOLUT ", ptr%isolut
    print *," ASSOCIATED VARNAME ", associated(ptr%varname)
    print *," ASSOCIATED SOLUT   ", associated(ptr%solut)
    print *," LBOUND SOLUT       ", lbound(ptr%solut,1), lbound(ptr%solut,2), lbound(ptr%solut,3)
    print *," UBOUND SOLUT P     ", ubound(ptr%solut,1), ubound(ptr%solut,2), ubound(ptr%solut,3)
    print *," ******************************************************************************"
    end if

    ptr%varname(1) =  'rho'
    ptr%varname(2) =  'u'
    ptr%varname(3) =  'v'
    ptr%varname(4) =  'w'
    ptr%varname(5) =  'cp'
    ptr%varname(6) =  'mach'
    ptr%varname(7) =  'pt-loss'
    if(eqntype == 'n' .or. eqntype == 'N') then
      ptr%varname(8) =  'energy'
      ptr%varname(9)  =  'gradn1'
      ptr%varname(10) =  'gradn2'
      ptr%varname(11) =  'gradn3'
    end if

if(.false.) then
    do m = 1,ptr%nmir
      icut = ptr%z(m)%ccid
      isol = ptr%z(m)%scid
      ! set value on the solid cell next to the cut cell
      ! take average of wsol from the cut cell neighbors
      if(ptr%blank(isol(1),isol(2),isol(3)) /= 0) cycle
      
      wtmp = ptr%z(m)%wsol
      ptmp = ptr%z(m)%psol
      ctr  = 1
      do jx = 1,m-1
        isol0 = ptr%z(jx)%scid
        if(all(isol==isol0)) then
          ctr = ctr + 1
          wtmp = wtmp + ptr%z(jx)%wsol
          ptmp = ptmp + ptr%z(jx)%psol
        end if
      end do
      if(ctr > 1) then
        wtmp = wtmp/ctr
        ptmp = ptmp/ctr
      end if
      ptr%w(isol(1),isol(2),isol(3),:) = wtmp(:)
      ptr%p(isol(1),isol(2),isol(3))   = ptmp
      !ptr%w(isol(1),isol(2),isol(3),:) = ptr%z(m)%wsol(:)
      !ptr%p(isol(1),isol(2),isol(3))   = ptr%z(m)%psol

    end do
end if

    do k = 2, ptr%kl
      do j = 2, ptr%jl
        do i = 2, ptr%il

          ! "variables=x,y,z,r,u,v,w,cp,mach,ptloss"
          rr = ptr%w(i,j,k,1)
          uu = ptr%w(i,j,k,2)/rr
          vv = ptr%w(i,j,k,3)/rr
          ww = ptr%w(i,j,k,4)/rr
          re = ptr%w(i,j,k,5)
          qq = 0.5d0*(uu**2 + vv**2 + ww**2)*rr
          pp = (gamma -1.0d0)*dim(re,qq)
          mm = sqrt((uu*uu+vv*vv+ww*ww)/gamma/pp*rr)
          ptloss = pp*(1.0d0+0.2d0*mm*mm)**3.5d0/p_tot - 1.0d0

          pp = (pp/p0 - 1.0d0)*2.0d0/(gamma*rm*rm)

          !! CELL INDEX IN COMPUTATONAL MESH 
          ie = (/i,j,k/)
          ptr%solut(ie(1),ie(2),ie(3),1) = rr
          ptr%solut(ie(1),ie(2),ie(3),2) = uu
          ptr%solut(ie(1),ie(2),ie(3),3) = vv
          ptr%solut(ie(1),ie(2),ie(3),4) = ww
          ptr%solut(ie(1),ie(2),ie(3),5) = pp
          ptr%solut(ie(1),ie(2),ie(3),6) = mm
          ptr%solut(ie(1),ie(2),ie(3),7) = ptloss
          if(eqntype == 'n' .or. eqntype == 'N') then
            ptr%solut(ie(1),ie(2),ie(3),8) = re/rr
            ptr%solut(ie(1),ie(2),ie(3),9)  = 0.0d0
            ptr%solut(ie(1),ie(2),ie(3),10) = 0.0d0
            ptr%solut(ie(1),ie(2),ie(3),11) = 0.0d0
          end if
        end do
      end do
    end do
    ptr=>ptr%nxt
  end do
  return
end subroutine upload_solut
