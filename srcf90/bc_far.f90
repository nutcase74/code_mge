subroutine bc_far(b,face)
! subroutine treating far field boundary conditions for the 
! most coare background mesh, with mesh level==1.
! six faces are treated independently, use 1-D Riemann invariables.
! I=1, inlet bc; I=IL+1, outlet bc;
! J=1, symmetric boundary condition;
! For 2-D test case, J=JL+1, symmetric bc is used.
! ===================================================================
  use typelib4
  use typepara
  implicit none
  integer, parameter :: dp=kind(1.0d0)
! input/output variables
  type(mgrid), intent(inout) :: b
  integer, intent(in) :: face
! local variables
  integer  :: i, j, k,phy1,dum1,dum2
  integer :: il,jl,kl,ie,je,ke,ib,jb,kb
  real(dp) :: rpo, rmi                                   ! Riemann invariables
  real(dp) :: rhoi, ui, vi, wi, prei, ci, si             ! variables interpoated from inner field
  real(dp) :: rhoe, ue, ve, we, pe, ce, se               ! variables at the exit
  real(dp) :: gm, cc, vel_norm
  real(dp) :: fnvec(3)                                   ! face normal vector
! *******************************************************************

  gm = 1.0d0/(gamma -1.0d0)

  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb

! for inlet bc, use subsonic inlet boundary condition
! for outlet bc, use Riemann invariables to treat boundary conditions

! on the left and right boundary
  if(face == 1 .or. face == 4) then
     ! left boundary
     if(face == 1) then
        phy1 = 2
        dum1 = 1
        dum2 = 0
     ! right boundary
     else if(face == 4) then
        phy1 = il
        dum1 = ie
        dum2 = ib
     end if
     fnvec(:) = b%facenvec(face,:)
     do k=0,kb
        do j=0,jb

           vel_norm = sum(b%w(phy1,j,k,2:4)*fnvec(1:3))/b%w(phy1,j,k,1)

           if(vel_norm <= 0.0d0) then ! inlet
              rhoe = b%w(phy1,j,k,1)
              ue   = b%w(phy1,j,k,2)/rhoe
              pe   = b%p(phy1,j,k)
              ce   = sqrt(gamma*pe/rhoe)
              
              rpo  = ue*fnvec(1) +2.0d0*ce/(gamma -1.0d0)
              rmi  = u0*fnvec(1) -2.0d0*c0/(gamma -1.0d0)
           
              ui   = 0.5d0*(rpo +rmi)*fnvec(1)
              ci   = 0.25d0*(gamma -1.0d0)*(rpo -rmi)
              vi   = v0
              wi   = w0
              si   = rho0**gamma/p0
              cc   = ci**2.0d0/gamma
              rhoi = (si*cc)**gm
              prei   = rhoi*cc
              
              b%w(dum1,j,k,1) = rhoi
              b%w(dum1,j,k,2) = rhoi*ui
              b%w(dum1,j,k,3) = rhoi*vi
              b%w(dum1,j,k,4) = rhoi*wi
              b%w(dum1,j,k,5) = prei*gm +0.5d0*rhoi*(ui*ui +vi*vi +wi*wi)
              b%p(dum1,j,k)   = prei
           else if(vel_norm > 0.0d0) then ! outlet

              rhoi = b%w(phy1,j,k,1)
              ui   = b%w(phy1,j,k,2)/rhoi
              vi   = b%w(phy1,j,k,3)/rhoi
              wi   = b%w(phy1,j,k,4)/rhoi
              prei = b%p(phy1,j,k)
              si   = rhoi**gamma/prei
              ci   = sqrt(gamma*prei/rhoi)          ! flow variables interpolated from inner field
              
              rpo  = ui*fnvec(1) +2.0d0*ci/(gamma -1.0d0)
              rmi  = u0*fnvec(1) -2.0d0*c0/(gamma -1.0d0)                     ! Riemann invariables
              
              ue   = 0.5d0*(rpo +rmi)*fnvec(1)
              ce   = 0.25d0*(gamma -1.0d0)*(rpo -rmi)
              ve   = vi
              we   = wi
              se   = si
              cc   = ce**2.0d0/gamma
              rhoe = (se*cc)**gm                                     ! variables at the exit
              pe   = rhoe*cc 
              
              b%w(dum1,j,k,1) = rhoe
              b%w(dum1,j,k,2) = rhoe*ue
              b%w(dum1,j,k,3) = rhoe*ve
              b%w(dum1,j,k,4) = rhoe*we
              b%w(dum1,j,k,5) = pe*gm +0.5d0*rhoe*(ue*ue +ve*ve+ we*we)
              b%p(dum1,j,k)   = pe
           end if

           b%w(dum2,j,k,:) = b%w(dum1,j,k,:)
           b%p(dum2,j,k)   = b%p(dum1,j,k)

        end do
     end do
  end if

! on the front and back boundary, use 1-D Riemann invariables for 3D case
  if(face == 2 .or. face ==5) then
     ! front boundary
     if(face == 2) then
        phy1 = 2
        dum1 = 1
        dum2 = 0
     ! back boundary
     else if(face == 5) then
        phy1 = jl
        dum1 = je
        dum2 = jb
     end if
     fnvec(:) = b%facenvec(face,:)
     do k=0,kb
        do i=0,ib
           vel_norm = sum(b%w(i,phy1,k,2:4)*fnvec(1:3))/b%w(i,phy1,k,1)
  
!           if (b%w(i,phy1,k,3)/b%w(i,phy1,k,1)<0) then       ! use inlet boundary condition 
           if (vel_norm <= 0.0d0) then       ! use inlet boundary condition 
              
              rhoe = b%w(i,phy1,k,1)
              ve   = b%w(i,phy1,k,3)/rhoe
              pe   = b%p(i,phy1,k)
              ce   = sqrt(gamma*pe/rhoe)
              
              rpo  = ve*fnvec(2) +2.0d0*ce/(gamma -1.0d0)
              rmi  = v0*fnvec(2) -2.0d0*c0/(gamma -1.0d0)
              
              vi   = 0.5d0*(rpo +rmi)*fnvec(2)
              ci   = 0.25d0*(gamma -1.0d0)*(rpo -rmi)
              ui   = u0
              wi   = w0
              si   = rho0**gamma/p0
              cc   = ci**2.0d0/gamma
              rhoi = (si*cc)**gm
              prei   = rhoi*cc
              
              b%w(i,dum1,k,1) = rhoi
              b%w(i,dum1,k,2) = rhoi*ui
              b%w(i,dum1,k,3) = rhoi*vi
              b%w(i,dum1,k,4) = rhoi*wi
              b%w(i,dum1,k,5) = prei*gm +0.5d0*rhoi*(ui*ui +vi*vi +wi*wi)
              b%p(i,dum1,k)   = prei
!          else if(b%w(i,phy1,k,3)/b%w(i,phy1,k,1)>=0) then             ! use outlet boundary condition
           else if(vel_norm > 0.0d0) then             ! use outlet boundary condition
              
              rhoi = b%w(i,phy1,k,1)
              ui   = b%w(i,phy1,k,2)/rhoi
              vi   = b%w(i,phy1,k,3)/rhoi
              wi   = b%w(i,phy1,k,4)/rhoi
              prei = b%p(i,phy1,k)
              si   = rhoi**gamma/prei
              ci   = sqrt(gamma*prei/rhoi)
              
              rpo  = vi*fnvec(2) +2.0d0*ci/(gamma -1.0d0)
              rmi  = v0*fnvec(2) -2.0d0*c0/(gamma -1.0d0)
              
              ve   = 0.5d0*(rpo +rmi)*fnvec(2)
              ce   = 0.25d0*(gamma -1.0d0)*(rpo -rmi)
              ue   = ui
              we   = wi
              se   = si
              cc   = ce**2.0d0/gamma
              rhoe = (se*cc)**gm
              pe   = rhoe*cc

              b%w(i,dum1,k,1) = rhoe
              b%w(i,dum1,k,2) = rhoe*ue
              b%w(i,dum1,k,3) = rhoe*ve
              b%w(i,dum1,k,4) = rhoe*we
              b%w(i,dum1,k,5) = pe*gm +0.5d0*rhoe*(ue*ue +ve*ve +we*we)
              b%p(i,dum1,k)   = pe
              
              
           end if
           
           b%w(i,dum2,k,:) = b%w(i,dum1,k,:)
           b%p(i,dum2,k)   = b%p(i,dum1,k)
        end do
     end do
  end if

! on the top and bottom boundary
  if(face == 3 .or. face == 6) then
     ! bottom boundary
     if(face == 3) then
        phy1 = 2
        dum1 = 1
        dum2 = 0
     ! top boundary
     else if(face == 6) then
        phy1 = kl
        dum1 = ke
        dum2 = kb
     end if
     fnvec(:) = b%facenvec(face,:)
     do j=0,jb
        do i=0,ib
           vel_norm = sum(b%w(i,j,phy1,2:4)*fnvec(1:3))/b%w(i,j,phy1,1)

!           if (b%w(i,j,phy1,4)/b%w(i,j,phy1,1)>0 ) then       ! use inlet boundary condition         
           if (vel_norm <=0.0d0 ) then       ! use inlet boundary condition         
              rhoe = b%w(i,j,phy1,1)
              we   = b%w(i,j,phy1,4)/rhoe
              pe   = b%p(i,j,phy1)
              ce   = sqrt(gamma*pe/rhoe)
              
              rpo  = we*fnvec(3) +2.0d0*ce/(gamma -1.0d0)
              rmi  = w0*fnvec(3) -2.0d0*c0/(gamma -1.0d0)
              
              wi   = 0.5d0*(rpo +rmi)*fnvec(3)
              ci   = 0.25d0*(gamma -1.0d0)*(rpo -rmi)
              ui   = u0
              vi   = v0
              si   = rho0**gamma/p0
              cc   = ci**2.0d0/gamma
              rhoi = (si*cc)**gm
              prei = rhoi*cc
              
              b%w(i,j,dum1,1) = rhoi
              b%w(i,j,dum1,2) = rhoi*ui
              b%w(i,j,dum1,3) = rhoi*vi
              b%w(i,j,dum1,4) = rhoi*wi
              b%w(i,j,dum1,5) = prei*gm +0.5d0*rhoi*(ui*ui +vi*vi +wi*wi)
              b%p(i,j,dum1)   = prei
!          else if(b%w(i,j,phy1,4)/b%w(i,j,phy1,1)<=0 ) then             
!          use outlet boundary condition
           else if(vel_norm > 0.0d0) then
              rhoi = b%w(i,j,phy1,1)
              ui   = b%w(i,j,phy1,2)/rhoi
              vi   = b%w(i,j,phy1,3)/rhoi
              wi   = b%w(i,j,phy1,4)/rhoi
              prei = b%p(i,j,phy1)
              si   = rhoi**gamma/prei
              ci   = sqrt(gamma*prei/rhoi)
              
              rpo  = wi*fnvec(3) +2.0d0*ci/(gamma -1.0d0)
              rmi  = w0*fnvec(3) -2.0d0*c0/(gamma -1.0d0)

              we   = 0.5d0*(rpo +rmi)*fnvec(3)
              ce   = 0.25d0*(gamma -1.0d0)*(rpo -rmi)
              ue   = ui
              ve   = vi
              se   = si
              cc   = ce**2.0d0/gamma
              rhoe = (se*cc)**gm
              pe   = rhoe*cc
              
              b%w(i,j,dum1,1) = rhoe
              b%w(i,j,dum1,2) = rhoe*ue
              b%w(i,j,dum1,3) = rhoe*ve
              b%w(i,j,dum1,4) = rhoe*we
              b%w(i,j,dum1,5) = pe*gm +0.5d0*rhoe*(ue*ue +ve*ve +we*we)
              b%p(i,j,dum1)   = pe
             
           end if
           
           b%w(i,j,dum2,:) = b%w(i,j,dum1,:)
           b%p(i,j,dum2)   = b%p(i,j,dum1)
           
        end do
     end do
  end if
  return
end subroutine bc_far

