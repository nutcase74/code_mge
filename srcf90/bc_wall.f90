subroutine bc_wall(b,face)

! subroutine for treating no slip wall boundary
! u,v,w = 0 ( u(0:1)=-u(3:2), rho(0:1)=rho(3:2), E(0:1)=E(3:2) )
! dp/dn = 0 ( p(0:1)= p(2) )
!
!--------------------------------------------------------------
  use typelib4
  implicit none
  type(mgrid), intent(inout) :: b 
  integer, intent(in) :: face
  integer :: i, j, k
  integer :: il, jl, kl
  integer :: ie,je,ke,ib,jb,kb
  integer :: dum1,dum2,phy1,phy2
!--------------------------------------------------------------

  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb

  if(face==1 .or. face==2 .or. face==3) then
     phy1 = 2
     phy2 = 3
     dum1 = 1
     dum2 = 0
  else if(face==4) then
     phy1 = il
     phy2 = il-1
     dum1 = ie
     dum2 = ib
  else if(face==5) then
     phy1 = jl
     phy2 = jl-1
     dum1 = je
     dum2 = jb
  else if(face==6) then
     phy1 = kl
     phy2 = kl-1
     dum1 = ke
     dum2 = kb
  end if

  if(face == 1 .or. face ==4) then
    do k=0,kb
      do j=0,jb
        b%w(dum1,j,k,1)   =  b%w(phy1,j,k,1)
        b%w(dum1,j,k,5)   =  b%w(phy1,j,k,5)
        b%w(dum1,j,k,2:4) = -b%w(phy1,j,k,2:4)
        b%w(dum2,j,k,1)   =  b%w(phy2,j,k,1)
        b%w(dum2,j,k,5)   =  b%w(phy2,j,k,5)
        b%w(dum2,j,k,2:4) = -b%w(phy2,j,k,2:4)
        b%p(dum1,j,k)     =  b%p(phy1,j,k)
        b%p(dum2,j,k)     =  b%p(phy1,j,k)
      end do
    end do
  else if(face == 2 .or. face ==5) then
    do k=0,kb
      do i=0,ib
        b%w(i,dum1,k,1)   =  b%w(i,phy1,k,1)
        b%w(i,dum1,k,5)   =  b%w(i,phy1,k,5)
        b%w(i,dum1,k,2:4) = -b%w(i,phy1,k,2:4)
        b%w(i,dum2,k,1)   =  b%w(i,phy2,k,1)
        b%w(i,dum2,k,5)   =  b%w(i,phy2,k,5)
        b%w(i,dum2,k,2:4) = -b%w(i,phy2,k,2:4)
        b%p(i,dum1,k)     =  b%p(i,phy1,k)
        b%p(i,dum2,k)     =  b%p(i,phy1,k)
      end do
    end do
  else if(face == 3 .or. face ==6) then
    do j=0,jb
      do i=0,ib
        b%w(i,j,dum1,1)   =  b%w(i,j,phy1,1)
        b%w(i,j,dum1,5)   =  b%w(i,j,phy1,5)
        b%w(i,j,dum1,2:4) = -b%w(i,j,phy1,2:4)
        b%w(i,j,dum2,1)   =  b%w(i,j,phy2,1)
        b%w(i,j,dum2,5)   =  b%w(i,j,phy2,5)
        b%w(i,j,dum2,2:4) = -b%w(i,j,phy2,2:4)
        b%p(i,j,dum1)     =  b%p(i,j,phy1)
        b%p(i,j,dum2)     =  b%p(i,j,phy1)
      end do
    end do
  end if

  return
end subroutine bc_wall
