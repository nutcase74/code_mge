subroutine bc_periodic(b,face)
! subroutine treating periodic boundary conditions 
!
!
! ===================================================================
  use typelib4
  use typepara
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b
! local variables
  integer, intent(in) :: face
  integer :: i, j, k
  integer :: il,jl,kl,ie,je,ke,ib,jb,kb
  integer :: dum1,dum2,phy1,phy2
! *******************************************************************
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
     dum1 = 1
     dum2 = 0
     if(face==1) then
        phy1 = il
        phy2 = il-1
     else if(face==2) then
        phy1 = jl
        phy2 = jl-1
     else if(face==3) then
        phy1 = kl
        phy2 = kl-1
     end if
  else if(face==4 .or. face==5 .or. face==6) then
     phy1 = 2
     phy2 = 3
     if(face==4) then
        dum1 = ie
        dum2 = ib
     else if(face==5) then
        dum1 = je
        dum2 = jb
     else if(face==6) then
        dum1 = ke
        dum2 = kb
     end if
  end if

  !print *, face,dum1,dum2,phy1,phy2,jl,je,jb
!------------(translational) periodic boundary condition---------------------
! i-direction
  if(face==1 .or. face==4) then  
     do k=0,kb
        do j=0,jb
           b%w(dum2,j,k,:) = b%w(phy2,j,k,:)
           b%w(dum1,j,k,:) = b%w(phy1,j,k,:)
           b%p(dum2,j,k)   = b%p(phy2,j,k)
           b%p(dum1,j,k)   = b%p(phy1,j,k)
        end do
     end do

! j-direction
  else if(face==2 .or. face==5) then
     do k=0,kb
        do i=0,ib
           b%w(i,dum2,k,:) = b%w(i,phy2,k,:)
           b%w(i,dum1,k,:) = b%w(i,phy1,k,:)
           b%p(i,dum2,k)   = b%p(i,phy2,k)
           b%p(i,dum1,k)   = b%p(i,phy1,k)
        end do
     end do

! k-direction
  else if(face==3 .or. face==6) then
     do j=0,jb
        do i=0,ib
           b%w(i,j,dum2,:) = b%w(i,j,phy2,:)
           b%w(i,j,dum1,:) = b%w(i,j,phy1,:)
           b%p(i,j,dum2)   = b%p(i,j,phy2)
           b%p(i,j,dum1)   = b%p(i,j,phy1)
        end do
     end do
  else
     print *, "(bcperiodic) face number not specified : ", face
  end if

  return
end subroutine bc_periodic
