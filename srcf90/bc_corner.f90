subroutine bc_corner(b)
! set values on block corner cells
! ===================================================================
  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), intent(inout) :: b
! local variables
  integer  :: i,j,k,il,jl,kl,ie,je,ke,ib,jb,kb
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

  ! print *,"before :",b%w(1 ,1 ,1 ,:)
  ! 12 corner edges, each edge patch has four edges.
  do k = 2,kl
    b%w(1 ,1 ,k,:) = b%w(2 ,1 ,k,:) + b%w(1 ,2 ,k,:) - b%w(2 ,2 ,k,:)
    !b%w(0 ,1 ,k,:) = b%w(1 ,1 ,k,:) + b%w(0 ,2 ,k,:) - b%w(1 ,2 ,k,:)
    !b%w(1 ,0 ,k,:) = b%w(2 ,0 ,k,:) + b%w(1 ,1 ,k,:) - b%w(2 ,1 ,k,:)
    !b%w(0 ,0 ,k,:) = b%w(1 ,0 ,k,:) + b%w(0 ,1 ,k,:) - b%w(1 ,1 ,k,:)

    call depvars(b%w(1 ,1 ,k,:), b%p(1 ,1 ,k))
    !call depvars(b%w(0 ,1 ,k,:), b%p(0 ,1 ,k))
    !call depvars(b%w(1 ,0 ,k,:), b%p(1 ,0 ,k))
    !call depvars(b%w(0 ,0 ,k,:), b%p(0 ,0 ,k))
     
    b%w(ie,1 ,k,:) = b%w(il,1 ,k,:) + b%w(ie,2 ,k,:) - b%w(il,2 ,k,:)
    !b%w(ib,1 ,k,:) = b%w(ie,1 ,k,:) + b%w(ib,2 ,k,:) - b%w(ie,2 ,k,:)
    !b%w(ie,0 ,k,:) = b%w(il,0 ,k,:) + b%w(ie,1 ,k,:) - b%w(il,1 ,k,:)
    !b%w(ib,0 ,k,:) = b%w(ie,0 ,k,:) + b%w(ib,1 ,k,:) - b%w(ie,1 ,k,:)

    call depvars(b%w(ie,1 ,k,:), b%p(ie,1 ,k))
    !call depvars(b%w(ib,1 ,k,:), b%p(ib,1 ,k))
    !call depvars(b%w(ie,0 ,k,:), b%p(ie,0 ,k))
    !call depvars(b%w(ib,0 ,k,:), b%p(ib,0 ,k))

    b%w(1 ,je,k,:) = b%w(2 ,je,k,:) + b%w(1 ,jl,k,:) - b%w(2 ,jl,k,:)
    !b%w(1 ,jb,k,:) = b%w(2 ,jb,k,:) + b%w(1 ,je,k,:) - b%w(2 ,je,k,:)
    !b%w(0 ,je,k,:) = b%w(1 ,je,k,:) + b%w(0 ,jl,k,:) - b%w(1 ,jl,k,:)
    !b%w(0 ,jb,k,:) = b%w(1 ,jb,k,:) + b%w(0 ,je,k,:) - b%w(1 ,je,k,:)
    call depvars(b%w(1 ,je,k,:), b%p(1 ,je,k))
    !call depvars(b%w(1 ,jb,k,:), b%p(1 ,jb,k))
    !call depvars(b%w(0 ,je,k,:), b%p(0 ,je,k))
    !call depvars(b%w(0 ,jb,k,:), b%p(0 ,jb,k))

    b%w(ie,je,k,:) = b%w(il,je,k,:) + b%w(ie,jl,k,:) - b%w(il,jl,k,:)
    !b%w(ie,jb,k,:) = b%w(il,jb,k,:) + b%w(ie,je,k,:) - b%w(il,je,k,:)
    !b%w(ib,je,k,:) = b%w(ie,je,k,:) + b%w(ib,jl,k,:) - b%w(ie,jl,k,:)
    !b%w(ib,jb,k,:) = b%w(ie,jb,k,:) + b%w(ib,je,k,:) - b%w(ie,je,k,:)
    call depvars(b%w(ie,je,k,:), b%p(ie,je,k))
    !call depvars(b%w(ie,jb,k,:), b%p(ie,jb,k))
    !call depvars(b%w(ib,je,k,:), b%p(ib,je,k))
    !call depvars(b%w(ib,jb,k,:), b%p(ib,jb,k))
  end do

  do j = 2,jl
    b%w(1 ,j ,1,:) = b%w(2 ,j ,1,:) + b%w(1 ,j ,2,:) - b%w(2 ,j ,2,:)
    !b%w(0 ,j ,1,:) = b%w(1 ,j ,1,:) + b%w(0 ,j ,2,:) - b%w(1 ,j ,2,:)
    !b%w(1 ,j ,0,:) = b%w(2 ,j ,0,:) + b%w(1 ,j ,1,:) - b%w(2 ,j ,1,:)
    !b%w(0 ,j ,0,:) = b%w(1 ,j ,0,:) + b%w(0 ,j ,1,:) - b%w(1 ,j ,1,:)
    call depvars(b%w(1 ,j ,1,:), b%p(1 ,j ,1))
    !call depvars(b%w(0 ,j ,1,:), b%p(0 ,j ,1))
    !call depvars(b%w(1 ,j ,0,:), b%p(1 ,j ,0))
    !call depvars(b%w(0 ,j ,0,:), b%p(0 ,j ,0))
     
    b%w(ie,j ,1,:) = b%w(il,j ,1,:) + b%w(ie,j ,2,:) - b%w(il,j ,2,:) 
    !b%w(ib,j ,1,:) = b%w(ie,j ,1,:) + b%w(ib,j ,2,:) - b%w(ie,j ,2,:)
    !b%w(ie,j ,0,:) = b%w(il,j ,0,:) + b%w(ie,j ,1,:) - b%w(il,j ,1,:)
    !b%w(ib,j ,0,:) = b%w(ie,j ,0,:) + b%w(ib,j ,1,:) - b%w(ie,j ,1,:)
    call depvars(b%w(ie,j ,1,:), b%p(ie,j ,1))
    !call depvars(b%w(ib,j ,1,:), b%p(ib,j ,1))
    !call depvars(b%w(ie,j ,0,:), b%p(ie,j ,0))
    !call depvars(b%w(ib,j ,0,:), b%p(ib,j ,0))

    b%w(1 ,j,ke,:) = b%w(2 ,j,ke,:) + b%w(1 ,j,kl,:) - b%w(2 ,j,kl,:)
    !b%w(1 ,j,kb,:) = b%w(2 ,j,kb,:) + b%w(1 ,j,ke,:) - b%w(2 ,j,ke,:)
    !b%w(0 ,j,ke,:) = b%w(1 ,j,ke,:) + b%w(0 ,j,kl,:) - b%w(1 ,j,kl,:)
    !b%w(0 ,j,kb,:) = b%w(1 ,j,kb,:) + b%w(0 ,j,ke,:) - b%w(1 ,j,ke,:)
    call depvars(b%w(1 ,j ,ke,:), b%p(1 ,j ,ke))
    !call depvars(b%w(1 ,j ,kb,:), b%p(1 ,j ,kb))
    !call depvars(b%w(0 ,j ,ke,:), b%p(0 ,j ,ke))
    !call depvars(b%w(0 ,j ,kb,:), b%p(0 ,j ,kb))

    b%w(ie,j,ke,:) = b%w(il,j,ke,:) + b%w(ie,j,kl,:) - b%w(il,j,kl,:)
    !b%w(ie,j,kb,:) = b%w(il,j,kb,:) + b%w(ie,j,ke,:) - b%w(il,j,ke,:)
    !b%w(ib,j,ke,:) = b%w(ie,j,ke,:) + b%w(ib,j,kl,:) - b%w(ie,j,kl,:)
    !b%w(ib,j,kb,:) = b%w(ie,j,kb,:) + b%w(ib,j,ke,:) - b%w(ie,j,ke,:)
    call depvars(b%w(ie,j ,ke,:), b%p(ie,j ,ke))
    !call depvars(b%w(ie,j ,kb,:), b%p(ie,j ,kb))
    !call depvars(b%w(ib,j ,ke,:), b%p(ib,j ,ke))
    !call depvars(b%w(ib,j ,kb,:), b%p(ib,j ,kb))
  end do

  do i = 0,ib
    b%w(i,1 ,1 ,:) = b%w(i,1 ,2 ,:) + b%w(i,2 ,1 ,:) - b%w(i,2 ,2 ,:)
    !b%w(i,0 ,1 ,:) = b%w(i,1 ,1 ,:) + b%w(i,0 ,2 ,:) - b%w(i,1 ,2 ,:)
    !b%w(i,1 ,0 ,:) = b%w(i,2 ,0 ,:) + b%w(i,1 ,1 ,:) - b%w(i,2 ,1 ,:)
    !b%w(i,0 ,0 ,:) = b%w(i,1 ,0 ,:) + b%w(i,0 ,1 ,:) - b%w(i,1 ,1 ,:)
    call depvars(b%w(i ,1 ,1 ,:), b%p(i ,1 ,1))
    !call depvars(b%w(i ,0 ,1 ,:), b%p(i ,0 ,1))
    !call depvars(b%w(i ,1 ,0 ,:), b%p(i ,1 ,0))
    !call depvars(b%w(i ,0 ,0 ,:), b%p(i ,0 ,0))
     
    b%w(i,1 ,ke,:) = b%w(i,2 ,ke,:) + b%w(i,1 ,kl,:) - b%w(i,2 ,kl,:)
    !b%w(i,1 ,kb,:) = b%w(i,2 ,kb,:) + b%w(i,1 ,ke,:) - b%w(i,2 ,ke,:)
    !b%w(i,0 ,ke,:) = b%w(i,1 ,ke,:) + b%w(i,0 ,kl,:) - b%w(i,1 ,kl,:)
    !b%w(i,0 ,kb,:) = b%w(i,1 ,kb,:) + b%w(i,0 ,ke,:) - b%w(i,1 ,ke,:)
    call depvars(b%w(i ,1 ,ke,:), b%p(i ,1 ,ke))
    !call depvars(b%w(i ,1 ,kb,:), b%p(i ,1 ,kb))
    !call depvars(b%w(i ,0 ,ke,:), b%p(i ,0 ,ke))
    !call depvars(b%w(i ,0 ,kb,:), b%p(i ,0 ,kb))

    b%w(i,je,1 ,:) = b%w(i,je,1 ,:) + b%w(i,jl,2 ,:) - b%w(i,jl,2 ,:)
    !b%w(i,jb,1 ,:) = b%w(i,jb,1 ,:) + b%w(i,je,2 ,:) - b%w(i,je,2 ,:)
    !b%w(i,je,0 ,:) = b%w(i,je,0 ,:) + b%w(i,jl,0 ,:) - b%w(i,jl,1 ,:)
    !b%w(i,jb,0 ,:) = b%w(i,jb,0 ,:) + b%w(i,je,0 ,:) - b%w(i,je,1 ,:)
    call depvars(b%w(i ,je,1 ,:), b%p(i ,je,1))
    !call depvars(b%w(i ,jb,1 ,:), b%p(i ,jb,1))
    !call depvars(b%w(i ,je,0 ,:), b%p(i ,je,0))
    !call depvars(b%w(i ,jb,0 ,:), b%p(i ,jb,0))

    b%w(i,je,ke,:) = b%w(i,jl,ke,:) + b%w(i,je,kl,:) - b%w(i,jl,kl,:)
    !b%w(i,jb,ke,:) = b%w(i,je,ke,:) + b%w(i,jb,kl,:) - b%w(i,je,kl,:)
    !b%w(i,je,kb,:) = b%w(i,jl,kb,:) + b%w(i,je,ke,:) - b%w(i,jl,ke,:)
    !b%w(i,jb,kb,:) = b%w(i,je,kb,:) + b%w(i,jb,ke,:) - b%w(i,je,ke,:)
    call depvars(b%w(i ,je,ke,:), b%p(i ,je,ke))
    !call depvars(b%w(i ,jb,ke,:), b%p(i ,jb,ke))
    !call depvars(b%w(i ,je,kb,:), b%p(i ,je,kb))
    !call depvars(b%w(i ,jb,kb,:), b%p(i ,jb,kb))
  end do
  !print *,"after  :",b%w(1 ,1 ,1 ,:)

  ! 8 corner points
  b%w(1 ,1 ,1 ,:) = (b%w(2 ,1 ,1 ,:) + b%w(1 ,2 ,1 ,:) + b%w(1 ,1 ,2 ,:))*2.0d0/3.0d0 - b%w(2 ,2 ,2 ,:)
  b%w(1 ,1 ,ke,:) = (b%w(2 ,1 ,ke,:) + b%w(1 ,2 ,ke,:) + b%w(1 ,1 ,kl,:))*2.0d0/3.0d0 - b%w(2 ,2 ,kl,:)
  b%w(1 ,je,1 ,:) = (b%w(2 ,je,1 ,:) + b%w(1 ,jl,1 ,:) + b%w(1 ,je,2 ,:))*2.0d0/3.0d0 - b%w(2 ,jl,2 ,:)
  b%w(ie,1 ,1 ,:) = (b%w(il,1 ,1 ,:) + b%w(ie,2 ,1 ,:) + b%w(ie,1 ,2 ,:))*2.0d0/3.0d0 - b%w(il,2 ,2 ,:)
  b%w(ie,1 ,ke,:) = (b%w(il,1 ,ke,:) + b%w(ie,2 ,ke,:) + b%w(ie,1 ,kl,:))*2.0d0/3.0d0 - b%w(il,2 ,kl,:)
  b%w(ie,je,1 ,:) = (b%w(il,je,1 ,:) + b%w(ie,jl,1 ,:) + b%w(ie,je,2 ,:))*2.0d0/3.0d0 - b%w(il,jl,2 ,:)
  b%w(1 ,je,ke,:) = (b%w(2 ,je,ke,:) + b%w(1 ,jl,ke,:) + b%w(1 ,je,kl,:))*2.0d0/3.0d0 - b%w(2 ,jl,kl,:)
  b%w(ie,je,ke,:) = (b%w(il,je,ke,:) + b%w(ie,jl,ke,:) + b%w(ie,je,kl,:))*2.0d0/3.0d0 - b%w(il,jl,kl,:)

  call depvars( b%w(1 ,1 ,1 ,:), b%p(1 ,1 ,1 ))
  call depvars( b%w(1 ,1 ,ke,:), b%p(1 ,1 ,ke))
  call depvars( b%w(1 ,je,1 ,:), b%p(1 ,je,1 ))
  call depvars( b%w(ie,1 ,1 ,:), b%p(ie,1 ,1 ))
  call depvars( b%w(ie,1 ,ke,:), b%p(ie,1 ,ke))
  call depvars( b%w(ie,je,1 ,:), b%p(ie,je,1 ))
  call depvars( b%w(1 ,je,ke,:), b%p(1 ,je,ke))
  call depvars( b%w(ie,je,ke,:), b%p(ie,je,ke))

  return
end subroutine bc_corner
