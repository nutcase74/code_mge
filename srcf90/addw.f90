subroutine addw(blk,mgl,mlvl,sb)

! prolongation from coarse grid to finer grid
! use volume weighting, transfer the variations of flow
! variables in a multigrid cycle
! =================================================================
!
! INTERPOLATION FROM COARSE GRID TO FINE GRID USING TRILINEAR INTERPOLATION
!
!!
!! N.B. DOES IT INCLUDE GHOST CELLS
!!
!*****************************************************************************
!
! SUBROUTINES CALLED : eval_gc2x  cell_center  eval_x2cell  trilinear_factor
!                      vertex_octanode  shft_cell
!
! FUNCTIONS CALLED   :
!
!*****************************************************************************

  use typelib4
  implicit none

  integer, parameter :: dp=kind(1.0d0), nval=5

! input/output variables
  type(mgrid), pointer :: blk
  integer, intent(in)  :: mlvl, mgl, sb

! local variables
  integer :: il,jl,kl,ie,je,ke,ib,jb,kb
  type(mgrid), pointer :: b, ptr
  real(dp), dimension(:,:,:,:), pointer :: store_w, dw
  integer :: i, j, k, p

  integer  :: m, ndone, nerror, nonmatch, found, levc, sr
  integer  :: sublowc(3), subuppc(3), icel(3), tarcel(3), shft(3), err_cel, err_shft
  integer  :: myupp(3), mylow(3)
  real(dp) :: subori(3), eta(3), xi(3), approx(3), del(3), pt(3), octa(8,3), fac(8), myori(3)
  real(dp) :: val(8,nval), vsal(8,nval)

  include 'trilinear.h'
  include 'cell_metric.h'
  include 'debug.inc'

  nullify(store_w, dw)

  b => get_mgrid(blk,mgl,sb)

  if (b%iparent==0) return

  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb

  allocate(store_w(0:ib,0:jb,0:kb,nval), dw(0:ib,0:jb,0:kb,nval)  )

  store_w = 0.0d0
  dw = 0.0d0

  do k=0,kb
    do j=0,jb
      do i=0,ib
        store_w(i,j,k,:) = b%w(i,j,k,:)
      end do
    end do
  end do

  nerror = 0
  ndone  = 0
  nonmatch = 0
  found = 0

 ! print *,"in ADDW :", mgl-1, mlvl

  ilop:&
  do m=1, size(b%emparent)

    levc = b%emparent(m)%level
    sr   = b%emparent(m)%serial

    ptr => get_mgrid(blk, levc, sr)

    sublowc = b%emparent(m)%low
    subuppc = b%emparent(m)%upp
    call eval_gc2x(ptr%ori,ptr%gcori,ptr%cm,real(sublowc,dp),subori)

    mylow = b%emparent(m)%mylow
    myupp = b%emparent(m)%myupp

    call eval_gc2x(b%ori,b%gcori,b%cm,real(mylow,dp),myori)

!    print *," ADDW PARENT GRID ",ptr%code, " DOMAIN ORI@ ", real(subori), "MYORI@",real(myori)

!    print *, "movfin from level: ", mgl-1," to level: ",mgl, "indices :", &
!         b%lowcell(1)-1,b%lowcell(2)-1,b%lowcell(3)-1,b%uppcell(1)+1,b%uppcell(2)+1,b%uppcell(3)+1

!    print *, "addw    mylow /myupp :",mylow,myupp, mgl-1,mlvl
    do i=mylow(1), myupp(1) 

      !! 0 layer GHOST CELL REQUIRED
      if ( i<b%lowcell(1) .or. i>b%uppcell(1))  cycle

      do j=mylow(2), myupp(2)

        if ( j<b%lowcell(2) .or. j>b%uppcell(2))  cycle

        do k=mylow(3), myupp(3)

          if ( k<b%lowcell(3) .or. k>b%uppcell(3))   cycle
      

          ndone = ndone + 1

          icel = (/i,j,k/)
          call cell_center(b%ori, b%gcori, b%cm, icel, pt)

          !! FIND THE CELL ON PTR THAT CONTAINS pt
          call eval_x2cell(pt, ptr%ori,ptr%gcori, ptr%icm, ptr%lowcell, ptr%uppcell, tarcel, eta, err_cel)

          if (err_cel/=0) nonmatch = nonmatch + 1
          if (err_cel/=0) cycle

          found = found +1

          !!************************************************************
          !! OPTION: COMPARE CELL-CENTER COORDIANTES
          call trilinear_factor(eta(1), eta(2), eta(3), fac)
          call vertex_octanode(ptr%ori, ptr%gcori, ptr%cm, tarcel, octa)
          approx = (/ (sum(octa(:,p)*fac),p=1,3) /)
          del = approx - pt
          if ( sqrt(sum(del*del))>epsilon(1.0)) then
            print *,"ERR:",sqrt(sum(del*del))
            nerror = nerror +1
          end if
          !!************************************************************

          !! FOR CELL-CENTER SCHEME, GET THE VALID CELL-INDEX "SHFT", and LOCAL COORDINATES "XI"
          call shft_cell(shft, xi, ptr%lowcell, ptr%uppcell, tarcel, eta, err_shft)

          call trilinear_factor(xi(1), xi(2), xi(3), fac)
       !   if(i==2.and. j==2 .and. k==2) then
       !      print *,"addw :",i,j,k,shft,xi,tarcel
       !   end if   
       !   if(i==2.and. j==2 .and. k==b%uppcell(3)) then
       !      print *,"addw :",i,j,k,shft,xi,tarcel
       !   end if   
          
          val = reshape( ptr%w(shft(1):shft(1)+1,&
                               shft(2):shft(2)+1,&
                               shft(3):shft(3)+1,1:nval), (/8,nval/))

          if (mgl-1 /= mlvl) then
            vsal= reshape(ptr%ws(shft(1):shft(1)+1,&
                                 shft(2):shft(2)+1,&
                                 shft(3):shft(3)+1,1:nval), (/8,nval/))
            val = val -vsal
          end if

          dw(i,j,k,:) = (/ ( sum(val(:,p)*fac), p=1,nval ) /)
         
          if(mgl-1 == mlvl) then
            if(b%blank(i,j,k)/=0) b%w(i,j,k,:) = dw(i,j,k,:)
          elseif(mgl-1 /= mlvl) then
            if(b%blank(i,j,k)/=0) b%w(i,j,k,:) = store_w(i,j,k,:) +dw(i,j,k,:)
          end if

        end do
      end do
    end do

!    i=2
!    j=2
!    p=26
!    do k=2,myupp(3),2
!       print *, i,j,k,0.125d0*(sum(b%w(i:i+1,j:j+1,k:k+1,1)))
!       print *, 26,2,p,ptr%w(26,2,p,1)
!       p = p+1
!    end do

  end do ilop

!  if (debug) then
!    print *," **********  ADDW_EM  ********** "
!    print *," B%id     ", b%id
!    print *," NDONE    ", ndone
!    print *," FOUND    ", found
!    print *," NERROR   ", nerror
!    print *," NONMATCH ", nonmatch
!  end if

  deallocate(store_w,dw)
  return
end subroutine addw
