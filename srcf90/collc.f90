subroutine collc(blk,mgl,sb)
!
! restrict flow variables and residuals
! Because the region of coarse grid is larger than that of finer grid,
! so the restriction are just implemented on part of coarse grid that 
! under the finer grid.
! ===============================================================
!
!  INTERPOLATION FROM FINE GRID TO COARSE GRID
!
!
! N.B. DOES IT INCLUDE GHOST CELLS 
!
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
! input variables
  type(mgrid), pointer :: blk
  integer, intent(in) :: mgl, sb
! local variables
  integer :: i, j, k, m, p
  type(mgrid), pointer :: b, ctr

  integer  :: ndone, nerror, nonmatch, found, levf, sr
  integer  :: icel(3), tarcel(3), shft(3), err_cel, err_shft
  integer  :: myupp(3), mylow(3)
  real(dp) :: eta(3), xi(3), approx(3), del(3), pt(3), octa(8,3), fac(8) !, myori(3)
  real(dp) :: val(8,nval) !, vsal(8,nval)
  integer  :: low(4), upp(4)

  b=>get_mgrid(blk,mgl,sb)

!  if (b%ichild==0) return

  nerror = 0
  ndone  = 0
  nonmatch = 0
  found    = 0

  low(:) = (/(lbound(b%w,i),i=1,4)/)
  upp(:) = (/(ubound(b%w,i),i=1,4)/)

  if (.not. associated(b%marker)) allocate( b%marker(low(1):upp(1),low(2):upp(2),low(3):upp(3)) )
  b%marker = 0

  !print *, "in COLLC, low /upp :", low,upp
  !print *, "in COLLC, lowc/uppc:", b%lowcell,b%uppcell

  ilop:&
  do m=1, b%ichild

    levf = b%emchild(m)%level
    sr   = b%emchild(m)%serial

    ctr => get_mgrid(blk, levf, sr)

    mylow = b%emchild(m)%mylow
    myupp = b%emchild(m)%myupp

    !print *, "collc mylow/myupp", mylow,myupp
    !print *, "collc bllow/blupp", b%lowcell, b%uppcell
    !print *, "collc cllow/clupp", ctr%lowcell, ctr%uppcell

    do i=mylow(1), myupp(1) 

      !! 0 layer GHOST CELL REQUIRED
      if ( i<b%lowcell(1).or. i>b%uppcell(1)) cycle

      do j=mylow(2), myupp(2)

        if ( j<b%lowcell(2).or. j>b%uppcell(2)) cycle

        do k=mylow(3), myupp(3) 

          if ( k<b%lowcell(3).or. k>b%uppcell(3)) cycle

          ndone = ndone + 1

          icel = (/i,j,k/)
          call cell_center(b%ori, b%gcori, b%cm, icel, pt)

          !! FIND THE CELL ON CTR THAT CONTAINS pt
          call eval_x2cell(pt, ctr%ori, ctr%gcori, ctr%icm, ctr%lowcell, &
                           ctr%uppcell, tarcel, eta, err_cel)

          if (err_cel/=0) nonmatch = nonmatch + 1
          if (err_cel/=0) cycle

          found = found +1

          !!************************************************************
          !! OPTION: COMPARE CELL-CENTER COORDIANTES
          call trilinear_factor(eta(1), eta(2), eta(3), fac)
          call vertex_octanode(ctr%ori, ctr%gcori, ctr%cm, tarcel, octa)
          approx = (/ (sum(octa(:,p)*fac),p=1,3) /)
          del = approx - pt
          if ( sqrt(sum(del*del))>epsilon(1.0)) then
            print *,"ERR:",sqrt(sum(del*del))
            nerror = nerror +1
          end if
          !!************************************************************

          !! FOR CELL-CENTER SCHEME, GET THE VALID CELL-INDEX "SHFT", and LOCAL COORDINATES "XI"
          call shft_cell(shft, xi, ctr%lowcell, ctr%uppcell, tarcel, eta, err_shft)
          call trilinear_factor(xi(1), xi(2), xi(3), fac)

          !print *,icel
          !if(j==2 .and. k==mylow(3)) print *,icel
          !if(i==26.and. j==2 .and. k==26) then
          !   print *,"collc :",i,j,k,shft,xi,tarcel
          !end if   
          !if(i==b%uppcell(1).and. j==2 .and. k==2) then
          !   print *,"addw :",i,j,k,shft,xi,tarcel
          !end if   

          b%wr(i,j,k,:) = 0.0d0
          !! COLLECT RESIDUAL FROM FINE GRID TO COARSE GRID
          val = reshape(ctr%dw(shft(1):shft(1)+1,&
                               shft(2):shft(2)+1,&
                               shft(3):shft(3)+1,1:nval), (/8,nval/))
          val = val/ctr%cellvolume
          b%wr(i,j,k,:) = (/ ( sum(val(:,p)*fac), p=1,nval ) /) *b%cellvolume

          !! INTERPOLATE CONS. VARIABLES FROM FINE GRID TO COARSE GRID
          val = reshape(ctr%w( shft(1):shft(1)+1,&
                               shft(2):shft(2)+1,&
                               shft(3):shft(3)+1,1:nval), (/8,nval/))

          b%w(i,j,k,:) = (/ ( sum(val(:,p)*fac), p=1,nval ) /)
!
! ----    SET MARKER TO FLAG CELL FOR HAVING %wr VALUE, USED TO DISTINGUISH CELLS 
!         IN GENERATING FORCE TERM
!
          b%marker(i,j,k) = 1
        end do
      end do
    end do

!    i=2
!    j=2
!    p=26
!    do k=2,myupp(3),2
!       print *, i,j,k,0.125d0*(sum(ctr%w(i:i+1,j:j+1,k:k+1,1)))
!       print *, 26,2,p,b%w(26,2,p,1)
!       p = p+1
!    end do


  end do ilop

!    print *," **********  COLLC_EM  ********** "
!    print *," B%id     ", b%id
!    print *," NDONE    ", ndone
!    print *," FOUND    ", found
!    print *," NERROR   ", nerror
!    print *," NONMATCH ", nonmatch

  return
end subroutine collc

