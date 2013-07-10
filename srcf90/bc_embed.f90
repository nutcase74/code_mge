subroutine bc_embed(blk, lv)
!
!  CALC B.C. AT GHOST CELL FROM EMBEDDED NEXT COARSE GRID
!
!  DEVELOPMENT HISTORY: 
!
!  07/2011  : openMP added
! 
!***************************************************************
  use typelib4
  implicit none

  integer, parameter :: dp=kind(1.0d0), nval=5

  integer, intent(in) :: lv
  type(mgrid), pointer :: blk

  !  LOCAL VARIABLES

  type(mgrid), pointer   :: gridf, gridc
  type(tembed), dimension(:), pointer :: pdata
  integer, dimension(:), pointer :: popul

  integer  :: mylow(3), myupp(3)
  integer  :: i, j, k, p, levc
  integer  :: sr, nn, n, nerror, nwithin, ndone
  real(dp) :: pt(3), eta(3), xi(3), fac(8)
  integer  :: tarcel(3), err_cel, shft(3), err_shft, icel(3), low(3),upp(3)
  integer  :: lowdom(3), uppdom(3)
  real(dp) :: val(8,nval), qq, eqq

  include 'debug.inc'

  nullify(gridf, gridc, pdata, popul)
!
! GET BLOCK POPULATION
!
  call eval_popul(blk,popul)

!  if (bcond_overlap_off) then
!    print *, " NCOND OVERLAP IS OFF "
!    return
!  endif
!
! ----- GET COARSE GRID LEVEL FROM WHICH VALUES ARE INTERPOLATED TO CURRENT (FINE GRID) LEVEL
!
  if(lv .le. 1) return
  levc = lv - 1

  do sr = 1, popul(lv)

    gridf => get_mgrid(blk,lv,sr)
    if (.not. associated(gridf%emparent)) cycle

    pdata => gridf%emparent

    do nn = 1, gridf%iparent
        
      gridc => get_mgrid(blk, levc, pdata(nn)%serial )

     ! including two ghost cells if they fall on the physical domain of 
      ! the block with which it overlaps

      low   = pdata(nn)%low
      upp   = pdata(nn)%upp

      !! SET FINE GRID RANGE SUCH THAT IT INCLUDES DUMMY CELLS
      mylow = gridf%lowcell-2
      myupp = gridf%uppcell+2

      lowdom= gridc%lowcell
      uppdom= gridc%uppcell

      ! ---- KEEP TRACK OF ERROR
      nerror  = 0
      nwithin = 0
      ndone   = 0

      !print *,"bc_embed lv, levc  :",sr, lv, levc
!      print *,"bc_embed mylow /myupp :",mylow,myupp
      !print *,"lowdom/uppdom:",lowdom,uppdom
      !print *,"physical     :",gridf%lowcell,gridf%uppcell

      do i=mylow(1), myupp(1)
        do j=mylow(2), myupp(2)
          do k=mylow(3), myupp(3)

            icel = (/i,j,k/)  !! mycell index

            !! SKIP PHYSICAL CELLS OF FINE GRID
            if (all(icel>=gridf%lowcell.and.icel<=gridf%uppcell)) cycle

            call cell_center(gridf%ori, gridf%gcori, gridf%cm, icel, pt)
            
            !! FIND THE CELL ON PTR THAT CONTAINS pt
            call eval_x2cell(pt, gridc%ori,gridc%gcori, gridc%icm, lowdom, uppdom,&
                             tarcel, eta, err_cel)

            if (err_cel /= 0) cycle
            if (err_cel /= 0) cycle

            ndone = ndone +1
            if (.not. all((tarcel>=lowdom).and.(tarcel<=uppdom)) ) then
              nerror = nerror+1
              write(*,"(I3,3I4,a,3I4,a,3F10.4,a,3I4,a,3I4)") gridf%code,icel, &
                      " BC_EMBED OUTSIDE SEL DOM, tar",tarcel,&
                      " eta ",eta," LO ",low," HI ",upp
            else
              nwithin = nwithin + 1
            end if

            !! FOR CELL-CENTER SCHEME, GET THE VALID CELL-INDEX "SHFT", and LOCAL COORDINATES "XI"
            call shft_cell(shft, xi, lowdom, uppdom, tarcel, eta, err_shft)

            call trilinear_factor(xi(1), xi(2), xi(3), fac)
            !if(i==gridf%uppcell(1)+2.and. j==2 .and. k==4) then
            !   print *,"bc_em :",i,j,k,shft,xi,tarcel
            !end if   
            !if(i==gridf%uppcell(1)+2.and. j==2 .and. k==2) then
            !   print *,"bc_em :",i,j,k,shft,xi,tarcel
            !end if   
            
            val = reshape( gridc%w(shft(1):shft(1)+1,&
                                   shft(2):shft(2)+1,&
                                   shft(3):shft(3)+1,1:nval), (/8,nval/))

            if (gridf%blank(i,j,k) /= 0) then
              gridf%w(i,j,k,:) = (/ ( sum(val(:,p)*fac), p=1,nval ) /) 
            end if

            call depvars(gridf%w(i,j,k,:), gridf%p(i,j,k))

          end do
        end do
      end do


      if (nerror/=0) then
      write(*,"(a,3I5,3(a,I5))") " LEVEL=",lv,gridf%code,gridc%code," BC_EMBED, NDONE=",ndone, &
                           " NWITHIN=",nwithin, " NERROR=",nerror
      end if
    end do ! enddo nn loop
  end do ! enddo b loop
  
  return
end subroutine bc_embed
          


