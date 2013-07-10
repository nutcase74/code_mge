subroutine bc_overlap(blk,lv)

  use typelib4
  implicit none
  integer, parameter :: dp=kind(1.0d0), nval=5
  type(mgrid), pointer :: blk
  integer, intent(in) :: lv
!
! LOCAL VARIABLES
!
  type(mgrid), pointer   :: gtr, ptr
  type(tembed), dimension(:), pointer :: pdata
  integer, dimension(:), pointer :: popul
  integer  :: tarcel(3), err_cel, shft(3), err_shft, icel(3)
  integer  :: sublowc(3), subuppc(3), mylow(3), myupp(3)
  integer  :: i, j, k, p
  integer  :: b, nn, n, sr
  real(dp) :: cenpt(3), eta(3), xi(3), fac(8)

  real(dp) :: val(8,nval),valp(8)

  nullify(gtr, ptr, pdata, popul)
!
! GET BLOCK POPULATION
!
  call eval_popul(blk,popul)
!
! HANDLE ERROR
!
  if (lv>size(popul)) return

  do b=1, popul(lv)

    gtr   => get_mgrid(blk, lv, b )

    if (.not. associated(gtr%emoverlap)) cycle

    pdata => gtr%emoverlap

    do nn = 1, gtr%ioverlap   !! ioverlap is the number of block gtr overlaps with
        
      sr      =  pdata(nn)%serial
      ptr     => get_mgrid(blk, lv, sr)

      ! index of the overlapping part of the partner block (only physical cell)
      sublowc = pdata(nn)%low   
      subuppc = pdata(nn)%upp

      ! index of overlapping part of the current block      
      ! including two ghost cells if they fall on the physical domain of 
      ! the block with which it has overlapping

      mylow = pdata(nn)%mylow  
      myupp = pdata(nn)%myupp  
     
      do i=mylow(1), myupp(1)
        do j=mylow(2), myupp(2)
          do k=mylow(3), myupp(3)
            icel = (/i,j,k/)  !! mycell index
            
            if (all(icel>=gtr%lowcell.and.icel<=gtr%uppcell)) cycle

            call cell_center(gtr%ori, gtr%gcori, gtr%cm, icel, cenpt)
            
            !! FIND THE CELL ON PTR THAT CONTAINS pt
            call eval_x2cell(cenpt, ptr%ori,ptr%gcori, ptr%icm, ptr%lowcell, ptr%uppcell,&
                          tarcel, eta, err_cel)
            if(err_cel /= 0) cycle

            !! FOR CELL-CENTER SCHEME, GET THE VALID CELL-INDEX "SHFT", and LOCAL COORDINATES "XI"
            call shft_cell(shft, xi, ptr%lowcell, ptr%uppcell, tarcel, eta, err_shft)

            call trilinear_factor(xi(1), xi(2), xi(3), fac)
            
            val = reshape( ptr%w(shft(1):shft(1)+1,&
                                 shft(2):shft(2)+1,&
                                 shft(3):shft(3)+1,1:nval), (/8,nval/))
            valp = reshape(ptr%p(shft(1):shft(1)+1,&
                                 shft(2):shft(2)+1,&
                                 shft(3):shft(3)+1), (/8/))
            gtr%w(i,j,k,:) = (/ ( sum(val(:,p)*fac), p=1,nval ) /)
            !gtr%p(i,j,k)   =  sum(valp(:)*fac)

            call depvars(gtr%w(i,j,k,:), gtr%p(i,j,k))

          end do
        end do
      end do
    end do
  end do ! enddo b loop

  if (associated(popul)) deallocate(popul)
  return
end subroutine bc_overlap



