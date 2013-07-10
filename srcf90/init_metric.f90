subroutine init_metric(blk)
!
! Computation related to grid and body surface vector are conducted;
! Re-assign the number of parent and child mesh for information exchange
! between embedding grids; Find the overlapping region for 3d uniform grids;
! For 2D case, mesh correlation are computed cause missing in the preprocessor;
! For symmetric airfoil, values of symmetric intersection points are not the 
! same from the pre-process, so they are forced equivalent to get the symmetric
! results in steady calculation when angle of attack is 0 degree.
!
! ===================================================================================

  use iolib
  use typelib4
  implicit none
  integer, parameter   :: dp=kind(1.0d0)
  type(mgrid), pointer :: blk
!  integer, intent(in)  :: imcyc

  integer :: i, j, k, n, m
  integer :: ip, jp, kp
  integer :: ie,je,ke,ib,jb,kb

  real(dp), dimension(3) :: x1, x2 
  type(mgrid), pointer :: b

  integer :: icel(3), nerror, qerror
  real(dp) :: lref, gc(3), xtmp(3), del, mxerr, mxdif

  include 'debug.inc'
  include 'mgridtopo.h'
!  include 'plot3d.h'
  include 'cell_metric.h'

!!  print *, "Compute cell surface vectors for non-uniform cartesian grids"

  b=>blk
  do while (associated(b)) 

    ib = b%ib
    jb = b%jb
    kb = b%kb
    ie = b%ie
    je = b%je
    ke = b%ke
!
! ----- GRID POINTS COORDINATES
!
    lref   = maxval(b%cm)
    nerror = 0
    mxerr  = -huge(1.0d0)

    qerror = 0
    mxdif  = -huge(1.0d0)
    
    if (associated(b%x)) deallocate(b%x)
    allocate(b%x(0:b%ib,0:b%jb,0:b%kb,3))

    b%x = 0.0d0
    do k=1,b%kb
      do j=1,b%jb
        do i=1,b%ib
          icel = (/i,j,k/) 
          gc = real( (icel-2), dp )
          b%x(i,j,k,:) = matmul(b%cm,gc) + b%ori(:)
          call eval_gc2x(b%ori,b%gcori,b%cm, real(icel,dp), xtmp)
          del = sqrt(sum((b%x(i,j,k,:) - xtmp)**2))
          if (del>0.001d0*lref) qerror = qerror +1
          mxdif = max (del, mxdif)
        end do
      end do
    end do

    if (qerror/=0) then
      print *," WARNING:: GRID ERROR ", qerror," MXDIF ", real(mxdif)," LREF ", real(lref)
    end if
!
! ----- CELL-CENTER COORDINATES
!
    if (associated(b%xcen)) deallocate(b%xcen)
    allocate(b%xcen(0:b%ib,0:b%jb,0:b%kb,3))

    b%xcen = 0.0d0
    do k=1,b%kl+1
      do j=1,b%jl+1
        do i=1,b%il+1

          b%xcen(i,j,k,:) = .125d0*(b%x(i,  j,k,  :) +b%x(i,  j+1,k+1,:) +&
                                    b%x(i+1,j,k+1,:) +b%x(i+1,j+1,k,  :) +&
                                    b%x(i,  j,k+1,:) +b%x(i,  j+1,k,  :) +&
                                    b%x(i+1,j,k,  :) +b%x(i+1,j+1,k+1,:) )

        !!\BEGIN OPTION: CHECK ACCURACY OF XCEN
        !xtmp(:) = (/ ( (0.125d0*sum(b%x(i:i+1,j:j+1,k:k+1,d))),d=1,3) /)

          icel = (/ i,j,k /)
          call cell_center(b%ori, b%gcori, b%cm, icel, xtmp)
          del = sqrt(sum((b%xcen(i,j,k,:) - xtmp)**2))
          if (del>0.001d0*lref) nerror = nerror +1
          mxerr = max (del, mxerr)
        !!\END OPTION
        end do
      end do
    end do

!    print *," XCEN ERROR ", nerror," MXERR ", real(mxerr)," LREF ", real(lref)
!    print *," GRID ERROR ", qerror," MXDIF ", real(mxdif)," LREF ", real(lref)
!
! ---- FACE VECTORS
!
    allocate(b%si(0:ib,0:jb,0:kb,3))
    allocate(b%sj(0:ib,0:jb,0:kb,3))
    allocate(b%sk(0:ib,0:jb,0:kb,3))

    b%si   = 0.0d0
    b%sj   = 0.0d0
    b%sk   = 0.0d0

    do k=0,ke
      kp = k +1
      do j=0,je
        jp = j +1
        do i=0,ie
          ip = i +1
          x1(:) = b%x(i,jp,k,:) - b%x(i,j,kp,:)
          x2(:) = b%x(i,j ,k,:) - b%x(i,jp,kp,:)
          b%si(i,j,k,1) = 0.5d0*(x1(2)*x2(3) -x1(3)*x2(2))
          b%si(i,j,k,2) = 0.5d0*(x1(3)*x2(1) -x1(1)*x2(3))
          b%si(i,j,k,3) = 0.5d0*(x1(1)*x2(2) -x1(2)*x2(1))   

          if ((i==1) .and. (j==1) .and. (k==1)) then
            if (screen_verbose) print *,"METRIC :: BLK",b%level, b%serial, "  SI ",b%si(i,j,k,:)
          end if

        end do
      end do
    end do

    do k=0,ke
      kp = k +1
      do j=0,je
        jp = j +1
        do i=0,ie
          ip = i +1
          x1(:) = b%x(ip,j,k ,:) - b%x(i,j,kp,:)
          x2(:) = b%x(ip,j,kp,:) - b%x(i,j,k ,:)
          b%sj(i,j,k,1) = 0.5d0*(x1(2)*x2(3) -x1(3)*x2(2))
          b%sj(i,j,k,2) = 0.5d0*(x1(3)*x2(1) -x1(1)*x2(3))
          b%sj(i,j,k,3) = 0.5d0*(x1(1)*x2(2) -x1(2)*x2(1))

          if ((i==1) .and. (j==1) .and. (k==1)) then
            if (screen_verbose) print *,"METRIC :: BLK",b%level, b%serial," SJ ",b%sj(i,j,k,:)
          end if

        end do
      end do
    end do

    do k=0,ke
      kp = k +1
      do j=0,je
        jp = j +1
        do i=0,ie
          ip = i +1
          x1(:) = b%x(ip,jp,k,:) - b%x(i,j ,k,:)
          x2(:) = b%x(ip,j ,k,:) - b%x(i,jp,k,:)
          b%sk(i,j,k,1) = 0.5d0*(x1(2)*x2(3) -x1(3)*x2(2))
          b%sk(i,j,k,2) = 0.5d0*(x1(3)*x2(1) -x1(1)*x2(3))
          b%sk(i,j,k,3) = 0.5d0*(x1(1)*x2(2) -x1(2)*x2(1))

          if ((i==1) .and. (j==1) .and. (k==1)) then
            if (screen_verbose) print *,"METRIC :: BLK",b%level, b%serial," SK ",b%sk(i,j,k,:)
          end if

        end do
      end do
    end do
!
! ---- SET BLANK VALUE
!
    if (associated(b%blank)) deallocate(b%blank)
    allocate(b%blank(0:b%ib,0:b%jb,0:b%kb))
    b%blank(:,:,:) = 1

    do k = 2,b%kl
      do j = 2,b%jl
        do i = 2,b%il
          !! N.B. PREPROCESSOR PROVIDE CELL INDEX THAT'S CONSISTENT WITH FLOW SOLVER, 
          !!      NO NEED TO OFFSET CELL INDEX
          if (b%celltype(i,j,k)==csolid) b%blank(i,j,k) = 0
        end do
      end do
    end do

!  if(ibm_on) then
!  do m = 1,b%nmir
!    i = b%z(m)%ccid(1)
!    j = b%z(m)%ccid(2)
!    k = b%z(m)%ccid(3)
!    b%blank(i,j,k) = 1
!    b%celltype(i,j,k) = 88
!  end do
!  end if
!
! ---  MARK OUT CUT CELLS FOR GHOST CELLS. (HOW? WHY?)
!
    b%blank(0,:,:)    = b%blank(2,:,:)
    b%blank(1,:,:)    = b%blank(2,:,:)
    b%blank(b%ie,:,:) = b%blank(b%il,:,:)
    b%blank(b%ib,:,:) = b%blank(b%il,:,:)
    b%blank(:,0,:)    = b%blank(:,2,:)
    b%blank(:,1,:)    = b%blank(:,2,:)
    b%blank(:,b%je,:) = b%blank(:,b%jl,:)
    b%blank(:,b%jb,:) = b%blank(:,b%jl,:)
    b%blank(:,:,0)    = b%blank(:,:,2)
    b%blank(:,:,1)    = b%blank(:,:,2)
    b%blank(:,:,b%ke) = b%blank(:,:,b%kl)
    b%blank(:,:,b%kb) = b%blank(:,:,b%kl)

    print *,count(b%blank==0), " NUMBER OF BLANK==0 in BLOCK of LEVEL/SERIAL", b%level,b%serial

    b=>b%nxt
  end do
 
  return
end subroutine init_metric
