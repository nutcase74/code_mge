program main

  ! test the forward and backward interpolation on a 2x2x2 domain
  call test

end program main


subroutine test
  implicit none
  integer, parameter :: dp=kind(1.0d0), dof=3, nmarker=32, ncell=8

  real(dp) :: x(3,3,3,3), marker(nmarker,3), res(3,3,3), res_new(3,3,3), vs(8,3), vs_m(3)
  real(dp) :: dsize, xi(3), factor8(8), res_m(nmarker)
  real(dp) :: matxm(nmarker,8),matx(8,8),invx(8,8),force(8),res1(8)
  integer :: i,j,k,m,n,id, idx, nn(8,8),shft(3),shft1(3)
  integer :: cell(8,3), mcell(nmarker)
  integer :: ic,jc,kc
 
  nn = 1          ! number of overlapping cell center points

  ! cell coordinate and residual on each cell
  n = 0
  dsize = 1.0d0
  !print *, "residual on each node:"
  do i= 1,3
     do j=1,3
        do k=1,3
           n = n+1
           x(i,j,k,1) =(i-1)*dsize 
           x(i,j,k,2) =(j-1)*dsize 
           x(i,j,k,3) =(k-1)*dsize 
           !res(n) = dsize
           res(i,j,k) = (3-i)*dsize+1
           !print*, n, res(i,j,k)
        end do
     end do
  end do
  !res_new = res
  res_new = 0.0d0

  ! marker coordinate
!  marker(1,:) = (/1.7d0, 0.3d0, 0.7551d0 /)
!  marker(2,:) = (/1.5d0, 0.5d0, 0.6258d0 /)
!  marker(3,:) = (/1.8d0, 0.6d0, 0.3367d0 /)
!  marker(4,:) = (/1.3d0, 0.7d0, 0.7551d0 /)

!  marker(5,:) = (/1.6d0, 0.7d0, 0.3755d0 /)
!  marker(6,:) = (/1.3d0, 0.9d0, 0.6258d0 /)
!  marker(7,:) = (/1.6d0, 0.9d0, 0.3144d0 /)
!  marker(8,:) = (/1.9d0, 0.9d0, 0.2126d0 /)
  
  marker(1,:) = (/1.7d0, 0.1d0, 0.7551d0 /)
  marker(2,:) = (/1.7d0, 0.2d0, 0.7551d0 /)
  marker(3,:) = (/1.7d0, 0.3d0, 0.7551d0 /)
  marker(4,:) = (/1.7d0, 0.4d0, 0.7551d0 /)

  marker(5,:) = (/1.7d0, 0.5d0, 0.7551d0 /)
  marker(6,:) = (/1.7d0, 0.6d0, 0.7551d0 /)
  marker(7,:) = (/1.7d0, 0.7d0, 0.7551d0 /)
! marker(8,:) = (/1.7d0, 0.8d0, 0.7551d0 /)
  marker(8,:) = (/1.6d0, 1.3d0, 0.3755d0 /)

!  marker(9,:) = (/1.9d0, 0.9d0, 0.2126d0 /)
  marker(9 ,:) = (/1.7d0, 1.7d0, 0.7551d0 /)
  marker(10,:) = (/1.5d0, 1.5d0, 0.6258d0 /)
  marker(11,:) = (/1.8d0, 1.4d0, 0.3367d0 /)
  marker(12,:) = (/1.3d0, 1.3d0, 0.7551d0 /)

  marker(13,:) = (/1.6d0, 1.3d0, 0.3755d0 /)
  marker(14,:) = (/1.3d0, 1.1d0, 0.6258d0 /)
  marker(15,:) = (/1.6d0, 1.1d0, 0.3144d0 /)
  marker(16,:) = (/1.9d0, 1.1d0, 0.2126d0 /)

  marker(17,:) = (/1.7d0, 1.7d0, 1.2449d0 /)
  marker(18,:) = (/1.5d0, 1.5d0, 1.3742d0 /)
  marker(19,:) = (/1.8d0, 1.4d0, 1.6633d0 /)
  marker(20,:) = (/1.3d0, 1.3d0, 1.2449d0 /)

  marker(21,:) = (/1.6d0, 1.3d0, 1.6245d0 /)
  marker(22,:) = (/1.3d0, 1.1d0, 1.3742d0 /)
  marker(23,:) = (/1.6d0, 1.1d0, 1.6856d0 /)
  marker(24,:) = (/1.9d0, 1.1d0, 1.7874d0 /)

  marker(25,:) = (/1.7d0, 0.3d0, 1.2449d0 /)
  marker(26,:) = (/1.5d0, 0.5d0, 1.3742d0 /)
  marker(27,:) = (/1.8d0, 0.6d0, 1.6633d0 /)
  marker(28,:) = (/1.3d0, 0.7d0, 1.2449d0 /)

  marker(29,:) = (/1.6d0, 0.7d0, 1.6245d0 /)
  marker(30,:) = (/1.3d0, 0.9d0, 1.3742d0 /)
  marker(31,:) = (/1.6d0, 0.9d0, 1.6856d0 /)
  marker(32,:) = (/1.9d0, 0.9d0, 1.7874d0 /)


  ! cell index
!  cell(1,:) = (/ 1,2,4,5,10,11,13,14 /)
!  cell(2,:) = (/ 2,3,5,6,11,12,14,15 /)
!  cell(3,:) = (/ 4,5,7,8,13,14,16,17 /)
!  cell(4,:) = (/ 5,6,8,9,14,15,17,18 /)
!  cell(5,:) = (/ 10,11,13,14,19,20,22,23 /)
!  cell(6,:) = (/ 11,12,14,15,20,21,23,24 /)
!  cell(7,:) = (/ 13,14,16,17,22,23,25,26 /)
!  cell(8,:) = (/ 14,15,17,18,23,24,26,27 /)

  cell(1,:) = (/ 1,1,1 /)
  cell(2,:) = (/ 2,1,1 /)
  cell(3,:) = (/ 1,2,1 /)
  cell(4,:) = (/ 2,2,1 /)
  cell(5,:) = (/ 1,1,2 /)
  cell(6,:) = (/ 2,1,2 /)
  cell(7,:) = (/ 1,2,2 /)
  cell(8,:) = (/ 2,2,2 /)

  ! cell id which contains the marker
  mcell(1 :8 ) = 2
  mcell(9 :16) = 4
  mcell(17:24) = 8
  mcell(25:32) = 6


           
   ! count the number of markers associated with common node 
   do m = 2,8,2
     shft = cell(m,:)
     id = 0
     do k=shft(3),shft(3)+1
        do j=shft(2),shft(2)+1
           do i=shft(1),shft(1)+1
              id = id + 1
              do n = m+2,8,2
                 shft1 = cell(n,:)
                 idx = 0
                 do kc=shft1(3),shft1(3)+1
                    do jc=shft1(2),shft1(2)+1
                       do ic=shft1(1),shft1(1)+1
                          idx = idx+1
                          if(i==ic.and. j==jc.and. k==kc) then
                             nn(n,idx) = nn(n,idx)+1  
                             nn(m,id ) = nn(m,id )+1  
                          end if
                       end do
                    end do
                 end do
                 
              end do
            end do
        end do
     end do


   end do

  print *,
  !print *, "residual on each marker:"
  do i=1,nmarker

     id = mcell(i)
     n=0
     do m = cell(id,3),cell(id,3)+1
       do j = cell(id,2),cell(id,2)+1
         do k = cell(id,1),cell(id,1)+1
           n=n+1
           vs(n,:) = x(k,j,m,:)
           res1(n) = res(k,j,m)
           !print *,n,k,j,m,vs(n,:)
         end do
       end do
     end do

     !do j = 1,8
     !  idx     = cell(j,id)
     !  vs(j,:) = x(idx,:)
     !end do
     
     call trilinear_approx(vs, marker(i,:), xi)
     print *, i,xi,marker(i,:)

     call trilinear_factor(xi(1), xi(2), xi(3), factor8)
     print *, i,factor8
     
     res_m(i) = sum( factor8(1:8)*res1(1:8) )

     matxm(i,:) = factor8(:)
     !print *, i,matxm(i,:)
  end do

  print *, "force on 8 markers:"
  do i=2,ncell,2
     do j=1,8
        id = (i/2-1)*8+j
        matx(j,:) = matxm(id,:)
        res1(j)=res_m((i/2-1)*8+j)
     end do

     call svdinv(matx,8,8)
     !if (nmark(m)/=8) then
     !   invx = transpose(matx)
     !else
     !   invx = matx
     !end if
     invx = matx

     print *, i,res1
     !print *, invx
     force(:) = matmul(invx,res1)
     ! force interpolated backward
     print *, i, force

     idx=0
     do k = cell(i,3),cell(i,3)+1
        do j = cell(i,2),cell(i,2)+1
           do n = cell(i,1),cell(i,1)+1
              idx = idx+1
              res_new(n,j,k) = res_new(n,j,k)-force(idx)/real(nn(i,idx),dp)
              !print *, i,n,j,k,nn(i,idx),res_new(n,j,k)
           end do
        end do
     end do
     
  end do

  print *, "Final result:"
  do i = 1,3
     do j = 1,3
        do k = 1,3
           print *, i,j,k,res(i,j,k),res_new(i,j,k)
        end do
     end do
  end do

  return
end subroutine test
