subroutine euler(b,itime,ncyc,mgl,mlvl,lres,lfir,lwork)
!**********************************************************************************************
!
! USE RUNGE-KUTTA SIMPLIFIED 5-STAGE TIME STEPPING SCHEME TO SOLVE THE
! ORDINARY DIFFERENTIAL EQUATIONS DISCRETED BY THE FINITE VOLUME METHOD.
!
! ACCELERATION TECHNIQUES INCLUDING LOCAL TIME STEP, RESIDUAL SMOOTH AND
! MULTIGRID ARE USED. 
!
! IN MULTIGRID PART, FORCING TERM ARE ADDED TO THE RESIDUAL ON THE RIGHT HAND SIDE. 
!
! FOR UNSTEADY CALCULATION, TIME DERIVATIVES ARE INCORPORATED INTO THE
! RESIDUAL ON THE RIGHT HAND SIDE TO USE THE SAME TIME STEPPING SCHEME FOR
! STEADY/UNSTEADY CALCULATION.
!
!**********************************************************************************************
!
! DEVELOPEMNT HISTORY:
!
! 2011-08-31  added lwork for MPI operations
!
!**********************************************************************************************

  use typelib
  implicit none

  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), pointer :: b
  integer, intent(in) :: itime 
  integer, intent(in) :: ncyc,mlvl,mgl
  logical, intent(in) :: lres, lfir, lwork
!
! LOCAL variables
!
  integer  :: i, j, k, n, m, id
  real(dp) :: fn, de, rt, lam, rhoq
  real(dp) :: fcoll, uu, vv, ww, ee
  integer  :: nx,ny,nz,il,jl,kl,ib,jb,kb,ie,je,ke

  real(dp), dimension(:,:,:), pointer :: tem
  real(dp), dimension(:,:,:,:), pointer :: vel
  real(dp), dimension(:,:,:,:), pointer :: dui,duj,duk
!  real(dp), dimension(:,:,:,:), pointer :: gradfi,gradfj,gradfk
!  real(dp), dimension(:,:,:,:), pointer :: gradti,gradtj,gradtk
! ******************************************************************

  real*8 :: tmp_finish, tmp_start

  include 'debug.inc'
  include 'timer.inc'
!  include 'mgridtopo.h'
!
! ----- RETURN IF BLOCK DOES NOT BELONG TO CURRENT PID 
!
  if (.not. lwork) return

  nx = b%ncell(1)
  ny = b%ncell(2)
  nz = b%ncell(3)
  il = b%il
  jl = b%jl
  kl = b%kl
  ie = b%ie
  je = b%je
  ke = b%ke
  ib = b%ib
  jb = b%jb
  kb = b%kb

  fcoll = 1.0d0

! save the initial values on first entry to a coarser
! mesh during the multigrid cycle

  if (lfir) then
    b%ws = b%w
  end if

! save the flow variables for the R-K scheme
  do k=2,kl
    do j=2,jl
      do i=2,il
        b%wn(i,j,k,:) = b%w(i,j,k,:)
      end do
    end do
  end do

  b%dw = 0.0d0
  b%fw = 0.0d0

  if(distype == 'r2' .or. distype == 'R2') then
    allocate(dui(0:ib,0:jb,0:kb,5),duj(0:ib,0:jb,0:kb,5),duk(0:ib,0:jb,0:kb,5))
    dui = 0.0d0
    duj = 0.0d0
    duk = 0.0d0
  end if

  if(eqntype == 'n' .or. eqntype == 'N') then
!    allocate(gradfi(0:ib,0:jb,0:kb,1:9),gradfj(0:ib,0:jb,0:kb,1:9), &
!             gradfk(0:ib,0:jb,0:kb,1:9),gradti(0:ib,0:jb,0:kb,1:3), &
!             gradtj(0:ib,0:jb,0:kb,1:3),gradtk(0:ib,0:jb,0:kb,1:3) )
    allocate(   vel(0:ib,0:jb,0:kb,1:3),   tem(0:ib,0:jb,0:kb) )
 
!    gradfi = 0.0d0
!    gradfj = 0.0d0
!    gradfk = 0.0d0
!    gradti = 0.0d0
!    gradtj = 0.0d0
!    gradtk = 0.0d0
    vel    = 0.0d0
    tem    = 0.0d0
  end if

  m = 1

! ---  MULTISTAGE SCHEME

  multistage: do while (m<=mstage)

    !print *, "stage :", m, "mgl :", mgl
    fn = cstp(m)*cfl

    !
    ! blending of the new and old artificial dissipation  
    ! initialize dissipation
    !  
    if(m > 1) then
      do k=2,kl
        do j=2,jl
          do i=2,il
            b%fw(i,j,k,:) = b%fw(i,j,k,:)*(1.0d0-cdis(m))
          end do
        end do
      end do
    end if

    if(eqntype == 'n' .or. eqntype == 'N') then

      ! compute velocity and temperature from conservative variables
      do i=0,ib
        do j=0,jb
          do k=0,kb
            vel(i,j,k,1:3)=b%w(i,j,k,2:4)/b%w(i,j,k,1)
            !rhoq = sum(b%w(i,j,k,2:4)**2.0d0)
            !t(i,j,k) = (gamma-1.0d0)*(b%w(i,j,k,5)-0.5d0*rhoq)/b%w(i,j,k,1)/gas_tmean
            tem(i,j,k) = b%p(i,j,k)/(b%w(i,j,k,1)*gas_const)

            if(b%blank(i,j,k) == 0) then
              vel(i,j,k,:) = 0.0d0
              tem(i,j,k)   = 0.0d0
            end if
          end do
        end do
      end do


      ! --- CALCULATE VISCID FLUX
      ! it is added to dissipation and applied with blending. eq(6.6)
      if (cdis(m) > 0.0d0) then

        !call wall_time(tmp_start)
        call fgradi(b,vel,tem)
        call fgradj(b,vel,tem)
        call fgradk(b,vel,tem)

!        call fgradi(b,gradfi,gradti,vel,tem)
!        call fgradj(b,gradfj,gradtj,vel,tem)
!        call fgradk(b,gradfk,gradtk,vel,tem)

        !call wall_time(tmp_finish)
        !op_gradface = op_gradface + (tmp_finish - tmp_start)
        !print *, "OP_GRADFACE : ", tmp_finish - tmp_start

        !call wall_time(tmp_start)
        call vflux (b,m,vel,tem)
!        call vflux (b,m,gradfi,gradfj,gradfk,gradti,gradtj,gradtk,vel,tem)
        !call wall_time(tmp_finish)
        !op_fluxvis = op_fluxvis + (tmp_finish - tmp_start)
        !print *, "OP_FLUXVIS  :", tmp_finish - tmp_start

      end if
    end if

    if(distype == 'c' .or. distype == 'C') then
      ! --- COMPUTE THE ARTIFICIAL DISSIPATION
       
      if (cdis(m) > 0.0d0) then
         !call wall_time(tmp_start)
         call dflux(b,m)
         !call wall_time(tmp_finish)
         !print *, "OP_FLUXdis  :", tmp_finish - tmp_start
      end if
      ! --- CALCULATE INVISCID FLUX
      !call wall_time(tmp_start)
      call eflux(b)
      !call wall_time(tmp_finish)
      !print *, "OP_FLUXCEN  :", tmp_finish - tmp_start
      !print *

    ! r2 (2nd order roe's scheme) and r1 (1st order) are not validated yet
    else if(distype == 'r2' .or. distype == 'R2') then
      !! ------- differences of primitive variables
      call var_diffs (b,dui,duj,duk)
      if (cdis(m) > 0.0d0) call dflux_roe2(b,m,dui,duj,duk)
      call eflux_roe2(b,dui,duj,duk)

    ! r1 (1st order roe's scheme))
    else if(distype == 'r1' .or. distype == 'R1') then
      if (cdis(m) > 0.0d0) call dflux_roe1(b,m)
      call eflux_roe1(b)

    else
      print *, "ERROR: distype should be c/r1/r2, here distype is :", distype
      stop
    end if

    do k=2,kl
      do j=2,jl
        do i=2,il
          b%dw(i,j,k,:) = b%dw(i,j,k,:) +b%fw(i,j,k,:)
        end do
      end do
    end do

   !if(ibm_on) call ibforce(b)
   !if(ibm_on) call force_int(b)
    
   !debug
   !do k = 2,kl
   !  do j = 2,jl
   !    do i = 2,il
   !      if (b%w(i,j,k,2) .ne. b%w(i,j,k,2) ) then
   !      !if (abs(b%w(i,j,k,1)) .lt. 10e-4 ) then
   !        print *,"euler1 ijk",ncyc,mgl,mlvl,b%x(i,j,k,:)
   !        print *,i,j,k,"dw",b%dw(i,j,k,1),"fw",b%fw(i,j,k,1)
   !        print *,i,j,k," w",b%w(i,j,k,:)
   !        stop
   !      end if
   !    end do
   !  end do
   !end do



! for unsteady calculation, incorporate transient terms to the residual on 
! the right hand side. Time integration use second-order backward difference,
! for the first physical time step, first-order difference is used.

    if (execmode==1) then
      !
      ! ---  TIME ACCURATE SOLUTION FOR UNSTEADY FLOW
      !

      if (itime>1) then
        do k=2,kl
          do j=2,jl
            do i=2,il
              b%dw(i,j,k,:) = b%dw(i,j,k,:) + 0.5d0*b%cellvolume/dtstep&
                   *(3.0d0*b%w(i,j,k,:)-4.0d0*b%w1(i,j,k,:)+b%w2(i,j,k,:))
            end do
          end do
        end do
      else if(itime==1) then

        do k=2,kl
          do j=2,jl
            do i=2,il
              b%dw(i,j,k,:) = b%dw(i,j,k,:) + 0.5d0*b%cellvolume/dtstep &
                                              *(b%w(i,j,k,:) -b%w1(i,j,k,:))
            end do
          end do
        end do
      end if

    end if

   !debug
   !do k = 2,kl
   !  do j = 2,jl
   !    do i = 2,il
   !      if (b%dw(i,j,k,1) .ne. b%dw(i,j,k,1) ) then
   !        print *,"euler2 ijk",i,j,k,"dw",b%dw(i,j,k,1),"fw",b%fw(i,j,k,1)
   !        stop
   !      end if
   !    end do
   !  end do
   !end do

! ********************************************************
! *                                                      *
! *   modification of the residues on any coarser mesh   *
! *                                                      *
! ********************************************************
    if (mgl .lt. mlvl) then            ! mgl>0 means in the coarser grid
     
      if (lfir.and.m.eq.1) then    ! for the coarse grid enter RK first time
        do k=2,kl
          do j=2,jl
            do i=2,il
              b%qf(i,j,k,:) = 0.0d0
              if (b%marker(i,j,k)==1) then
                b%qf(i,j,k,:) = fcoll*b%wr(i,j,k,:) -b%dw(i,j,k,:)
              end if
              b%dw(i,j,k,:) = b%dw(i,j,k,:) +b%qf(i,j,k,:)
            end do
          end do
        end do
      else 

        do k=2,kl
          do j=2,jl
            do i=2,il
              b%dw(i,j,k,:) = b%dw(i,j,k,:) +b%qf(i,j,k,:)
            end do
          end do
        end do

      end if

      if (execmode ==1 ) then 
        do k=2,kl
          do j=2,jl
            do i=2,il
              de = fn*b%vdtim(i,j,k)
              if (vt == 1) de = fn*b%dtim
              lam = 1.0d0 + de*3.0d0/2.0d0/dtstep
              b%dw(i,j,k,:) = b%dw(i,j,k,:)/lam
            end do
          end do
        end do 
      end if
    end if


! if it's solid cell, let the variables unchange
    do k=2,kl
      do j=2,jl
        do i=2,il
          b%dw(i,j,k,:) = b%dw(i,j,k,:)*real(b%blank(i,j,k),dp)
        end do
      end do
    end do

   if(ibm_on) call ibforce(b)

! for multigrid algorithm, just calculate the residual, do not update the flow variables
    if (lres) then
      if(associated(dui  ))  deallocate(dui)
      if(associated(duj  ))  deallocate(duj)
      if(associated(duk  ))  deallocate(duk)
!      if(associated(gradfi)) deallocate(gradfi)
!      if(associated(gradfj)) deallocate(gradfj)
!      if(associated(gradfk)) deallocate(gradfk)
!      if(associated(gradti)) deallocate(gradti)
!      if(associated(gradtj)) deallocate(gradtj)
!      if(associated(gradtk)) deallocate(gradtk)
      if(associated(vel))    deallocate(vel)
      if(associated(tem))    deallocate(tem)

      return
    end if

! multiple the time step
    do k=2,kl
      do j=2,jl
        do i=2,il
          de = fn*b%vdtim(i,j,k)/b%cellvolume
          if(vt==1) de = fn*b%dtim/b%cellvolume
          b%dw(i,j,k,:) = de*b%dw(i,j,k,:)
        end do
      end do
    end do

! Residual smoothing, only on odd stages
    if(mod(m,2)/=0) then
      call psmoo_old(b)
    end if

! update flow variables

    do k=2,kl
      do j=2,jl
        do i=2,il

          lam = 1.0d0

          if (execmode==1 .and. mgl==mlvl) then                  

            de = fn*b%vdtim(i,j,k)
            if (vt == 1) de = fn*b%dtim
            lam = 1.0d0+de*3.0d0/2.0d0/dtstep
          end if
          b%w(i,j,k,:) = b%wn(i,j,k,:) -b%dw(i,j,k,:)/lam
          call depvars(b%w(i,j,k,:), b%p(i,j,k))
          !if(b%x(i,j,k,1) > 0.51 .and. b%x(i,j,k,1)<0.53 .and. j==jl/2 .and. k==kl/2) &
          !     print *, i,j,k, b%w(i,j,k,:),b%p(i,j,k)
        end do
      end do
    end do

! compute for the next stage
    m = m + 1

!! 
!! update boundary condition, only bc_cutcell & bcond
!!
!!!---------------------------------------------------
   if(ibm_on) then
      !call bc_IBM1(b) 
      !call bc_pressure(b) 
      !call bc_velocity(b)
   else
      call bc_cutcell(b)
   end if

!!**********************************************************
!! removed EMBED & OVERLAPP IN EULER AS THEY ARE OBSOLETE
!!    call bc_embed(b,mgl)   ! 
!!    call bc_overlap(b,mgl)
!************************************************************
    call bcond(b)

  end do multistage
!
! to calculate the average residual
! The measure of convergence is the residual for the density
! defined as the root mean square value of delta_rho/delta_t
!
  if (debug) print *,"LEVEL SERIAL ", b%level, b%serial, " SIZE ", size(b%rtrms), "NCYC ", ncyc

  b%rtrms(ncyc) = 0.0d0

  do k=2,b%kl
    do j=2,b%jl
      do i=2,b%il
        de = cfl*b%vdtim(i,j,k)
        if (vt==1) de = cfl*b%dtim
        rt = 0.0d0
        if(b%blank(i,j,k) /= 0) then
           rt = b%dw(i,j,k,1)/de
        end if
        b%rtrms(ncyc) = b%rtrms(ncyc) + rt**2.0d0
      end do
    end do
  end do

  b%rtrms(ncyc) = sqrt(b%rtrms(ncyc)/real(nx*ny*nz,dp))

  if(associated(dui  ))  deallocate(dui)
  if(associated(duj  ))  deallocate(duj)
  if(associated(duk  ))  deallocate(duk)
!  if(associated(gradfi)) deallocate(gradfi)
!  if(associated(gradfj)) deallocate(gradfj)
!  if(associated(gradfk)) deallocate(gradfk)
!  if(associated(gradti)) deallocate(gradti)
!  if(associated(gradtj)) deallocate(gradtj)
!  if(associated(gradtk)) deallocate(gradtk)
  if(associated(vel))    deallocate(vel)
  if(associated(tem))    deallocate(tem)

  return
end subroutine euler
