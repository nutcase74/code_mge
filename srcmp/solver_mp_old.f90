subroutine solver_mp(blk,project,job,model)


! This is the steady/unsteady Cartesian grid Euler solver subroutine.
! V-cycle multigrid technique is implemented on levels of embedded 
! Cartesian grids of different resolution.
! nperiod cycles are calculated with nstep time intervals in each cycle.
!   For steady case, nperiod=1 and nstep=1. mcyc(i) multigrid cycles
! are calculated for one steady state. 
!   For unsteady case, Jameson's dual time stepping method is implemented;
! finally, force coefficient are transformed to fourier series representation.
! ===============================================================================
!
!!! COPY FROM NSSOLVER FOR REFERENCE
!
!     multigrid ::
!
!     nlevel  the number of level in the multigrid cycle
!
!     mgtype    1     v cycle
!               2     w cycle
!
!     lres  true    compute only residual 
!           false   compute residual and update solution
!
!     lfir  true    save the initial values on first entry to a coarser in euler
!           false   
!
!     fcoll         relaxation factor for residual collection
!                   recommended value 1.
!
!     fadd          relaxation factor for multigrid corrections
!                   recommended value 1.
!
!     fbc    1.     update the far field on the coarse meshes
!            0.     freeze the far field on the coarse meshes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use typepara
  use typelib4
  use ctrlib
  use runlib
  use typelib5

  use codepak

! ----- MPI OPS
  use mplib
  use mpigrid

  implicit none

  include 'mpif.h'

  integer, parameter :: dp=kind(1.0d0)
  type(mgrid), pointer         :: blk
  type(family3d), pointer      :: project
  character(len=*), intent(in) :: job, model

  logical :: savenow
  logical :: lres, lfir

! local variables
  integer :: mgl, mlvl
  integer :: m,n,err
  integer :: ninner, ncyc
  integer :: period,rstep,itime

  type(mgrid), pointer :: btr
  integer :: nlevel, ierr
  integer, dimension(:), pointer :: mblks,icyc

!!! VARIABLE ADDED 07/12 for BLK2MGRID
  real*8   :: timer_on, timer_off

  integer  :: atcellcenter
  real(dp) :: p_tot

  include 'fconst.h'
  include 'mgridtopo.h'
  include 'debug.inc'

  call eval_popul(blk,mblks)
  nlevel = size(mblks)
  allocate(icyc(nlevel))

  itime = 0

  alop: do period=1, nperiod

    do rstep=1,nstep

    print *, "PERIOD, RSTEP ", period, rstep, " ASSO POSTGRID: ", associated(project%mgrida)

    if (runtype/=0) then

! -----  OBTAIN KINEMATICS OF MOVING SURFACE, BCAST FROM MASTER TO ALL WORKERS

      if (myrank==imaster) call mvsurf_coupled(blk, project%mgrida)
      call bcast_movbc(blk)
    end if

    do mlvl = 1, nlevel

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
! --- NUMBER OF RELAXATION ON EACH GRID
      ninner = 1

! --- TIME STEP ON LEVEL mlvl USING MULTIGRID
      mgl    = mlvl
      lop0: do ncyc = 1, mcyc(mlvl)             ! mcyc(mlvl) == number of iteration on current level

! ----- SOLUTION OF CURRENT LEVEL (mlvl) USING MULTIGRID CYCLC 

        icyc = 0
        lopmg:do

! ------  RESTRICTION AND INJECTION STARTS FROM INITIAL LEVEL UNTIL COARSEST GRID IS REACHED

          lopr:do

            icyc(mgl) = icyc(mgl) + 1

            lres = .FALSE.
            if (mgl==mlvl) then
              lfir = .FALSE.
            else
              lfir = .TRUE.
            end if

! --------- COMPUTE RESIDUAL AND UPDATE SOLUTION

            do n=1,ninner
              do m=1, mblks(mgl)
                btr=>get_mgrid(blk,mgl,m)
                call step(btr)
                call euler_mproc(btr, itime, ncyc, mgl, mlvl, lres, lfir)

!                call sendrecv_grid2(blk, xoverlap, cvar_dat, 1, mgl, nprocs, ierr)
!                call bcond_overlap(btr, mgl, nlevel, mblks)

              end do
!! change from 0 to 1
!! transfer partial grid causes diff with serial code results
!              call MPI_BARRIER(MPI_COMM_WORLD,ierr)

              timer_on = MPI_WTIME()
              call reset_up2date(blk)
              call sendrecv_grid(blk, xoverlap, cvar_dat,    mgl, 0, nprocs, ierr)
              !call sendrecv_grid(blk, xoverlap, phantom_dat, mgl, 0, nprocs, ierr)
              timer_off = MPI_WTIME()
              data_overlap = data_overlap + (timer_off-timer_on)

              timer_on = MPI_WTIME()
              call bcond_overlap(btr, mgl, nlevel, mblks)
              timer_off = MPI_WTIME()
              time_overlap = time_overlap + (timer_off-timer_on)

            end do

!*************************************************************************************
! OPTIONAL FOR DEBUG ::
! COMPARE CV OF A GRID ON DIFFERENT WORKERS
! 
!            call compare_grid(blk,  4, cvar_dat, ierr)
!            call compare_grid(blk,  5, cvar_dat, ierr)
!            call compare_grid(blk,  7, cvar_dat, ierr)
!            call compare_grid(blk, 14, cvar_dat, ierr)
!*************************************************************************************


            if (mgl == mlvl) then

! ----------- OUPUT CONVERGENCE INFO TO SCREEN ONLY ON BASE GRID

              do m=1,mblks(mgl)
                btr=> get_mgrid(blk, mgl, m)
                !print *, "Level= ", mgl,"ID=",btr%id, " NCYC=", ncyc, " RESID", btr%rtrms(ncyc)
                write(*,fmt="(4(a,I5),a,e16.10)") " IM=",myrank," Level=", mgl," ID=",btr%id, " NCYC=", ncyc,&
                                                 " RESID=", btr%rtrms(ncyc)

              end do
            end if
          
            if (ncyc >= mcyc(mlvl)) exit lop0

            if (mlvl == 1) cycle lop0
            
            ! ----- new residual (with W^(n+1))

            lres = .TRUE.
            lfir = .FALSE.

            do m=1,mblks(mgl)
              btr=>get_mgrid(blk,mgl,m)
              call euler_mproc(btr,itime,ncyc,mgl,mlvl,lres,lfir)

!              call sendrecv_grid2(blk, xoverlap, cvar_dat, mgl, nprocs, ierr)
!              call bcond_overlap(btr, mgl, nlevel, mblks)
            end do

!!            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!*************************************************************************************
! REMOVE TO SAVE MPI OPS
!
!            call sendrecv_grid(blk, xoverlap, cvar_dat, mgl, 1, nprocs, ierr)
!            call bcond_overlap(btr, mgl, nlevel, mblks)
!*************************************************************************************

! -----     RESTRICTION: transfer solution/residual to next coarse grid   
!           => movco => collc (w, dw)

            mgl = mgl - 1
            timer_on = MPI_WTIME()
            call sendrecv_grid(blk, xchild, restrict_dat, mgl, 0, nprocs, ierr)
            timer_off = MPI_WTIME()
            data_movco = data_movco + (timer_off-timer_on)

            timer_on = MPI_WTIME()
            call movco(blk,mgl,nlevel,mblks)
            timer_off = MPI_WTIME()
            time_movco = time_movco + (timer_off-timer_on)


            if (mgl == 1) exit lopr
          end do lopr

! ------  TIME STEP ON THE COARSEST GRID

!          if (myrank==imaster) print *," #########################  TIME STEP ON COARSEST GRID ", mgl

!          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

          lres = .FALSE.
          lfir = .TRUE.
          do n=1, ninner

            do m=1,mblks(mgl)
              btr=>get_mgrid(blk,mgl,m)
              call step(btr)
              call euler_mproc(btr,itime,ncyc,mgl,mlvl,lres,lfir)

!              call sendrecv_grid2(blk, xoverlap, cvar_dat, 0, mgl, nprocs, ierr)
!              call bcond_overlap(btr,mgl,nlevel,mblks)
            end do

!*************************************************************************************
! REMOVE TO SAVE MPI OPS
!
!            call sendrecv_grid(blk, xoverlap, cvar_dat, mgl, 1, nprocs, ierr)
!            call bcond_overlap(btr,mgl,nlevel,mblks)
!*************************************************************************************

          end do

          lopp:do
!
! -----     PROLOGATION: transfer solution correction from a coarse grid to next fine grid 
!           => movfin => addw  (w, ws) => emparent
!
            mgl = mgl + 1

            timer_on = MPI_WTIME()
            call sendrecv_grid(blk, xparent, prolong_dat, mgl, 0, nprocs, ierr)
            timer_off = MPI_WTIME()
            data_movfin = data_movfin + (timer_off-timer_on)

            timer_on = MPI_WTIME()
            call movfin(blk, mgl, mlvl, nlevel, mblks)
            timer_off = MPI_WTIME()
            time_movfin = time_movfin + (timer_off-timer_on)

            if (mgl == mlvl) cycle lop0
            if (icyc(mgl) /= mgtype) cycle lopmg
            icyc(mgl) = 0
          end do lopp
        end do lopmg
      end do lop0 

      if (mlvl < nlevel) then
! ----- PROLONGATION:  TRANSFER SOLUTION TO THE NEXT FINER GRID
!       MOVFIN => ADDW => emparent
        mgl = mlvl +1

        timer_on = MPI_WTIME()
        call sendrecv_grid(blk, xparent, prolong_dat, mgl, 0, nprocs, ierr)
        timer_off = MPI_WTIME()
        data_movfin = data_movfin + (timer_off-timer_on)

        timer_on = MPI_WTIME()
        call movfin(blk,mgl,mlvl,nlevel,mblks) 
        timer_off = MPI_WTIME()
        time_movfin = time_movfin + (timer_off-timer_on)

      end if

    end do ! do mlvl=1,nlevel


!! -----  OPTION: COMPARE BLK'S VALUES ********************
!    call compare_grid(blk,  5, cvar_dat, ierr)
!    call compare_grid(blk,  7, cvar_dat, ierr)
!    call compare_grid(blk, 14, cvar_dat, ierr)
!!*********************************************************

    if (runtype>0) then
      !!
      !! update the flow variables of the previous two physical time step
      !!
      do mgl=1, nlevel
        do m = 1, mblks(mgl)
          btr=>get_mgrid(blk,mgl,m)
          btr%w2 = btr%w1
          btr%w1 = btr%w
        end do
      end do
    end if

! ------  ALL WORKERS SEND DATA TO MASTER
    call sendto_grid(blk, imaster, cvar_dat, ierr)
    if (myrank==imaster) then
      p_tot = gas_pdyna
      call blk2mgrid(blk, project%mgrida, nlevel, mblks, p_tot, atcellcenter)
      call step_fsijob(project%mgrida, project%groupa, project%csgroupa, project%ctrlvol, &
                       cfd_dof, csd_dof, &
                       fsi_tt, fsi_dt, fsi_f2s, job, model, savenow, runtype, err)
    end if 

    itime = itime +1
    if (runtype==0) exit alop

    end do
  end do alop

  if (associated(mblks)) deallocate(mblks)
  if (associated(icyc)) deallocate(icyc)
  return
end subroutine solver_mp
