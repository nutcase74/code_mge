subroutine solver_mp(blk,project,job,model)

! This is the steady/unsteady Cartesian grid Euler solver subroutine.
! V-cycle multigrid technique is implemented on levels of embedded 
! Cartesian grids of different resolution.
! nperiod cycles are calculated with nstep time intervals in each cycle.
!   For steady case, nperiod=1 and nstep=1. mcyc(i) multigrid cycles
! are calculated for one steady state. 
!   For unsteady case, Jameson's dual time stepping method is implemented;
! finally, force coefficient are transformed to fourier series representation.
!
! ===============================================================================
!
!!! COPY FROM NSSOLVER FOR REFERENCE
!
!     multigrid ::
!
!     nmesh         the number of meshes in the multigrid cycle
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
  use mpilib


  implicit none
  include 'mpif.h'

  integer, parameter :: dp=kind(1.0d0)
  type(mgrid),      pointer :: blk
  type(family3d),   pointer :: project
  character(len=*), intent(in) :: job, model
  logical :: savenow
  logical :: lres, lfir

! local variables
  integer :: mgl, mlvl
  integer :: m,n,err
  integer :: ninner, ncyc
  integer :: period,rstep,itime

!!! VARIABLE ADDED 07/12 for BLK2MGRID
  integer :: atcellcenter, ierr, nlevel
  type(mgrid), pointer :: btr
  integer, dimension(:), pointer :: mblks,icyc
!
! VARIABLES FOR KEEPING TIME
!
  real*8   :: tstart, tfinish
  integer  :: save_interval
  logical  :: overwrt

  include 'fconst.h'
  include 'mgridtopo.h'
  include 'timer.inc'

  nullify(mblks,icyc)
  save_interval = 1000000

  call eval_popul(blk, mblks)
  nlevel = size(mblks)
  allocate(icyc(nlevel))

!
! --- START CLOCK
!
  call wall_time(tstart)

  itime = 0
  alop: &
  do period=1, nperiod
    do rstep=1,nstep

      print *,"IM",myrank," PERIOD, RSTEP ", period, rstep, " ASSO POSTGRID: ", associated(project%mgrida), &
              " APPROXBC ", approxbc

      savenow = ( ((period==nperiod) .and. (rstep==nstep)) .or. execmode==0)
!
! -----  OBTAIN KINEMATICS OF MOVING SURFACE, BCAST FROM MASTER TO ALL WORKERS
!
      if (approxbc) then
        if (myrank==imaster) call movsurf(blk)
        call bcast_movbc(blk)
      end if

      do mlvl=1,nlevel
    
        ninner = 1                         ! number of relaxation on each grid
        mgl    = mlvl

        lop0: do ncyc = 1,mcyc(mlvl)             ! number of multigrid cycles

          icyc = 0

          lopmg:do
            lopr:do

              icyc(mgl) = icyc(mgl) + 1

              !------restriction-----------------------------------------------------------
              lres = .FALSE.
              if (mgl == mlvl) then
                lfir = .FALSE.
              else
                lfir = .TRUE.
              end if

              ! compute residual and update solution

              do n=1,ninner
                do m=1, mblks(mgl)
                  btr=>get_mgrid(blk, mgl, m)
                  call timestep  (btr, (count(myload==btr%id)/=0) ) !, myrank )
                  call euler     (btr, itime, ncyc, mgl, mlvl, lres, lfir, (count(myload==btr%id)/=0) ) !, myrank )

                  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

                  call sendrecv_grid(blk, xoverlap, cvar_dat, mgl, ierr)
                  call bc_overlap(blk, mgl)
                end do
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

            if (mgl .eq. mlvl) then

! ----------- OUPUT CONVERGENCE INFO TO SCREEN ONLY ON BASE GRID

              do m=1,mblks(mgl)
                btr=> get_mgrid(blk, mgl, m)
                !print *, "Level= ", mgl, " NCYC = ", ncyc, "residual", btr%rtrms(ncyc)
                ! check NaN
                !if (btr%rtrms(ncyc) .ne. btr%rtrms(ncyc)) stop
                if(count(myload==btr%id)/=0) &
                write(*,fmt="(4(a,I5),a,e16.10)")" IM=",myrank, " Level=", mgl," ID=",btr%id, " NCYC=", ncyc,&
                                                 " RESID= ", btr%rtrms(ncyc)

              end do
            end if
            
            if (ncyc >= mcyc(mlvl)) exit lop0
            if (mlvl == 1) cycle lop0
            
            ! ----- new residual (with W^(n+1))
            lres = .TRUE.
            lfir = .FALSE.

            do m=1,mblks(mgl)
              btr=>get_mgrid(blk,mgl,m)
              call euler     (btr,itime,ncyc,mgl,mlvl,lres,lfir,(count(myload==btr%id)/=0)) !, myrank)

              !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

              call sendrecv_grid(blk, xoverlap, cvar_dat, mgl, ierr)
              call bc_overlap(blk, mgl)
            end do

            mgl = mgl - 1

            !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            call sendrecv_grid(blk, xchild, restrict_dat, mgl, ierr)
            call movco(blk, mgl, nlevel, mblks)   ! transfer solution to the coarser grid
            
            if(mgl == 1) exit lopr
          end do lopr

          ! ------prolongation------------------------------------------------------
          lres = .FALSE.
          lfir = .TRUE.

          lopp:do
            do n=1,ninner
              do m=1,mblks(mgl)
                btr=>get_mgrid(blk,mgl,m)

                call timestep( btr, (count(myload==btr%id)/=0) ) !, myrank )
                call euler   ( btr, itime, ncyc, mgl, mlvl, lres, lfir, (count(myload==btr%id)/=0)) !, myrank )
                call sendrecv_grid(blk, xoverlap, cvar_dat, mgl, ierr)
                call bc_overlap(blk, mgl)

              end do
            end do

            lres = .FALSE.
            lfir = .FALSE.

            mgl = mgl + 1

            call sendrecv_grid(blk, xparent, prolong_dat, mgl, ierr)
            call movfin(blk,mgl,mlvl,nlevel, mblks)
            if (mgl.eq.mlvl) exit lopmg
            if (icyc(mgl) /=mgtype) cycle lopmg
            icyc(mgl) = 0
          end do lopp
        end do lopmg

!
! SAVE INTERIM RESULT FOR DEBUG
!
        if (myrank==imaster) then
          if (mlvl==nlevel.and. modulo(ncyc,save_interval)==0) then
            overwrt = .false. !! NOT OVERWRITE MESH SOL
            call upload_solut(blk, gas_pdyna, atcellcenter)
            call step_fsijob(project%mgrida, project%groupa, project%csgroupa, project%ctrlvol, &
               cfd_dof, csd_dof, fsi_tt, fsi_dt, fsi_f2s, job, model, savenow, execmode, err, overwrt)
          end if
        end if

      end do lop0 

      if (mlvl.lt.nlevel) then
        mgl = mlvl +1
        call sendrecv_grid(blk, xparent, prolong_dat, mgl, ierr)
        call movfin(blk,mgl,mlvl,nlevel,mblks) ! transfer solution to the next finer grid
      end if

    end do ! do mlvl=1,nmesh

    if (execmode==1) then
      !!
      !! update the flow variables of the previous two physical time step
      !!
      do mgl = 1,nlevel
        do m = 1,mblks(mgl)

          btr=>get_mgrid(blk,mgl,m)

          btr%w2 = btr%w1
          btr%w1 = btr%w

        end do
      end do
    end if

! ------  ALL WORKERS SEND DATA TO MASTER
      call sendto_grid(blk, imaster, cvar_dat, ierr)

      if (myrank==imaster) then
        call upload_solut(blk, gas_pdyna, atcellcenter)
        call step_fsijob(project%mgrida, project%groupa, project%csgroupa, project%ctrlvol, &
                       cfd_dof, csd_dof, &
                       fsi_tt, fsi_dt, fsi_f2s, job, model, savenow, execmode, err)
      end if

      itime = itime +1
      if (execmode==0) exit alop

    end do
  end do alop

  if (associated(mblks)) deallocate(mblks)
  if (associated(icyc)) deallocate(icyc)

  call wall_time(tfinish)
  op_solver = op_solver + (tfinish - tstart)

  return
end subroutine solver_mp
