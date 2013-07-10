program main_mproc

  use typelib
  use jobpak
  use runlib
  use ctrlib
  use mgridlib
  use filelib
!
! CODELIB provides interface between block_type & mgrid
!
  use codepak
!
! MPI-RELATED MODULES
!
  use mpilib
  use mproclib

  implicit none

  integer, parameter :: dp=kind(1.0d0)
  character(len=64)  :: job, model
  type(jobsheet)     :: jb

  type(family3d), pointer :: project
  type(mgrid), pointer    :: blk
  integer, dimension(:), pointer :: popul

  integer :: err_job, err_cloud, error
  integer :: mxmcyc
  integer :: i, j, itmp
  logical :: preponly, rom_on

!! MPI VARIABLES
  real*8   :: tstart, tend
  integer  :: ierr, taskid, nproc
  character(len=32) :: text

  include 'debug.inc'
  include 'mgridtopo.h'
  include 'mpif.h'

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr)

  myrank      = taskid
  nprocs      = nproc
  debug       = .false.
  debug_omp   = .false.
  silent_mode = (myrank/=imaster)  !! MGRID PREPROCESS GOES SILENT FOR ALL WORKER EXCEPT MASTER

  nullify(project, blk,  popul)
  call make_jobsheet(jb)

!  print *, "  CARTESIAN GRID EULER SOLVER   "
!  print *, "  taskid = ", taskid, " nproc = ", nproc

  if (myrank == imaster) then
    !
    ! --- READ INPUT FROM JOBSHEET, PREPARE GEOMETRY & GRID 
    !
    call make_fsijob(project, "jobsheet", job, model, preponly, rom_on, jb, err_job)
    mxmcyc = maxval(jb%mcyc)

    if (err_job/=0) then
      !
      ! --- ABORT FULL-ORDER SIMULATIONS WHEN ROM IS REQUIRED
      !
      print *," ERROR IN MAKE_FSIJOB, CHECK JOBSHEET SETTINGS FOR CONFLICTS/ERRORS  "
      call clear_fsijob(myproj)
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif

    if (rom_on) then
      !
      ! --- ABORT FULL-ORDER SIMULATIONS WHEN ROM IS REQUIRED
      !
      call clear_fsijob(myproj)
      call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    endif
  else
    print *," Iam ", myrank
    call make_family(project) 
    print *, " Associated(FM) ",  associated(project)
  end if

  !! DEBUG OPTION: ASSIGN ALL MESH TO ALL WORKERS, all2all=.false. by default
  !! all2all = .true.

  call MPI_BCAST(job,    len(job)  ,MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(model , len(model),MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mxmcyc, 1,         MPI_INTEGER,  imaster, MPI_COMM_WORLD, ierr)

  !!
  !! DISTRIBUTE MGRID TOPOLOGIES
  !!
  call bcast_gridtopo(project%mgrida)

!  call MPI_FINALIZE(ierr)
!  stop
  !!
  !! PERFORM CPU-LOAD DISTRIBUTION 
  !!
  call mpidistr_grid (project%mgrida, nprocs)

  print *," **********************************************"
!  print *," Im ",myrank," MYLOAD=", myload
  print *," Im ",myrank," IFLOW =", iflow
  print *," **********************************************"

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  call op_create_maxvec

  blk    => project%mgrida
  call eval_popul(blk, popul)
!
! --- SET MPLIB VARIABLES
!
  num_levels = size(popul)
  num_grids  = sum(popul)
!
! ----- HANDLE PARAMETERS INPUT
!
  if (myrank == imaster) call input_param
  call bcast_input
!
! --- ALLOCATE STACK FOR FLOW VARIABLES
!
  call make_flow(blk, mxmcyc,myrank)
!
! --- CALC METRICS
!
  call init_metric(blk)
!
! ----- INITIALISE FLOW 
!
  call init_flow_mp(blk, (iflow==1), job)
!
  print *,"====================================================="
  print *
  print *," (fsi_job)    ", trim(job)
  print *," (fsi_model)  ", trim(model)
  print *," TIME STEP IN FSI COMPUTATION, dtstep :",dtstep
  print *," TIME STEP RATIO F/S                  :",fsi_f2s
  print *," FSI_TT                               :",fsi_tt
  print *," FSI_dT                               :",fsi_dt
  print *," FSI_f2s                              :",fsi_f2s
  print *," ASSOCIATED POSTGRID                  :", associated(project%mgrida)
  print *," EXECMMODE                            :",execmode
  print *," EQNTYPE                              :",eqntype
  print *," DISTYPE                              :",distype
  print *," DISSIP                               :",dissip
  print *," LIMFAC                               :",limfac
  print *," VOLREF                               :",volref
  print *," LIMREF                               :",limref
  print *," EPSENTR                              :",epsentr
  print *," VT                                   :",vt
  print *," GAMMA                                :",gamma
  print *," rho0                                 :",rho0
  print *," GAS_TMEAN                            :",gas_tmean
  print *," PRLMEAN                              :",prlmean
  print *
  print *,"====================================================="


!! OPTION:: TURN OFF bcond_overlap for debug
!  bcond_overlap_off=.false.

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

! --- SET TIMER
!  call reset_mpi_timer
  tstart = MPI_WTIME()

  itmp = count_mgrid(blk)
  call reset_timer(itmp)

  call solver_mp(blk,project,job,model)

  tend = MPI_WTIME()

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
! ----- SEND SOLUTION TO MASTER FOR FILES OUTPUT
!
  call sendto_grid(blk,imaster,cvar_dat,error)
  if (myrank==imaster)  call write_solut(blk, job, error)

! --- OPTION: WRITE SOLUTION OF EACH EVERY FLOW BLOCK IN INDIVIDUAL FILES
!  if (myrank==imaster)  call write_solut_one(blk, job, error)

!
! --- PREPARE TO QUIT
!
  if (associated(popul)) deallocate(popul)
  call clear_jobsheet(jb)
  call clear_fsijob(project)
  call clear_ctrl
  call clear_param
!
! ----- EXIT MPI ENVIRONMENT
!
  print *," TIME ELASPED =",tend - tstart
!  call display_mpi_timer

  call print_timer
  call mpi_finalize(ierr)

end program main_mproc 
