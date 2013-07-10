program main

!***********************************************************************
!
! SOLVE EULER SET OF EQUATIONS OF MOTION 
!
!***********************************************************************
  use typelib
  use jobpak
  use runlib
  use ctrlib
  use mgridlib
!
! CODELIB provides interface between block_type & mgrid
!
  use codepak

  implicit none

  integer, parameter :: dp=kind(1.0d0)
  character(len=64)  :: job, model, fname, fnmesh
  type(jobsheet)     :: jb

  type(family3d), pointer :: project
  type(mgrid), pointer    :: blk

  integer  :: err_job, err_cloud, error
  integer  :: mxmcyc
  integer  :: itmp
  logical  :: preponly, rom_on
  real     :: finish, start

  logical :: savenow
  integer :: err

  include 'debug.inc'
  include 'mgridtopo.h'
  include 'timer.inc'


  debug = .false.
  nullify(project, blk)

  call make_fsijob(project, "jobsheet", job, model, preponly, rom_on, jb, err_job)

  if (rom_on) then
    call clear_fsijob(project)
    stop
  end if

  blk => project%mgrida
!  call mgrid_display(blk)
!
! GET INTERPOLATION COEFFICIENT FOR CLOUD POINTS
!
  call cpu_time(start)
  call solid_ipolcoef(blk, default_iorder, default_ncloud, default_where2find, err_cloud)
  call cpu_time(finish)
  print '("Time = ",f16.8," seconds  SOLID IPOLCOEF ")',finish-start

  mxmcyc = maxval(jb%mcyc)
!
! --- RETRIEVE PARAMETERS
!
  call input_param
!
! --- ALLOCATE STACK FOR FLOW VARIABLES
!
  call make_flow(blk, mxmcyc)
!
! --- GRID METRIC
!
  call init_metric(blk)
!
! --- INITIALIZE FLOW FIELD
!
!$ZXQ  call init_flow(blk, job)
  call init_flow(blk, (iflow==1), job)

  print *,"====================================================="
  print *
  print *," (fsi_job)    ", trim(job)
  print *," (fsi_model)  ",trim(model)
  print *," TIME STEP IN FSI COMPUTATION, dtstep :",dtstep
  print *," TIME STEP RATIO F/S                  :",fsi_f2s
  print *," FSI_TT                               :",fsi_tt
  print *," FSI_dT                               :",fsi_dt
  print *," FSI_f2s                              :",fsi_f2s
  print *," ASSOCIATED POSTGRID                  :", associated(project%mgrida)
  print *
  print *,"====================================================="

!
! RESET TIMER FOR itmp NUMBER OF GRIDS
!
  itmp = count_mgrid(blk)
  call reset_timer(itmp)

  call cpu_time(start)
  call solver(blk, project, job, model)
  call cpu_time(finish)
  print '("Time = ",f16.8," seconds  SOLVER ")',finish-start
  call print_timer

! ---  WRITE SOLUTION TO FILE FOR RESTARTING
!
  call write_solut(blk, job, error)
!
! --- RELEASE STACK
!
  call clear_jobsheet(jb)
  call clear_fsijob(project)
  call clear_param

end program main 
