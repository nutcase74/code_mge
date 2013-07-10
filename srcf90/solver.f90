subroutine solver(blk,project,job,model)

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
!     nlevel         number of level in the multigrid cycle
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
  use typelib
  use ctrlib
  use runlib
!  use typelib5

  implicit none
  integer, parameter :: dp=kind(1.0d0)
  type(mgrid),         pointer :: blk
  type(family3d),      pointer :: project
  character(len=*), intent(in) :: job, model
!
! local variables
!
  logical :: savenow, lwork
  logical :: lres, lfir

  integer :: mgl, mlvl
  integer :: m,n,err,io
  integer :: ninner, ncyc, nlevel
  integer :: period,rstep,itime
  
  integer, dimension(:), pointer :: mblks,icyc

!!! VARIABLE ADDED 07/12 for BLK2MGRID
  integer :: atcellcenter
  real(dp) :: max_res, ini_res
  type(mgrid), pointer :: btr
  character(len=64) :: fname
!
! VARIABLES FOR KEEPING TIME
!
  real(dp) :: waltime
  real*8   :: tstart, tfinish
  integer  :: save_interval
  logical  :: overwrt

  include 'fconst.h'
  include 'timer.inc'

  nullify(mblks,icyc)
!
! GET MULTI-LEVEL EMBEDDED MESH TOPOLOGY
!
  call eval_popul(blk, mblks)
  nlevel = size(mblks)
  allocate(icyc(nlevel))

! 
! SET LWORK TRUE FOR SERIAL CODE
!
  lwork = .true.
!
! SET TIME (cycle) INTERVAL TO SAVE THE SOLUTION ON THE FINEST LEVEL
!
  save_interval = 1000000

  fname = trim(job)//"_res.dat"
  call channel(io)
  open (unit=io, file=fname, status="unknown", action="write", form="formatted")
!
! --- START CLOCK
!
  call wall_time(tstart)

!
! START TIME MARCHING
!
  itime = 0
  alop: do period=1, nperiod

    do rstep=1, nstep

    print *, "PERIOD, RSTEP ", period, rstep , " ASSO POSTGRID: ", associated(project%mgrida)

    !! OBTAIN KINEMATICS OF MOVING SURFACE 

    if (approxbc) call movsurf(blk) !!, project%mgrida)

    do mlvl=1, nlevel
    
      ninner = 1                         ! number of relaxation on each grid
      mgl    = mlvl

      lop0: do ncyc = 1,mcyc(mlvl)       ! number of multigrid cycles

        icyc = 0                         ! type of multigrid cycle, 1-v;2-w cycle
        lopmg:do
          lopr:do

            icyc(mgl) = icyc(mgl) + 1

            !if(mgl.eq.mlvl-1) then
            !  cfl = abs(cflf)
            !else
            !  cfl = abs(cflc)
            !end if

!            print *, "r ncyc/icyc/mgl/mlvl",ncyc, icyc,mgl,mlvl

            !------RESTRICTION-----------------------------------------------------------
            !
            lres = .FALSE.
            if(mgl .eq. mlvl) then
               lfir = .FALSE.
            else
               lfir = .TRUE.
            end if

            ! -- COMPUTE RESIDUAL AND UPDATE SOLUTION

            do n=1,ninner
              do m=1, mblks(mgl)
                btr=>get_mgrid(blk, mgl, m)
                call timestep  (btr, lwork )
                call euler     (btr, itime, ncyc, mgl, mlvl, lres, lfir, lwork)
                call bc_overlap(blk, mgl)
              end do
            end do

            if (mgl .eq. mlvl) then

              ! --- OUTPUT RESIDUAL TO SCREEN ONLY ON THE FINEST GRID

              max_res = 1000.0d0
              do m=1,mblks(mgl)

                btr=> get_mgrid(blk, mgl, m)

                write(*,fmt="(3(a,I6),a,e18.10)") " Level=", mgl,"   ID=",btr%id, &
                            "   NCYC=", ncyc,"   RESID= ", btr%rtrms(ncyc)
                ! save the maxmum residual for each cycle 
                ! save the first residual for each level
                if(ncyc == 1)  ini_res = btr%rtrms(ncyc)
                if(btr%rtrms(ncyc) < max_res ) then
                  max_res = btr%rtrms(ncyc)
                end if

                ! check NaN
                !
                if (btr%rtrms(ncyc) .ne. btr%rtrms(ncyc)) stop
              end do
 
              ! --- WRITE MAXIMUM RESIDUAL TO FILE

              write(io, fmt="(2i6, e18.10)") mgl, ncyc, max_res
            end if
            
!!@            if(ncyc >= mcyc(mlvl)) exit lop0

!            if(mlvl == 1 .and. ini_res/max_res > 1000000d0) exit lop0
            if(mlvl == 1) cycle lop0
    
            ! ----- new residual (with W^(n+1))
            lres = .TRUE.
            lfir = .FALSE.

            do m=1,mblks(mgl)
              btr=>get_mgrid(blk,mgl,m)
              call euler(btr,itime,ncyc,mgl,mlvl,lres,lfir, lwork)
              call bc_overlap(blk, mgl)
            end do

            mgl = mgl - 1
            call movco(blk, mgl, nlevel, mblks)   ! transfer solution to the coarser grid
           
            if (mgl == 1) exit lopr
          end do lopr

          ! ------prolongation------------------------------------------------------
          lres = .FALSE.
          lfir = .TRUE.

          lopp:do
!            print *, "p ncyc/icyc/mgl/mlvl",ncyc, icyc,mgl,mlvl
            do n=1,ninner
              do m=1,mblks(mgl)
                btr=>get_mgrid(blk,mgl,m)
                call timestep( btr, lwork )
                call euler   ( btr,itime,ncyc,mgl,mlvl,lres,lfir, lwork)
                call bc_overlap(blk, mgl)
              end do
            end do
            lres = .FALSE.
            lfir = .FALSE.


            mgl = mgl + 1
            call movfin(blk, mgl, mlvl, nlevel, mblks)
            if (mgl.eq.mlvl) exit lopmg

            if (icyc(mgl) /=mgtype) cycle lopmg
            icyc(mgl) = 0
          end do lopp
        end do lopmg
!stop
        ! if the convergence criterion is reached, end the loop 
        if (mgl .eq. mlvl) then
           if(ini_res/max_res > 1000000d0) exit lop0
        end if

!
! SAVE INTERIM RESULT FOR DEBUG
!
        if (mlvl==nlevel.and. modulo(ncyc,save_interval)==0) then
        !if (modulo(ncyc,save_interval)==0) then
          overwrt = .false. !! NOT OVERWRITE MESH SOL
          fsi_tt = itime*dtstep
          fsi_dt = dtstep
          call upload_solut(blk, gas_ptot, atcellcenter)
          call step_fsijob(project%mgrida, project%groupa, project%csgroupa, project%ctrlvol, &
               cfd_dof, csd_dof, fsi_tt, fsi_dt, fsi_f2s, job, model, savenow, execmode, err, overwrt)
        end if
 
      end do lop0 

      if (mlvl .lt. nlevel ) then
        mgl = mlvl +1
        call movfin(blk,mgl,mlvl,nlevel,mblks) ! transfer solution to the next finer grid
      end if

    end do

    if (execmode==1) then
      !!
      !! FOR UNSTEADY FLOW :: UPDATE FLOW VARIABLES OF THE PREVIOUS TWO PHYSICAL TIME STEPS
      !!
      do mgl = 1,nlevel
        do m = 1,mblks(mgl)

          btr    => get_mgrid(blk,mgl,m)
          btr%w2 = btr%w1
          btr%w1 = btr%w

        end do
      end do
    end if

    overwrt = .false. !! NOT OVERWRITE MESH SOL
    call upload_solut(blk, gas_ptot, atcellcenter)
    fsi_tt = itime*dtstep
    fsi_dt = dtstep
    call step_fsijob(project%mgrida, project%groupa, project%csgroupa, project%ctrlvol, &
         cfd_dof, csd_dof, fsi_tt, fsi_dt, fsi_f2s, job, model, savenow, execmode, err, overwrt)


    itime = itime +1
    if (execmode==0) exit alop

    write(unit=io,fmt="(I6)") mgl

    end do
  end do alop

  close(unit=io, status="keep", iostat=err)

  if (associated(mblks)) deallocate(mblks)
  if (associated(icyc) ) deallocate(icyc)

  call wall_time(tfinish)
  op_solver = op_solver + (tfinish - tstart)

  return
end subroutine solver
