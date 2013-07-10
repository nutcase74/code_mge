subroutine bcast_input
!
! USE FOR TRANSFER INPUT PARAMETERS FOR OLD VERSION FLOW_SOLVER
!
! NOTE: A MUST TO HAVE  mstage = 2
!
!
  use ctrlib
  use typepara
  use mpilib

  implicit none

  include 'mpif.h'

  integer :: ierr, sz_cstp, sz_mcyc, sz_cdis

  sz_cstp=0
  sz_mcyc=0
  sz_cdis=0

  if (myrank==imaster) then
    print *, "MSTAGE  ",mstage,  "CSTP :",associated(cstp)
    print *, "mglevel ",mglevel, "MCYC :",associated(mcyc)
    print *, "MDISSIP ",mdissip, "CDIS :",associated(cdis)
    if (associated(cstp)) sz_cstp = size(cstp)
    if (associated(mcyc)) sz_mcyc = size(mcyc)
    if (associated(cdis)) sz_cdis = size(cdis)
  endif
!
! PARAMETERS DEFINED IN CTRLIB
!
  call MPI_BCAST(execmode,1,         MPI_INTEGER,  imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(iflow,   1,         MPI_INTEGER,  imaster, MPI_COMM_WORLD, ierr)
!  call MPI_BCAST(job,     len(job)  ,MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
!  call MPI_BCAST(model ,  len(model),MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(approxbc,1,         MPI_LOGICAL,  imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ibm_on,  1,         MPI_LOGICAL,  imaster, MPI_COMM_WORLD, ierr)
!
!
  call MPI_BCAST(mstage,  1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mglevel, 1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mdissip, 1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nperiod, 1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(nstep,   1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(i2d,     1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(psensor, 1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(mgtype,  1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(vt,      1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(gamma,   1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rm,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
!  call MPI_BCAST(al,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(rho0,    1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(p0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(c0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(ei0,     1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(u0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(v0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(w0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(h0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

! removed ptot 2012-08-06
!  call MPI_BCAST(ptot,    1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(vel_fs,  1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
!! OBSOLETE chord0   call MPI_BCAST(chord0,  1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(kappa,   1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

! removed am, a0 2012-08-06
!  call MPI_BCAST(am,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
!  call MPI_BCAST(a0,      1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(alpha,   1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(beta,    1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(cfl,     1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(vis2,    1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(vis4,    1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dtstep,  1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(smoo,    1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(sz_cstp,  1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(sz_mcyc,  1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(sz_cdis,  1, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)

  if (myrank/=imaster) then
    print *,"Im ",myrank," MSTAGE ",mstage 
    if (sz_cstp>0) allocate( cstp(sz_cstp) )
    if (sz_mcyc>0) allocate( mcyc(sz_mcyc) )
    if (sz_cdis>0) allocate( cdis(sz_cdis) )
  end if

  if (sz_cstp>0) &
  call MPI_BCAST(cstp,  sz_cstp, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

  if (sz_cdis>0) &
  call MPI_BCAST(cdis,  sz_cdis, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

  if (sz_mcyc>0) &
  call MPI_BCAST(mcyc,  sz_mcyc, MPI_INTEGER, imaster, MPI_COMM_WORLD, ierr)

!! IMPORTANT NOTE: OLD VERSION OF EUCART REQUIRE MDISSIP=2, AND NO CDIS is required
!  call MPI_BCAST(cdis,  size(cdis), MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

!
!  FOLLOWNG LINES ADDED ON 6-AUG-2012 
!
  call MPI_BCAST(eqntype,  len(eqntype), MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(distype,  len(distype), MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(dissip,   len(dissip),  MPI_CHARACTER,imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(limfac,   1,            MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(volref,   1,            MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(limref,   size(limref), MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(epsentr,  1,            MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(viscl,    1,            MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(prlmean,  1,            MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(prtmean,  1,            MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(gas_tmean,1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(gas_const,1, MPI_DOUBLE_PRECISION, imaster, MPI_COMM_WORLD, ierr)

  return
end subroutine bcast_input

!  real(dp), public, save  :: gas_mach, gas_reynolds, prlmean, prtmean
!  real(dp), public, save  :: gas_pdyna
!  real(dp), public, save  :: gas_tmean, gas_rhomean, gas_pmean
!  real(dp), public, save  :: gas_ttot, gas_rhotot, gas_ptot
!  real(dp), public, save  :: gas_mu, gas_const, gas_sos, gas_cv, gas_cp, gas_gamma
!  ! free-stream enthalpy, speed
!  real(dp), public, save  :: gas_h0, gas_q0
