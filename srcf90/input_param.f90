subroutine input_param
!
!  TRANSFER CONTROL PARAMETERS FROM ctrlib to typepara of new flow solver
!
  use typepara
  use ctrlib
  implicit none
  integer, parameter :: dp=kind(1.0d0)
!  real(dp) :: Ttmp, ptmp, refvisc, rat

! SET UP MULTI STAGE, NUMERICAL DISSIPATION, AND ITERATIONS  PARAMETERS
!
  if (associated(cstp)) deallocate(cstp)
  if (associated(cdis)) deallocate(cdis)
  if (associated(mcyc)) deallocate(mcyc)

  mstage  = imstage      
  mdissip = idissip
  mglevel = icycle
  
  psensor = pre_sensor
  mgtype  = multig_type

  allocate(cstp(mstage), mcyc(mglevel), cdis(mdissip) ) 
 
  cstp(1:mstage)  = rmstage(1:mstage)
  cdis(1:mdissip) = rdissip(1:mdissip)
  mcyc(1:mglevel) = rcycle(1:mglevel)

  cfl     = num_cfl
  cflf    = num_cflf
  cflc    = num_cflc
  vis2    = num_vis2
  vis4    = num_vis4
  smoo    = num_smoo
  
  vt = 1
  if(cfl<0.0d0) vt = 0    ! local time step
  cfl = abs(cfl)
  vis4 = vis4/32.0d0

  i2d=0

  gamma  = gas_gamma     !!1.40d0
  rm     = gas_mach      !! free stream Mach number
  alpha  = aoa_alpha     !! ANGLES OF ATTACK
  beta   = aoa_beta

!**********************************************************
! freestream velocity
  vel_fs = rm*sqrt(gamma*gas_pmean/gas_rhomean)
!**********************************************************

!-------------------- nondimensionalization---------------------
! rho0 = rho/rho_gas  = 1
! u0   = u/vel_fs     = 1
! p_dyna = 0.5*rho0*u0^2=0.5*gamma*p0*M^2===> p0 = 1/(gamma*M^2)
!---------------------------------------------------------------

!!!  rho0 = 1.27560704d0
  rho0 = 1.0d0
  p0   = 1.0d0/(gamma*rm**2.0d0)
  c0   = sqrt(gamma*p0/rho0)
  ei0  = p0/((gamma -1.0d0)*rho0)

  u0   = rm*c0 *cos(beta)*cos(alpha)
  v0   = rm*c0 *cos(beta)*sin(alpha)
  w0   = rm*c0 *sin(beta) 
  h0   = gamma*ei0 +0.5d0*(u0*u0 +v0*v0 +w0*w0)

! for unsteady case
!******************************************************************
!!!!  REMOVED FROM TYPEPARA
!!!!
!!!!  runtype = execmode !!jb%runtype             ! =0,steady solver; =1,unsteady solver.
!!!!  iflow   = i_flow                            ! =0 no initial flow; =1, input initial flow
!******************************************************************

  nperiod = fsi_period                        ! number of period 
  nstep   = fsi_step                          ! number of time interval in one period
  dtstep  = fsi_dt

!**********************************************************
!
! read from jobsheet
!
  eqntype = ceqntype
  dissip  = cdissip         ! matrix: matrix dissipation; scalar: scalar dissipation

!**********************************************************
! nondimensionalize time step
  dtstep = dtstep*vel_fs/ref_leng
  ! For transient flow, use a small time step size
  if(eqntype == 'n' .or. eqntype == 'N') then
    if( execmode == 1) dtstep = u0/ref_leng/5.0d0
  end if

! parameters used in upwind scheme
! read from jobsheet
  distype   = cdistype           ! r1: Roe scheme (1st order); r2: 2nd order; c: central difference
  limfac    = r2_limfac          ! limiter coefficient (Roe scheme)
  epsentr   = r2_epsentr         ! entropy correction coefficient (Roe's scheme)

  volref    = ref_leng**3.0d0  ! reference volume
  limref(1) = rho0
  limref(2) = sqrt(u0*u0+v0*v0+w0*w0)
  limref(3) = limref(2)
  limref(4) = limref(2)
  limref(5) = p0

  if(eqntype == 'n' .or. eqntype == 'N') then
    gas_reynolds = 500.0d0
    gas_mu       = rho0*u0*ref_leng/gas_reynolds
    viscl        = rho0*u0*ref_leng/gas_reynolds
  end if

!  ibm = ibm_on !.true.

!  print *," runtype,i2d :", runtype,i2d

  print *," mach number :", rm
  print *," gamma       :", gamma
  print *," alpha, beta :", alpha,beta
  print *," p0          :", p0
  print *," c0          :", c0
  print *," ei0         :", ei0
  print *," rho0        :", rho0
  print *," u, v, w     :", u0,v0,w0
  print *," h0          :", h0
  print *," dtstep      :", dtstep
  print *," PERIOD,nstep:", nperiod, nstep
  print *," MGLEVEL     :", mglevel
  print *," CFL         :", cfl
  print *," vt          :", vt
  print *," VIS2        :", vis2
  print *," VIS4        :", vis4
  print *," MCYC        :", mcyc
  print *," gas pmean   :", gas_pmean
  print *," gas rhomean :", gas_rhomean
  print *," gas tmean   :", gas_tmean
  print *," cstp        :", cstp
  print *," cdis        :", cdis
  print *," mstage, mdissip:", mstage, mdissip

  print *," multigrid cycle :", mgtype
  print *," pressure sensor :", psensor
  print *," smoothing coef  :", smoo
  print *," freestream vel  :", vel_fs
  print *," ref_leng        :", ref_leng

  print *," volref                      :", volref
  print *," limref                      :", limref
  print *," discretization type c/r2/r1 :", distype
  print *," dissip type scalar/matrix   :", dissip
  print *," limiter coeff (Roe)         :", limfac
  print *," entropy correction          :", epsentr
  print *," equation type, e/n          :", eqntype
  print *," laminar viscosity           :", viscl
  print *," Reynolds number             :", gas_reynolds
  print *," dynamic viscosity           :", gas_mu
  print *," laminar Prandtl number      :", Prlmean

  print *," Using IBM                   :", ibm_on
  print *," APPROXBC                    :", approxbc

  return

end subroutine input_param
 
