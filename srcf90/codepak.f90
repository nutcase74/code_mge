module codepak

  use roepak
  use viscpak
  use bclib
  use typepara
  use ctrlib

  implicit none
  public 

  public  :: input_param
  public  :: init_metric
  public  :: init_flow
  public  :: make_flow

  public  :: movsurf
  public  :: movco
  public  :: movfin
  public  :: addw
  public  :: collc

  public  :: dflux
  public  :: eflux
  public  :: dfluxb
  public  :: efluxb
  private :: ceflux

!  public  :: vflux
!  public  :: vfluxb
!  public  :: fgradi 
!  public  :: fgradj 
!  public  :: fgradk 

!  public  :: dflux_roe2
!  public  :: eflux_roe2
!  private :: roe2eflux
!  private :: roe2dflux

!  public  :: dflux_roe1
!  public  :: eflux_roe1
!  private :: roe1eflux
!  private :: roe1dflux

  public  :: psmoo
  public  :: psmoo_old

  public  :: solver
  public  :: timestep
  public  :: euler

  public  :: write_solut
  public  :: read_solut
  public  :: write_solut_one

!  public  :: check_cutcell

  public  :: upload_solut
  private :: progon3

  public  :: force_int
  public  :: bc_pressure

  logical, public, save :: bcond_overlap_off=.false.

contains

  include 'init_metric.f90'
  include 'init_flow.f90'
  include 'make_flow.f90'
  include 'input_param.f90'

  include 'solver.f90'
  include 'timestep.f90'
  include 'euler.f90'

  include 'upload_solut.f90'
  include 'psmoo.f90'
  include 'wsolut.f90'
  include 'rsolut.f90'

  include 'movco.f90'
  include 'movfin.f90'
  include 'addw.f90'
  include 'collc.f90'

  include 'ceflux.f90'
  include 'dflux.f90'
  include 'eflux.f90'
  include 'efluxb.f90'
  include 'dfluxb.f90'

!  include 'var_diffs.f90'
!  include 'roe2eflux.f90'
!  include 'eflux_roe2.f90'
!  include 'efluxb_roe2.f90'
!  include 'roe2dflux.f90'
!  include 'dflux_roe2.f90'
!  include 'dfluxb_roe2.f90'

!  include 'eflux_roe1.f90'
!  include 'efluxb_roe1.f90'
!  include 'dflux_roe1.f90'
!  include 'dfluxb_roe1.f90'

!  include 'vflux.f90'
!  include 'vfluxb.f90'
!  include 'fgradi.f90'
!  include 'fgradj.f90'
!  include 'fgradk.f90'

!  include 'check_cutcell.f90'

  include 'movsurf.f90'
  include 'progon3.f90'
  include 'force_int.f90'

end module codepak

