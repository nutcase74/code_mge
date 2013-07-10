module bclib

  use typepara
  implicit none
  private

  public  :: bcond
  public  :: bc_far
  public  :: bc_symm
  public  :: bc_periodic
  public  :: bc_overlap
  public  :: bc_cutcell
  public  :: bc_embed
  public  :: bc_corner
  public  :: bc_wall

  public  :: depvars
  private :: depvars_all, depvars_one

  interface depvars
    module procedure depvars_all, depvars_one
  end interface

contains

  include 'bcond.f90'
  include 'bc_overlap.f90'
  include 'bc_far.f90'
  include 'bc_symm.f90'
  include 'bc_periodic.f90'
  include 'bc_cutcell.f90'
  include 'bc_embed.f90'
  include 'depvars.f90'
  include 'bc_corner.f90'
  include 'bc_wall.f90'

end module bclib
