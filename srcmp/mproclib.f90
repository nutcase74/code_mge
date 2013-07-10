module mproclib

  use codepak
  use typepara
  implicit none
  private

  public :: init_flow_mp
  public :: bcast_input
  public :: solver_mp

contains

  include 'init_flow_mp.f90'
  include 'bcast_input.f90'
  include 'solver_mp.f90'

end module mproclib
