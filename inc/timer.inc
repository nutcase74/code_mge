  integer, parameter :: mxg=50
  integer :: whoami, num_of_grid, ops_grid(mxg), ops_level(mxg)
  real    :: time_grid(mxg), time_addw(mxg)
  real    :: op_addw, op_collc, op_euler, op_solver, op_step
  real    :: op_movfin, op_movco, op_overlap, op_overhead
  real    :: op_gradface, op_fluxvis
  real    :: dat_movfin, dat_movco, dat_overlap
  real    :: snd_movfin, snd_movco, snd_overlap
  real    :: snd_movfin_grid(mxg), snd_movco_grid(mxg), snd_overlap_grid(mxg)
  integer :: byte_movfin_grid(mxg), byte_movco_grid(mxg), byte_overlap_grid(mxg)
  integer :: count_movfin_grid(mxg), count_movco_grid(mxg), count_overlap_grid(mxg)

  common /timer_inc/whoami, num_of_grid, ops_grid, ops_level, time_grid, time_addw, &
                  op_addw,  op_collc, op_euler, &
                  op_solver, op_step, op_movfin, op_movco, &
                  op_overlap, op_overhead, &
                  op_gradface, op_fluxvis, &
                  dat_movfin, dat_movco, dat_overlap, &
                  snd_movfin, snd_movco, snd_overlap, &
                  snd_movfin_grid, snd_movco_grid, snd_overlap_grid, &
                  byte_movfin_grid, byte_movco_grid, byte_overlap_grid, &
                  count_movfin_grid, count_movco_grid, count_overlap_grid
