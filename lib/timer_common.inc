  integer :: level, onlevel(0:12)
  common/timer_private/ level, onlevel
  !$omp threadprivate(/timer_private/)
