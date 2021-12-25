## make targets
- run1/run2: runs the implementations for task1/task2 (run.c)
- validate: validates task1&task2 against MPI_Allreduce (validate.c)
- benchmark: runs benchmarking for task1&task2 (benchmark.c)
  - produces test1.log&test2.log
- benchmark-hydra: runs benchmarks on the hydra system (benchmark.c)
  - produces NODESxPPN.log files
- profile: creates performance profiles of task1&task2 implementations using [mpiP](https://github.com/LLNL/mpiP) (profile.c)
  - needs some setup to work
