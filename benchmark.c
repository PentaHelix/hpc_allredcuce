#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

// forward declaration
int blocksize;
void AllReduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

void run_validation() {
  double out0[100000] = {0.0};
  double in0[100000] = {0.0};
  double out1[100000] = {0.0};
  double in1[100000] = {0.0};

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int elements[] = {0, 1, 2, 8, 15, 21, 25, 87, 150, 212, 250, 875, 1500, 2125, 2500, 8750, 15000, 21250, 25000, 87500, 150000, 212500, 250000, 875000, 1500000, 2125000, 2500000, 4597152, 6694304, 8388608};
  // const int elements[] = {0, 1, 2, 8, 15, 21, 25, 87, 150, 212, 250, 875, 1500};
  const int blocksizes[] = {1, 31, 68, 100, 200, 500, 700, 1000, 2500};

  // vector size
  for (int size_i = 0; size_i < 30; size_i++) {
  // for (int size_i = 0; size_i < 12; size_i++) {
    int size = elements[size_i];
    for (int blocksize_i = 0; blocksize_i < 9; blocksize_i++) {
      blocksize = blocksizes[blocksize_i];
      for (int v = 1; v <= size; v++) {
        out0[v-1] = ((v + rank) % 10) + 1;
        out1[v-1] = ((v + rank) % 10) + 1;
      }

      double times[2] = {0.0, 0.0};

      for (int i = 0; i < 16; i++) {
        double measured[2];
        MPI_Barrier(MPI_COMM_WORLD);
        measured[0] = MPI_Wtime();
        AllReduce(out0, in0, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        measured[0] = MPI_Wtime() - measured[0];

        MPI_Barrier(MPI_COMM_WORLD);
        measured[1] = MPI_Wtime();
        MPI_Allreduce(out1, in1, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        measured[1] = MPI_Wtime() - measured[1];

        MPI_Allreduce(MPI_IN_PLACE, measured, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, measured, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        times[0] += measured[0];
        times[1] += measured[1];
      }
      if (rank == 0) printf("[%d, %d, %lf, %lf],\n", size, blocksize, times[0]/16*1000, times[1]/16*1000);
    }
  }

  if (rank == 0) printf("Benchmarks completed!\n");
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  run_validation();
  MPI_Finalize();
}