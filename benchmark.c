#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

// forward declaration
int blocksize;
void AllReduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

void run_benchmark(char* filename) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *fptr;
  if (rank == 0) {
    fptr = fopen(filename, "w");
    fprintf(fptr, "size;blocksize;custom;mpi\n");
  }

  static double out0[150000] = {0.0};
  static double in0[150000] = {0.0};
  static double out1[150000] = {0.0};
  static double in1[150000] = {0.0};

  const int elements[] = {0, 1, 2, 8, 15, 21, 25, 87, 150, 212, 250, 875, 1500, 2125, 2500, 8750, 15000, 21250, 25000, 87500, 150000};
  const int blocksizes[] = {8, 250, 500, 750, 2500};

  // vector size
  for (int size_i = 0; size_i < 21; size_i++) {
    int size = elements[size_i];
    for (int blocksize_i = 0; blocksize_i < 5; blocksize_i++) {
      blocksize = blocksizes[blocksize_i];

      for (int v = 1; v <= size; v++) {
        out0[v-1] = ((v + rank) % 10) + 1;
        out1[v-1] = ((v + rank) % 10) + 1;
      }

      double times[2] = {0.0, 0.0};

      for (int i = 0; i < 32; i++) {
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

        times[0] += measured[0];
        times[1] += measured[1];
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0) fprintf(fptr, "%d;%d;%lf;%lf\n", size, blocksize, times[0]/32.0*1000, times[1]/32.0*1000);
    }
  }

  if (rank == 0) printf("Benchmarks completed!\n");
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  run_benchmark(argv[1]);
  MPI_Finalize();
}