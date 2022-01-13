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
    fprintf(fptr, "size;blocksize;custom;allreduce;reduce+bcast\n");
  }

  static double input[10000000] = {0.0};
  static double out0[10000000] = {0.0};
  static double out1[10000000] = {0.0};
  static double out2[10000000] = {0.0};

  const int elements[] = {0, 1, 2, 8, 15, 21, 25, 87, 150, 212, 250, 875, 1500, 2125, 2500, 8750, 15000, 21250, 25000, 87500, 150000, 250000, 1000000, 2500000, 10000000};
  const int blocksizes[] = {250, 500, 750, 2500, 16000};

  // vector size
  for (int size_i = 0; size_i < 23; size_i++) {
    int size = elements[size_i];
    for (int blocksize_i = 0; blocksize_i < 5; blocksize_i++) {
      blocksize = blocksizes[blocksize_i];

      for (int v = 1; v <= size; v++) {
        input[v-1] = ((v + rank) % 10) + 1;
        input[v-1] = ((v + rank) % 10) + 1;
      }

      double times[3] = {0.0, 0.0, 0.0};

      for (int i = 0; i < 12; i++) {
        double measured[3] = {0.0, 0.0, 0.0};
        MPI_Barrier(MPI_COMM_WORLD);
        measured[0] = MPI_Wtime();
        AllReduce(input, out0, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        measured[0] = MPI_Wtime() - measured[0];

        if (blocksize_i == 0) {
          MPI_Barrier(MPI_COMM_WORLD);
          measured[1] = MPI_Wtime();
          MPI_Allreduce(input, out1, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          measured[1] = MPI_Wtime() - measured[1];

          MPI_Barrier(MPI_COMM_WORLD);
          measured[2] = MPI_Wtime();
          MPI_Reduce(input, out2, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
          MPI_Bcast(out2, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
          measured[2] = MPI_Wtime() - measured[2];

        }
        
        MPI_Allreduce(MPI_IN_PLACE, measured, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        times[0] += measured[0];
        times[1] += measured[1];
        times[2] += measured[2];
      }

      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == 0) fprintf(fptr, "%d;%d;%lf;%lf;%lf\n", size, blocksize, times[0]/12.0*1000, times[1]/12.0*1000, times[2]/12.0*1000);
    }
  }

  if (rank == 0) printf("Benchmarks completed!\n");
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  run_benchmark(argv[1]);
  MPI_Finalize();
}