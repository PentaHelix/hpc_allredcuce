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

  const int elements[] = {2500, 10, 25, 36, 99, 103, 455, 1023, 42311, 74589};

  // vector size
  for (int size_i = 0; size_i < 10; size_i++) {
    int size = elements[size_i];
    for (blocksize = 1; blocksize <= 100000; blocksize *= 10) {
      for (int v = 1; v <= size; v++) {
        out0[v-1] = ((v + rank) % 10) + 1;
        out1[v-1] = ((v + rank) % 10) + 1;
      }

      AllReduce(out0, in0, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(out1, in1, size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (int i = 0; i < size; i++) {
        if (in0[i] != in1[i]) {
          printf("%d: %lf != %lf\n", i, in0[i], in1[i]);
        }
      }
    }
  }

  if (rank == 0) printf("Tests completed!\n");
}


int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  run_validation();
  MPI_Finalize();
}