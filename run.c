#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

// forward declaration
int blocksize;
void AllReduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  double out[100000] = {0.0};
  double in[100000] = {[0 ... 99999] = 1.0};

  blocksize = 550;
  AllReduce(in, out, 100000, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("%d:\t%lf, %lf, %lf, %lf, ...\n", rank, out[0], out[1], out[2], out[3]);
  MPI_Finalize();
}