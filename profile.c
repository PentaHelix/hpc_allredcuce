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
  double in[100000] = {0.0};

  blocksize = 550;
  AllReduce(in, out, 100000, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Finalize();
}