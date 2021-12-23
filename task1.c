#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))

int blocksize = 1;

int depth(int val) {
  // __builtin_ctz: index of the least significant set bit
  return __builtin_ctz(~val);
}

int parent (int node) {
  int D = depth(node);
  int p = node >> (D+2) << (D+2);
  return p | (1 << (D+1)) - 1;
} 

void nodeinfo(int node, int comm_size, int* ancestor, int* left, int* right, int* layer) {
  *layer = depth(node);
  *ancestor = parent(node);
  int treeheight = floor(log2(comm_size));
  while (*ancestor >= comm_size && depth(*ancestor) <= treeheight) *ancestor = parent(*ancestor);
  if (*ancestor >= comm_size) *ancestor = -1;
  *left =  *layer == 0 ? -1 : node - (1 << *layer-1);
  *right = *layer == 0 ? -1 : node + (1 << *layer-1);

  int d = *layer-2;
  while (*right >= comm_size) {
    if (d == -1) {
      *right = -1;
      break;
    }
    *right = *right - (1 << d);
    d--;
  }
}

int _Reduce (const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  // --- SETUP ---
  int sentblocks = 0;
  int recvdblocks = 0;
  int blockcount = ceil((double) count/blocksize);

  int comm_size;
  MPI_Comm_size(comm, &comm_size);
  int treeheight = floor(log2(comm_size));

  int node;
  MPI_Comm_rank(comm, &node);
  int ancestor, left, right, layer;
  nodeinfo(node, comm_size, &ancestor, &left, &right, &layer);

  int typesize;
  MPI_Type_size(datatype, &typesize);
  void* value = malloc(typesize*count);
  memcpy(value, sendbuf, typesize*count);

  int round = -1;
  int rounds = 2*(treeheight + blockcount);

  while (++round < rounds) {
    // --- SENDING ---
    // is node on the right side of parent to send up (even round => left child sends, odd round => right child sends)
    bool matchesEvenOdd = ancestor < node && round % 2 == 1 || ancestor > node && round % 2 == 0;
    // is there data to be sent
    bool hasBlocksToSend = sentblocks != blockcount && (recvdblocks >= 2 || layer == 0);
    // in case nodes between this and ancestor are missing, make sure the ancestor is ready to receive
    bool isAncestorReceiving = ancestor != -1 && round/2 >= depth(ancestor) - 1;

    bool shouldSend = matchesEvenOdd && hasBlocksToSend && isAncestorReceiving;

    // --- RECEIVING ---
    int descendent = round % 2 == 0 ? left : right;
    // has node already started recieving
    bool isReceiving = round/2 >= layer - 1;
    // are there unrecieved blocks left for this node
    bool hasUnreceived = recvdblocks/2 != blockcount;
    // does descendent exist
    bool hasDescendent = descendent != -1;

    bool shouldReceive = isReceiving && hasUnreceived && hasDescendent;

    int chunk_to_send = min(count - (sentblocks)*blocksize, blocksize);
    int chunk_to_recv = min(count - (recvdblocks/2)*blocksize, blocksize);

    // --- MPI CALLS ---
    if (shouldSend && shouldReceive) {
      MPI_Sendrecv(value + typesize*blocksize*sentblocks,       chunk_to_send, datatype, ancestor, 0,
                  recvbuf + typesize*blocksize*(recvdblocks/2), chunk_to_recv, datatype, descendent, 0,
                  comm, MPI_STATUS_IGNORE);
    } else if (shouldSend) {
      MPI_Send(value + typesize*blocksize*sentblocks, chunk_to_send, datatype, ancestor, 0, comm);
    } else if (shouldReceive) {
      MPI_Recv(recvbuf + typesize*blocksize*(recvdblocks/2), chunk_to_recv, datatype, descendent, 0, comm, MPI_STATUS_IGNORE);
    }

    // --- LOCAL REDUCTION --- 
    if (shouldReceive) {
      int offset = typesize*blocksize*(recvdblocks/2);
      // local reduction
      if (round % 2 == 0) {
        MPI_Reduce_local(recvbuf + offset, value + offset, chunk_to_recv, datatype, op);
      } else {
        MPI_Reduce_local(value + offset, recvbuf + offset, chunk_to_recv, datatype, op);
        memcpy(value+offset, recvbuf+offset, chunk_to_recv*typesize);
      }

    }

    // increment recvdblocks even if one of the ancestors does not exists, to keep steps in sync
    if (isReceiving && hasUnreceived) recvdblocks++;
    if (shouldSend) {
      sentblocks++;
    }
  }

  memcpy(recvbuf, value, typesize*count);

  return (1 << treeheight) - 1;
}

void _Broadcast (void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
  int sentblocks = 0;
  int recvdblocks = 0;
  int blockcount = ceil((double)count/blocksize);

  int comm_size;
  MPI_Comm_size(comm, &comm_size);
  int treeheight = floor(log2(comm_size));

  int node;
  MPI_Comm_rank(comm, &node);
  int ancestor, left, right, layer;
  nodeinfo(node, comm_size, &ancestor, &left, &right, &layer);

  int typesize;
  MPI_Type_size(datatype, &typesize);

  int round = -1;
  int rounds = 2*(treeheight + blockcount);

  while (++round < rounds) {
    // --- SENDING ---
    int descendent = round % 2 == 0 ? left : right;
    // does descendent exist
    bool hasDescendent = descendent != -1;
    // is there data left to be send
    bool hasBlocksLeft = sentblocks/2 != blockcount && (recvdblocks >= 1 || layer == treeheight);

    bool shouldSend =  hasDescendent && hasBlocksLeft;

    // --- RECEIVING ---
    // is node on the right side of parent to send up (even round => left child receives, odd round => right child receives)
    bool matchesEvenOdd = ancestor < node && round % 2 == 1 || ancestor > node && round % 2 == 0;
    // are there unreceived blocks left for this node
    bool hasUnreceived = recvdblocks != blockcount;
    // has the ancestor sent a block this round
    bool hasSendingAncestor = ancestor != -1 && round/2 >= (treeheight-depth(ancestor));

    bool shouldReceive = matchesEvenOdd && hasUnreceived && hasSendingAncestor;

    int chunk_to_send = min(count - (sentblocks/2)*blocksize, blocksize);
    int chunk_to_recv = min(count - (recvdblocks)*blocksize, blocksize);

    // --- MPI CALLS ---
    if (shouldSend && shouldReceive) {
      MPI_Sendrecv(buffer + typesize*blocksize*(sentblocks/2), chunk_to_send, datatype, descendent, 0,
                  buffer + typesize*blocksize*recvdblocks,     chunk_to_recv, datatype, ancestor, 0,
                  comm, MPI_STATUS_IGNORE);
    } else if (shouldSend) {
      MPI_Send(buffer + typesize*blocksize*(sentblocks/2), chunk_to_send, datatype, descendent, 0, comm);
    } else if (shouldReceive) {
      MPI_Recv(buffer + typesize*blocksize*recvdblocks, chunk_to_recv, datatype, ancestor, 0, comm, MPI_STATUS_IGNORE);
    }

    if (hasBlocksLeft) sentblocks++;
    if (shouldReceive) recvdblocks++;
  }
}

void AllReduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  int root = _Reduce(sendbuf, recvbuf, count, datatype, op, comm);
  _Broadcast(recvbuf, count, datatype, root, comm);
}


void run_test() {
  double out0[100000] = {0.0};
  double in0[100000] = {0.0};
  double out1[100000] = {0.0};
  double in1[100000] = {0.0};

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  const int sizes[] = {1, 10, 25, 36, 99, 103, 455, 1023, 42311, 74589};

  // vector size
  for (int size_i = 0; size_i < 10; size_i++) {
    // blocksize (global var)
    int size = sizes[size_i];
    printf("%d\n", size);
    for (blocksize = 1; blocksize < 100000; blocksize *= 10) {
      for (int processes = 2; processes < 10; processes += 1) {
        MPI_Comm comm;
        MPI_Comm_split(MPI_COMM_WORLD, rank < processes ? 0 : MPI_UNDEFINED, 0, &comm);
        if (comm != MPI_COMM_NULL) {

          for (int v = 1; v <= size; v++) {
            out0[v-1] = ((v + rank) % 10) + 1;
            out1[v-1] = ((v + rank) % 10) + 1;
          }

          AllReduce(out0, in0, size, MPI_DOUBLE, MPI_SUM, comm);
          MPI_Allreduce(out1, in1, size, MPI_DOUBLE, MPI_SUM, comm);

          for (int i = 0; i < size; i++) {
            if (in0[i] != in1[i]) {
              printf("%d,%d: %lf != %lf\n", processes, i, in0[i], in1[i]);
            }
          }

          AllReduce(out0, in0, size, MPI_DOUBLE, MPI_PROD, comm);
          MPI_Allreduce(out1, in1, size, MPI_DOUBLE, MPI_PROD, comm);

          for (int i = 0; i < size; i++) {
            if (in0[i] != in1[i]) {
              printf("%d,%d: %lf != %lf\n", processes, i, in0[i], in1[i]);
            }
          }
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);

  // Get the rank of the calling process
  // int world_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // double out[] = {[0 ... 10000] = 2.124};
  // double in[10000] = {0.0};

  // AllReduce(out, in, 10000, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // // MPI_Allreduce(out, in, 10000, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // printf("result: (%lf, %lf, %lf, %lf, %lf)\n", in[0], in[240], in[3333], in[7489], in[9999]);

  run_test();
  // Finalize: Any resources allocated for MPI can be freed
  MPI_Finalize();
}