#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))
#define max(X, Y) (((X) > (Y)) ? (X) : (Y))

int blocksize = 1;

int get_layer(int val) {
  return __builtin_ctz(~val);
}

int get_parent (int node) {
  int D = get_layer(node);
  int p = node >> (D+2) << (D+2);
  return p | (1 << (D+1)) - 1;
} 

int do_skips (int rank, int overshoot) {
  int skips = 0;
  while (rank + skips >= overshoot + skips*2) skips++;
  return rank + skips;
}

int undo_skips (int node, int overshoot) {
  int skipped = max(0, (node - overshoot + 1)/2);
  return node - skipped;
}

void nodeinfo(int rank, int comm_size, int* node, int* parent, int* left, int* right, int* layer, int* treeheight) {
  *treeheight = floor(log2(comm_size));
  // no. of processes - size of n-1 tree
  int overshoot = (comm_size - (pow(2,*treeheight) - 1)) * 2;

  // skip any virtual nodes
  *node = do_skips(rank, overshoot);
  *layer = get_layer(*node);

  *parent = get_parent(*node);
  *left =  *layer == 0 ? -1 : *node - (1 << *layer-1);
  *right = *layer == 0 ? -1 : *node + (1 << *layer-1);

  if (*left % 2 == 0 && *left >= overshoot) *left = -1;
  if (*right % 2 == 0 && *right >= overshoot) *right = -1;

  *parent = undo_skips(*parent, overshoot);
  *left = undo_skips(*left, overshoot);
  *right = undo_skips(*right, overshoot);

  if (*parent >= comm_size) *parent = -1;

  // printf("rank %d, node %d, parent %d, left %d, right %d, layer %d\n", rank, *node, *parent, *left, *right, *layer);
}

void _Reduce (const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  // --- SETUP ---
  int sentblocks = 0;
  int recvdblocks = 0;
  int blockcount = ceil((double) count/blocksize);

  int comm_size;
  MPI_Comm_size(comm, &comm_size);

  int rank;
  MPI_Comm_rank(comm, &rank);
  int node, parent, left, right, layer, treeheight;
  nodeinfo(rank, comm_size, &node, &parent, &left, &right, &layer, &treeheight);

  int typesize;
  MPI_Type_size(datatype, &typesize);
  void* value = malloc(typesize*count);
  memcpy(value, sendbuf, typesize*count);

  int round = -1;
  int rounds = 2*(treeheight + blockcount);


  while (++round < rounds) {
    // --- SENDING ---
    // is node on the right side of parent to send up (even round => left child sends, odd round => right child sends)
    bool matchesEvenOdd = parent < rank && round % 2 == 1 || parent > rank && round % 2 == 0;
    // is there data to be sent
    bool hasBlocksToSend = sentblocks != blockcount && (recvdblocks >= 2 || layer == 0);
    bool isLayerSending = round/2 >= layer;
    bool hasParent = parent != -1;

    // if (rank == 7) printf("[%d] %d %d %d %d\n", round, matchesEvenOdd, hasBlocksToSend, isLayerSending, hasParent);

    bool shouldSend = matchesEvenOdd && hasBlocksToSend && isLayerSending && hasParent;

    // --- RECEIVING ---
    int child = round % 2 == 0 ? left : right;
    // has node already started recieving
    bool isReceiving = round/2 >= layer - 1;
    // are there unrecieved blocks left for this node
    bool hasUnreceived = recvdblocks/2 != blockcount;
    // does child exist
    bool hasDescendent = child != -1;

    bool shouldReceive = isReceiving && hasUnreceived && hasDescendent;

    int chunk_to_send = min(count - (sentblocks)*blocksize, blocksize);
    int chunk_to_recv = min(count - (recvdblocks/2)*blocksize, blocksize);

    // --- MPI CALLS ---
    if (shouldSend && shouldReceive) {
      // printf("[%d] send %d -> %d\n", round, rank, parent);
      // printf("[%d] recv %d -> %d\n", round, child, rank);
      MPI_Sendrecv(value + typesize*blocksize*sentblocks,       chunk_to_send, datatype, parent, 0,
                  recvbuf + typesize*blocksize*(recvdblocks/2), chunk_to_recv, datatype, child, 0,
                  comm, MPI_STATUS_IGNORE);
    } else if (shouldSend) {
      // printf("[%d] send %d -> %d\n", round, rank, parent);
      MPI_Send(value + typesize*blocksize*sentblocks, chunk_to_send, datatype, parent, 0, comm);
    } else if (shouldReceive) {
      // printf("[%d] recv %d -> %d\n", round, child, rank);
      MPI_Recv(recvbuf + typesize*blocksize*(recvdblocks/2), chunk_to_recv, datatype, child, 0, comm, MPI_STATUS_IGNORE);
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

    // increment recvdblocks even if one of the childs does not exists, to keep steps in sync
    if (isReceiving && hasUnreceived) recvdblocks++;
    if (shouldSend) {
      sentblocks++;
    }
  }

  memcpy(recvbuf, value, typesize*count);
}

void _Broadcast (void *buffer, int count, MPI_Datatype datatype, MPI_Comm comm) {
  int sentblocks = 0;
  int recvdblocks = 0;
  int blockcount = ceil((double)count/blocksize);

  int comm_size;
  MPI_Comm_size(comm, &comm_size);

  int rank;
  MPI_Comm_rank(comm, &rank);
  int node, parent, left, right, layer, treeheight;
  nodeinfo(rank, comm_size, &node, &parent, &left, &right, &layer, &treeheight);

  int typesize;
  MPI_Type_size(datatype, &typesize);

  int round = -1;
  int rounds = 2*(treeheight + blockcount);

  while (++round < rounds) {
    // --- SENDING ---
    int child = round % 2 == 0 ? left : right;
    bool hasDescendent = child != -1;
    bool hasBlocksLeft = sentblocks/2 != blockcount && (recvdblocks >= 1 || layer == treeheight);

    bool shouldSend =  hasDescendent && hasBlocksLeft;

    // --- RECEIVING ---
    bool matchesEvenOdd = parent < rank && round % 2 == 1 || parent > rank && round % 2 == 0;
    // are there unreceived blocks left for this node
    bool hasUnreceived = recvdblocks != blockcount;
    // has the parent sent a block this round
    bool hasParent = parent != -1 && round/2 >= (treeheight-layer-1);

    bool shouldReceive = matchesEvenOdd && hasUnreceived && hasParent;

    int chunk_to_send = min(count - (sentblocks/2)*blocksize, blocksize);
    int chunk_to_recv = min(count - (recvdblocks)*blocksize, blocksize);

    // --- MPI CALLS ---
    if (shouldSend && shouldReceive) {
      // printf("[%d] send %d -> %d\n", round, rank, child);
      // printf("[%d] recv %d -> %d\n", round, parent, rank);
      MPI_Sendrecv(buffer + typesize*blocksize*(sentblocks/2), chunk_to_send, datatype, child, 0,
                  buffer + typesize*blocksize*recvdblocks,     chunk_to_recv, datatype, parent, 0,
                  comm, MPI_STATUS_IGNORE);
    } else if (shouldSend) {
      // printf("[%d] send %d -> %d\n", round, rank, child);
      MPI_Send(buffer + typesize*blocksize*(sentblocks/2), chunk_to_send, datatype, child, 0, comm);
    } else if (shouldReceive) {
      // printf("[%d] recv %d -> %d\n", round, parent, rank);
      MPI_Recv(buffer + typesize*blocksize*recvdblocks, chunk_to_recv, datatype, parent, 0, comm, MPI_STATUS_IGNORE);
    }

    if (hasBlocksLeft) sentblocks++;
    if (shouldReceive) recvdblocks++;
  }
}

void AllReduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  _Reduce(sendbuf, recvbuf, count, datatype, op, comm);
  _Broadcast(recvbuf, count, datatype, comm);
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
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  double out[] = {1.0, 1.0, 1.0, 1.0, 1.0};
  double in[] = {0.0, 0.0, 0.0, 0.0, 0.0};

  AllReduce(out, in, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  printf("result: (%lf, %lf)\n", in[0], in[1]);

  // run_test();
  // Finalize: Any resources allocated for MPI can be freed
  MPI_Finalize();
}