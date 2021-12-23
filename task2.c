#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))

int blocksize = 1;

int get_layer(int val) {
  return __builtin_ctz(~val);
}

int parent (int node) {
  int D = get_layer(node);
  int p = node >> (D+2) << (D+2);
  return p | (1 << (D+1)) - 1;
} 

// constructing the tree 
void nodeinfo(int node, int comm_size, int *treeheight, int *treesize, int* ancestor, int* left, int* right, int* layer, int* dual) {
  *treeheight = floor(log2(ceil(comm_size/2.0)));
  *treesize = pow(2,*treeheight+1) - 1;

  int offset = node >= *treesize ? *treesize : 0;
  node -= offset;
  *layer = get_layer(node);
  *ancestor = *layer == *treeheight ? -1 : parent(node);

  *left =  *layer == 0 ? -1 : node - (1 << *layer-1);
  *right = *layer == 0 ? -1 : node + (1 << *layer-1);

  if (*ancestor != -1) *ancestor += offset;
  if (*left != -1) *left += offset;
  if (*right != -1) *right += offset;

  while (*ancestor >= comm_size && get_layer(*ancestor-offset) < *treeheight) *ancestor = parent(*ancestor-offset)+offset;
  if (*ancestor >= comm_size) *ancestor = -1;

  int d = *layer-2;
  while (*right >= comm_size) {
    if (d == -1) {
      *right = -1;
      break;
    }
    *right = (*right-offset) - (1 << d) + offset;
    d--;
  }

  if (*ancestor == -1) {
    // root node on right side
    if (offset != 0) *dual = *treesize/2;
    // root node on left side
    else {
      *dual = node + *treesize;
      d = *layer-1;
      while (*dual >= comm_size) {
        if (d == -1) {
          *dual = -1;
          break;
        }
        *dual = (*dual-offset) - (1 << d) + offset;
        d--;
      }
    }
  }
}

// allreduce implementation
void AllReduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
  int sentup = 0;
  int recvdup = 0;
  int sentdown = 0;
  int recvddown = 0;
  int swapped = 0;

  int blockcount = ceil((double) count/blocksize);

  int comm_size;
  MPI_Comm_size(comm, &comm_size);

  int node, ancestor, left, right, layer, dual, treesize, treeheight;
  MPI_Comm_rank(comm, &node);
  nodeinfo(node, comm_size, &treeheight, &treesize, &ancestor, &left, &right, &layer, &dual);
  int ancestorget_layer = get_layer(ancestor%treesize);

  int typesize;
  MPI_Type_size(datatype, &typesize);
  void* value = malloc(typesize*count);
  memcpy(value, sendbuf, typesize*count);

  void* bufup = malloc(typesize*count);

  int round = -1;
  int rounds = 2*(2*treeheight + blockcount);

  while (++round < rounds) {
    // is node on the right side of parent to send up/recv down (even round => left child, odd round => right child)
    bool matchesEvenOdd = ancestor < node && round % 2 == 1 || ancestor > node && round % 2 == 0;
    int descendent = round % 2 == 0 ? left : right;

    // --- SENDING UP---
    bool isAncestorReceivingUpYet = ancestor != -1 && round/2 >= ancestorget_layer - 1;
    bool hasUnsentUp = sentup != blockcount && (recvdup >= 2 || layer == 0);

    bool shouldSendUp = matchesEvenOdd && hasUnsentUp && isAncestorReceivingUpYet;


    // --- RECEIVING UP---
    bool isReceivingUpYet = round/2 >= layer - 1;
    bool hasUnreceivedUp = recvdup/2 != blockcount;
    bool hasDescendent = descendent != -1;

    bool shouldReceiveUp = isReceivingUpYet && hasUnreceivedUp && hasDescendent;

    // --- SENDING DOWN ---
    bool hasUnsentDown = sentdown/2 != blockcount && recvddown + swapped >= 1;

    bool shouldSendDown = hasDescendent && hasUnsentDown;

    // --- RECEIVING DOWN---
    bool hasUnreceivedDown = recvddown != blockcount;
    bool isAncestorSendingDownYet = ancestor != -1 && round/2 >= (2*treeheight-ancestorget_layer);

    bool shouldReceiveDown = matchesEvenOdd && hasUnreceivedDown && isAncestorSendingDownYet;

    int send_up_size   = min(count - (sentup)*blocksize,     blocksize);
    int recv_up_size   = min(count - (recvdup/2)*blocksize,  blocksize);
    int send_down_size = min(count - (sentdown/2)*blocksize, blocksize);
    int recv_down_size = min(count - (recvddown)*blocksize,  blocksize);
    int swap_size      = min(count - swapped*blocksize,      blocksize);

    void* send_up_from     = value   + typesize*blocksize*sentup;
    void* recv_up_into     = bufup   + typesize*blocksize*(recvdup/2);
    void* send_down_from   = recvbuf + typesize*blocksize*(sentdown/2);
    void* recv_down_into   = recvbuf + typesize*blocksize*recvddown;
    void* swap_from        = value   + typesize*blocksize*swapped;
    void* swap_into        = recvbuf + typesize*blocksize*swapped; 


    // --- MPI CALLS ---
    // there has to be a better way of pairing up send/recv...
    int ops = shouldSendUp << 3 | shouldReceiveUp << 2 | shouldSendDown << 1 | shouldReceiveDown;

    if ((ops & 0b1100) == 0b1100) {
      // printf("[%d] send up (%d) %d -> %d\n", round, sentup, node, ancestor);
      // printf("[%d] recv up (%d) %d -> %d\n", round, recvdup/2, descendent, node);
      MPI_Sendrecv(send_up_from, send_up_size, datatype, ancestor, 0,
                  recv_up_into,  recv_up_size, datatype, descendent, 0,
                  comm, MPI_STATUS_IGNORE);
      ops &= 0b0011;
    }

    if ((ops & 0b0011) == 0b0011) {
      // printf("[%d] send dn (%d) %d -> %d\n", round, sentdown/2, node, descendent);
      // printf("[%d] recv dn (%d) %d -> %d\n", round, recvddown, ancestor, node);
      MPI_Sendrecv(send_down_from, send_down_size, datatype, descendent, 0,
                  recv_down_into,  recv_down_size, datatype, ancestor, 0,
                  comm, MPI_STATUS_IGNORE);

      ops &= 0b1100;
    }

    if ((ops & 0b1001) == 0b1001) {
      // printf("[%d] send up (%d) %d -> %d\n", round, sentup, node, ancestor);
      // printf("[%d] recv dn (%d) %d -> %d\n", round, recvddown, ancestor, node);
      MPI_Sendrecv(send_up_from, send_up_size, datatype, ancestor, 0,
                  recv_down_into,  recv_down_size, datatype, ancestor, 0,
                  comm, MPI_STATUS_IGNORE);

      ops &= 0b0110;
    }

    if ((ops & 0b0110) == 0b0110) {
      // printf("[%d] send dn (%d) %d -> %d\n", round, sentdown/2, node, descendent);
      // printf("[%d] recv up (%d) %d -> %d\n", round, recvdup/2, descendent, node);
      MPI_Sendrecv(send_down_from, send_down_size, datatype, descendent, 0,
                  recv_up_into,  recv_up_size, datatype, descendent, 0,
                  comm, MPI_STATUS_IGNORE);
      ops &= 0b1001;
    }

    if ((ops & 0b1000) == 0b1000) {
      // printf("[%d] send up (%d) %d -> %d\n", round, sentup, node, ancestor);
      MPI_Send(send_up_from, send_up_size, datatype, ancestor, 0, comm);
    }

    if ((ops & 0b0100) == 0b0100) {
      // printf("[%d] recv up (%d) %d -> %d\n", round, recvdup/2, descendent, node);
      MPI_Recv(recv_up_into, recv_up_size, datatype, descendent, 0, comm, MPI_STATUS_IGNORE);
    }

    if ((ops & 0b0010) == 0b0010) {
      // printf("[%d] send dn (%d) %d -> %d\n", round, sentdown/2, node, descendent);
      MPI_Send(send_down_from, send_down_size, datatype, descendent, 0, comm);
    }

    if ((ops & 0b0001) == 0b0001) {
      // printf("[%d] recv dn (%d) %d -> %d\n", round, recvddown, ancestor, node);
      MPI_Recv(recv_down_into, recv_down_size, datatype, ancestor, 0, comm, MPI_STATUS_IGNORE);
    }

    // --- LOCAL REDUCTION --- 
    if (shouldReceiveUp) {
      void* reduce_into   = value + typesize*blocksize*(recvdup/2);
      // local reduction
      if (round % 2 == 0) {
        MPI_Reduce_local(recv_up_into, reduce_into, recv_up_size, datatype, op);
      } else {
        MPI_Reduce_local(reduce_into, recv_up_into, recv_up_size, datatype, op);
        memcpy(reduce_into, recv_up_into, recv_up_size*typesize);
      }
    }

    // increment recvdup/sentdown even if other process does not exists, to keep steps in sync
    if (isReceivingUpYet && hasUnreceivedUp) recvdup++;
    if (shouldSendUp) sentup++;
    if (shouldReceiveDown) recvddown++;
    if (hasUnsentDown) sentdown++;

    // --- DUAL ROOT SWAP ---
    bool isRoot = ancestor == -1;
    bool hasNewReduct = round % 2 == 1 && recvdup != 0 && swapped != blockcount;
    bool hasSwappingStarted = (round+1)/2 >= treeheight;
    bool hasDual = dual != -1;
    bool shouldSwap = isRoot && hasNewReduct && hasDual && hasSwappingStarted;

    if (shouldSwap) {
      MPI_Sendrecv(swap_from, swap_size, datatype, dual, 0,
                  swap_into, swap_size, datatype, dual, 0,
                  comm, MPI_STATUS_IGNORE);

      if (node < treesize) {
        MPI_Reduce_local(swap_from, swap_into, swap_size, datatype, op);
      } else {
        MPI_Reduce_local(swap_into, swap_from, swap_size, datatype, op);
        memcpy(swap_into, swap_from, swap_size*typesize);
      }

      swapped++;
    } else if (isRoot && hasNewReduct && hasSwappingStarted) {
      memcpy(swap_into, swap_from, swap_size*typesize);
      swapped++;
    }
  }
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

  if (rank == 0) printf("Tests completed\n");
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
