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

void nodeinfo(int rank, int comm_size, int* node, int* parent, int* left, int* right, int* layer, int* treeheight, int* dual) {
  // treeheight for half the nodes
  *treeheight = floor(log2(ceil(comm_size/2.0)));
  int treesize = pow(2,*treeheight+1) - 1;
  // no. of processes - size of n-1 tree
  int overshoot_l = (ceil(comm_size/2.0) - (pow(2,*treeheight) - 1)) * 2;
  int overshoot_r = (floor(comm_size/2.0) - (pow(2,*treeheight) - 1)) * 2;

  int offset = rank >= comm_size/2.0 ? ceil(comm_size/2.0) : 0;
  int overshoot = offset == 0 ? overshoot_l : overshoot_r;

  // skip any virtual nodes
  *node = do_skips(rank - offset, overshoot);
  *layer = get_layer(*node);

  *parent = get_parent(*node);
  *left =  *layer == 0 ? -1 : *node - (1 << *layer-1);
  *right = *layer == 0 ? -1 : *node + (1 << *layer-1);

  if (get_layer(*parent) > *treeheight) *parent = -1;
  if (*left % 2 == 0 && *left >= overshoot) *left = -1;
  if (*right % 2 == 0 && *right >= overshoot) *right = -1;

  if (*parent != -1) *parent = undo_skips(*parent, overshoot);
  if (*left != -1) *left = undo_skips(*left, overshoot);
  if (*right != -1) *right = undo_skips(*right, overshoot);

  if (offset != 0) {
    *node += treesize;
    if (*parent != -1) *parent += ceil(comm_size/2.0);
    if (*left != -1)   *left   += ceil(comm_size/2.0);
    if (*right != -1)  *right  += ceil(comm_size/2.0);
  }

  if (*left >= comm_size) *left = -1;
  if (*right >= comm_size) *right = -1;

  if (*parent == -1) {
    int root = (1 << *treeheight) - 1;
    *dual = undo_skips(root, offset == 0 ? overshoot_r : overshoot_l) + (rank >= comm_size/2.0 ? 0 : ceil(comm_size/2.0));
  }

  // printf("rank %d, node %d, overshoot %d, parent %d, left %d, right %d, layer %d, dual %d\n", rank, *node, overshoot, *parent, *left, *right, *layer, *dual);
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

  int typesize;
  MPI_Type_size(datatype, &typesize);
  void* value = malloc(typesize*count);
  memcpy(value, sendbuf, typesize*count);

  // trivial case
  if (comm_size == 1) {
    memcpy(recvbuf, sendbuf, typesize*count);
    return;
  }

  int rank, node, parent, left, right, layer, dual, treesize, treeheight;
  MPI_Comm_rank(comm, &rank);
  nodeinfo(rank, comm_size, &node, &parent, &left, &right, &layer, &treeheight, &dual);

  void* bufup = malloc(typesize*count);

  int round = -1;
  int rounds = 2*(2*treeheight + blockcount);

  while (++round < rounds) {
    // is node on the right side of parent to send up/recv down (even round => left child, odd round => right child)
    bool matchesEvenOdd = parent < rank && round % 2 == 1 || parent > rank && round % 2 == 0;
    int child = round % 2 == 0 ? left : right;

    // --- SENDING UP---
    bool isAncestorReceivingUpYet = parent != -1 && round/2 >= layer - 2;
    bool hasUnsentUp = sentup != blockcount && (recvdup >= 2 || layer == 0);

    bool shouldSendUp = matchesEvenOdd && hasUnsentUp && isAncestorReceivingUpYet;


    // --- RECEIVING UP---
    bool isReceivingUpYet = round/2 >= layer - 1;
    bool hasUnreceivedUp = recvdup/2 != blockcount;
    bool hasDescendent = child != -1;

    bool shouldReceiveUp = isReceivingUpYet && hasUnreceivedUp && hasDescendent;

    // --- SENDING DOWN ---
    bool hasUnsentDown = sentdown/2 != blockcount && recvddown + swapped >= 1;
    bool hasSenddownStarted = round/2 - treeheight >= treeheight - layer;

    bool shouldSendDown = hasDescendent && hasUnsentDown && hasSenddownStarted;

    // --- RECEIVING DOWN---
    bool hasUnreceivedDown = recvddown != blockcount;
    bool isAncestorSendingDownYet = parent != -1 && round/2 >= (2*treeheight-layer-1);

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

    if ((ops & 0b1001) == 0b1001) {
      // printf("[%d] send up (%d) %d -> %d\n", round, sentup, rank, parent);
      // printf("[%d] recv dn (%d) %d -> %d\n", round, recvddown, parent, rank);
      MPI_Sendrecv(send_up_from, send_up_size, datatype, parent, 0,
                  recv_down_into,  recv_down_size, datatype, parent, 0,
                  comm, MPI_STATUS_IGNORE);

      ops &= 0b0110;
    }

    if ((ops & 0b0110) == 0b0110) {
      // printf("[%d] send dn (%d) %d -> %d\n", round, sentdown/2, rank, child);
      // printf("[%d] recv up (%d) %d -> %d\n", round, recvdup/2, child, rank);
      MPI_Sendrecv(send_down_from, send_down_size, datatype, child, 0,
                  recv_up_into,  recv_up_size, datatype, child, 0,
                  comm, MPI_STATUS_IGNORE);
      ops &= 0b1001;
    }

    if ((ops & 0b1000) == 0b1000) {
      // printf("[%d] send up (%d) %d -> %d\n", round, sentup, rank, parent);
      MPI_Send(send_up_from, send_up_size, datatype, parent, 0, comm);
    }

    if ((ops & 0b0100) == 0b0100) {
      // printf("[%d] recv up (%d) %d -> %d\n", round, recvdup/2, child, rank);
      MPI_Recv(recv_up_into, recv_up_size, datatype, child, 0, comm, MPI_STATUS_IGNORE);
    }

    if ((ops & 0b0010) == 0b0010) {
      // printf("[%d] send dn (%d) %d -> %d\n", round, sentdown/2, rank, child);
      MPI_Send(send_down_from, send_down_size, datatype, child, 0, comm);
    }

    if ((ops & 0b0001) == 0b0001) {
      // printf("[%d] recv dn (%d) %d -> %d\n", round, recvddown, parent, rank);
      MPI_Recv(recv_down_into, recv_down_size, datatype, parent, 0, comm, MPI_STATUS_IGNORE);
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
    if (hasUnsentDown && hasSenddownStarted) sentdown++;

    // --- DUAL ROOT SWAP ---
    bool isRoot = parent == -1;
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
    }
  }
}