#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define min(X, Y) (((X) < (Y)) ? (X) : (Y))

int blocksize = 1;

int depth(int val) {
  return __builtin_ctz(~val);
}

int parent (int node) {
  int D = depth(node);
  int p = node >> (D+2) << (D+2);
  return p | (1 << (D+1)) - 1;
} 

void nodeinfo(int node, int comm_size, int *treeheight, int *treesize, int* ancestor, int* left, int* right, int* layer, int* dual) {
  *treeheight = floor(log2(ceil(comm_size/2.0)));
  *treesize = pow(2,*treeheight+1) - 1;
  int delta = node >= *treesize ? *treesize : 0;
  node -= delta;
  *layer = depth(node);
  *ancestor = *layer == *treeheight ? -1 : parent(node);

  *left =  *layer == 0 ? -1 : node - (1 << *layer-1);
  *right = *layer == 0 ? -1 : node + (1 << *layer-1);

  if (*ancestor != -1) *ancestor += delta;
  if (*left != -1) *left += delta;
  if (*right != -1) *right += delta;

  while (*ancestor >= comm_size && depth(*ancestor-delta) < *treeheight) *ancestor = parent(*ancestor-delta)+delta;
  if (*ancestor >= comm_size) *ancestor = -1;

  int d = *layer-2;
  while (*right >= comm_size) {
    if (d == -1) {
      *right = -1;
      break;
    }
    *right = (*right-delta) - (1 << d) + delta;
    d--;
  }

  // not strictly necessary to set to -1 for non-roots
  *dual = -1;

  if (*ancestor == -1) {
    // root node on right side
    if (delta != 0) *dual = *treesize/2;
    else {
      *dual = node + *treesize;
      d = *layer-1;
      while (*dual >= comm_size) {
        if (d == -1) {
          *dual = -1;
          break;
        }
        *dual = (*dual-delta) - (1 << d) + delta;
        d--;
      }
    }
  }
}

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

  // printf("node %d, ancestor %d, left %d, right %d, layer %d, dual %d\n", node, ancestor, left, right, layer, dual);
  int typesize;
  MPI_Type_size(datatype, &typesize);
  void* value = malloc(typesize*count);
  memcpy(value, sendbuf, typesize*count);

  void* bufup = malloc(typesize*count);

  int round = -1;
  int rounds = 2*(2*treeheight + blockcount);

  while (++round < rounds) {
    // --- SENDING UP---
    // is node on the right side of parent to send up (even round => left child sends, odd round => right child sends)
    bool matchesEvenOdd = ancestor < node && round % 2 == 1 || ancestor > node && round % 2 == 0;
    // is there data to be sent
    bool hasBlocksToSend = sentup != blockcount && (recvdup >= 2 || layer == 0);
    // in case nodes between this and ancestor are missing, make sure the ancestor is ready to receive
    bool isAncestorReceiving = ancestor != -1 && round/2 >= depth(ancestor) - 1;

    bool shouldSendUp = matchesEvenOdd && hasBlocksToSend && isAncestorReceiving;


    // --- RECEIVING UP---
    int descendent = round % 2 == 0 ? left : right;
    // has node already started recieving
    bool isReceiving = round/2 >= layer - 1;
    // are there unreceived blocks left for this node
    bool hasUnreceivedUp = recvdup/2 != blockcount;
    // does descendent exist
    bool hasDescendent = descendent != -1;

    bool shouldReceiveUp = isReceiving && hasUnreceivedUp && hasDescendent;

    // --- SENDING DOWN ---
    // is there data left to be sent
    bool hasBroadcastStarted = round/2 >= treeheight;
    bool hasBlocksToSendDown = sentdown/2 != blockcount && (recvddown >= 1 || ancestor == -1);
    bool isDescendentReceivingDown = round/2 >= (2*treeheight-layer);

    bool shouldSendDown = hasDescendent && hasBlocksToSendDown && hasBroadcastStarted && isDescendentReceivingDown;

    // --- RECEIVING DOWN---
    // are there unreceived blocks left for this node
    bool hasUnreceivedDown = recvddown != blockcount;
    // has the ancestor sent a block this round
    bool hasSendingAncestor = ancestor != -1 && round/2 >= (2*treeheight-depth(ancestor%treesize));

    bool shouldReceiveDown = matchesEvenOdd && hasUnreceivedDown && hasSendingAncestor && hasBroadcastStarted;

    int send_up_size   = min(count - (sentup)*blocksize, blocksize);
    int recv_up_size   = min(count - (recvdup/2)*blocksize, blocksize);
    int send_down_size = min(count - (sentdown/2)*blocksize, blocksize);
    int recv_down_size = min(count - (recvddown)*blocksize, blocksize);
    int swap_size      = min(count - swapped*blocksize, blocksize);

    void* send_up_from     = value + typesize*blocksize*sentup;
    void* recv_up_into     = bufup + typesize*blocksize*(recvdup/2);
    void* send_down_from   = recvbuf + typesize*blocksize*(sentdown/2);
    void* recv_down_into   = recvbuf + typesize*blocksize*recvddown;
    void* swap_from        = value + typesize*blocksize*swapped;
    void* swap_into        = recvbuf + typesize*blocksize*swapped; 

    int ops = shouldSendUp << 3 | shouldReceiveUp << 2 | shouldSendDown << 1 | shouldReceiveDown;

    // --- MPI CALLS ---
    // there has to be a better way of pairing up send/recv...
    if ((ops & 0b1100) == 0b1100) {
      // printf("[%d] send up (%d:%lf) %d -> %d\n", round, sentup, ((double*)send_up_from)[0], node, ancestor);
      // printf("[%d] recv up (%d:?) %d -> %d\n", round, (recvdup/2), node, descendent);
      MPI_Sendrecv(send_up_from, send_up_size, datatype, ancestor, 0,
                  recv_up_into,  recv_up_size, datatype, descendent, 0,
                  comm, MPI_STATUS_IGNORE);
      ops &= 0b0011;
    }

    if ((ops & 0b0011) == 0b0011) {
      // printf("[%d] send down (%d:%lf) %d -> %d\n", round, sentdown/2, ((double*)send_down_from)[0], node, descendent);
      // printf("[%d] recv down (%d:?) %d -> %d\n", round, recvddown, ancestor, node);
      MPI_Sendrecv(send_down_from, send_down_size, datatype, descendent, 0,
                  recv_down_into,  recv_down_size, datatype, ancestor, 0,
                  comm, MPI_STATUS_IGNORE);

      ops &= 0b1100;
    }

    if ((ops & 0b1001) == 0b1001) {
      // printf("[%d] send up (%d:%lf) %d -> %d\n", round, sentup, ((double*)send_up_from)[0], node, ancestor);
      // printf("[%d] recv down (%d:?) %d -> %d\n", round, (recvdup/2), node, ancestor);
      MPI_Sendrecv(send_up_from, send_up_size, datatype, ancestor, 0,
                  recv_down_into,  recv_down_size, datatype, ancestor, 0,
                  comm, MPI_STATUS_IGNORE);

      ops &= 0b0110;
    }

    if ((ops & 0b0110) == 0b0110) {
      // printf("[%d] send down (%d:%lf) %d -> %d\n", round, sentdown/2, ((double*)send_down_from)[0], node, descendent);
      // printf("[%d] recv up (%d:?) %d -> %d\n", round, recvddown, descendent, node);
      MPI_Sendrecv(send_down_from, send_down_size, datatype, descendent, 0,
                  recv_up_into,  recv_up_size, datatype, descendent, 0,
                  comm, MPI_STATUS_IGNORE);
      ops &= 0b1001;
    }

    if ((ops & 0b1000) == 0b1000) {
      // printf("[%d] send up (%d:%lf) %d -> %d\n", round, sentup, ((double*)send_up_from)[0], node, ancestor);
      MPI_Send(send_up_from, send_up_size, datatype, ancestor, 0, comm);
    }

    if ((ops & 0b0100) == 0b0100) {
      // printf("[%d] recv up (%d:?) %d -> %d\n", round, (recvdup/2), descendent, node);
      MPI_Recv(recv_up_into, recv_up_size, datatype, descendent, 0, comm, MPI_STATUS_IGNORE);
    }

    if ((ops & 0b0010) == 0b0010) {
      // printf("[%d] send down (%d:%lf) %d -> %d\n", round, sentdown/2, ((double*)send_down_from)[0], node, descendent);
      MPI_Send(send_down_from, send_down_size, datatype, descendent, 0, comm);
    }

    if ((ops & 0b0001) == 0b0001) {
      // printf("[%d] recv down (%d:?) %d -> %d\n", round, recvddown, ancestor, node);
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

    // increment recvdblocks even if one of the ancestors does not exists, to keep steps in sync
    if (isReceiving && hasUnreceivedUp) recvdup++;
    if (shouldSendUp) sentup++;
    if (shouldReceiveDown) recvddown++;
    if (hasBlocksToSendDown && hasBroadcastStarted && isDescendentReceivingDown) sentdown++;

    // --- DUAL ROOT SWAP ---
    bool isRoot = ancestor == -1;
    bool hasNewReduct = round % 2 == 1 && recvdup != 0 && swapped != blockcount;
    bool hasSwappingStarted = (round+1)/2 >= treeheight;
    bool hasDual = dual != -1;
    bool shouldSwap = isRoot && hasNewReduct && hasSwappingStarted && hasDual;

    if (shouldSwap) {
      // printf("[%d] swap %d <-> %d (%d:%lf)\n", round, node, dual, swapped, ((double*)swap_from)[0]);
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

int main(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);

  double in[] = { 1.0, 2.0, 3.0 };
  double out[] = { 0.0, 0.0, 0.0 };
  AllReduce(in, out, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  int node;
  MPI_Comm_rank(MPI_COMM_WORLD, &node);

  printf("%d: %lf, %lf, %lf\n", node, out[0], out[1], out[2]);

  // Get the rank of the calling process
  // int world_rank;
  // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // double out[] = {[0 ... 10000] = 2.124};
  // double in[10000] = {0.0};

  // AllReduce(out, in, 10000, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // // MPI_Allreduce(out, in, 10000, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // printf("result: (%lf, %lf, %lf, %lf, %lf)\n", in[0], in[240], in[3333], in[7489], in[9999]);

  // run_test();
  // Finalize: Any resources allocated for MPI can be freed
  MPI_Finalize();
}