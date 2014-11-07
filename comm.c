/*
**  Probabilistic latent semantic analysis (PLSA, multiprocessor version)
**  Copyright (C) 2009-2010  by Raymond Wan (r.wan@aist.go.jp)
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>

#include "config.h"
#if HAVE_MPI
#include <mpi.h>
#endif

#include "wmalloc.h"
#include "plsa-defn.h"
#include "comm.h"

#if HAVE_MPI
static void recvProbsFromOthers (INFO *info) {
  unsigned int k = 0;  /*  Cluster position on the main processor  */
  int result = 0;
  unsigned int tag = 0;
  unsigned int owner = 0;
  MPI_Status *status = wmalloc (sizeof (MPI_Status));

  if (info -> world_size == 1) {
    return;
  }

  for (k = 0; k < info -> num_clusters; k++) {
    owner = BLOCK_OWNER (k, info -> world_size, info -> num_clusters);
    if (owner != 0) {
      MSG_RECV_STATUS (info -> world_id, owner, info -> iter, TAG_PROBW1_Z, k);
      tag = MSG_TAG (info -> iter, TAG_PROBW1_Z, k);
      result = MPI_Recv (info -> probw1_z_curr + (k * info -> m), (info -> m), MPI_TYPE, owner, tag, MPI_COMM_WORLD, status);

      MSG_RECV_STATUS (info -> world_id, owner, info -> iter, TAG_PROBW2_Z, k);
      tag = MSG_TAG (info -> iter, TAG_PROBW2_Z, k);
      result = MPI_Recv (info -> probw2_z_curr + (k * info -> n), (info -> n), MPI_TYPE, owner, tag, MPI_COMM_WORLD, status);

      MSG_RECV_STATUS (info -> world_id, owner, info -> iter, TAG_PROBZ, k);
      tag = MSG_TAG (info -> iter, TAG_PROBZ, k);
      result = MPI_Recv (&info -> probz_curr[k], 1, MPI_TYPE, owner, tag, MPI_COMM_WORLD, status);
    }
  }

  wfree (status);

  return;
}

static void sendProbsToMain (INFO *info) {
  unsigned int k = 0;  /*  Cluster position on the main processor  */
  unsigned int p = 0;  /*  Cluster position on the non-main processors  */
  int result = 0;
  unsigned int tag = 0;

  for (p = 0; p < info -> block_size; p++) {
    k = info -> block_start + p;

    /*  Send p(i|z)  */
    MSG_SEND_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBW1_Z, k);
    tag = MSG_TAG (info -> iter, TAG_PROBW1_Z, k);
    result = MPI_Send (info -> probw1_z_curr + (p * info -> m), (info -> m), MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD);

    /*  Send p(j|z)  */
    MSG_SEND_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBW2_Z, k);
    tag = MSG_TAG (info -> iter, TAG_PROBW2_Z, k);
    result = MPI_Send (info -> probw2_z_curr + (p * info -> n), (info -> n), MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD);

    /*  Send p(z)  */
    MSG_SEND_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBZ, k);
    tag = MSG_TAG (info -> iter, TAG_PROBZ, k);
    result = MPI_Send (&info -> probz_curr[p], 1, MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD);
  }

  return;
}


static void sendProbsToOthers (INFO *info) {
  unsigned int k = 0;  /*  Cluster position on the main processor  */
  int result = 0;
  unsigned int tag = 0;
  unsigned int owner = 0;

  for (k = 0; k < info -> num_clusters; k++) {
    owner = BLOCK_OWNER (k, info -> world_size, info -> num_clusters);

    /*  If I am not the destination...  */
    if (owner != info -> world_id) {
      /*  Send p(i|z)  */
      MSG_SEND_STATUS (info -> world_id, owner, info -> iter, TAG_PROBW1_Z, k);
      tag = MSG_TAG (info -> iter, TAG_PROBW1_Z, k);
      result = MPI_Send (info -> probw1_z_curr + (k * info -> m), (info -> m), MPI_TYPE, owner, tag, MPI_COMM_WORLD);

      /*  Send p(j|z)  */
      MSG_SEND_STATUS (info -> world_id, owner, info -> iter, TAG_PROBW2_Z, k);
      tag = MSG_TAG (info -> iter, TAG_PROBW2_Z, k);
      result = MPI_Send (info -> probw2_z_curr + (k * info -> n), (info -> n), MPI_TYPE, owner, tag, MPI_COMM_WORLD);

      MSG_SEND_STATUS (info -> world_id, owner, info -> iter, TAG_PROBZ, k);
      tag = MSG_TAG (info -> iter, TAG_PROBZ, k);
      result = MPI_Send (&info -> probz_curr[k], 1, MPI_TYPE, owner, tag, MPI_COMM_WORLD);
    }
  }

  return;
}


static void recvProbsFromMain (INFO *info) {
  unsigned int k = 0;  /*  Cluster position on the main processor  */
  unsigned int p = 0;  /*  Cluster position on the non-main processors  */
  int result = 0;
  unsigned int tag = 0;
  unsigned int owner = 0;
  MPI_Status *status = wmalloc (sizeof (MPI_Status));

  for (k = 0; k < info -> num_clusters; k++) {
    owner = BLOCK_OWNER (k, info -> world_size, info -> num_clusters);

    /*  If the owner is me  */
    if (owner == info -> world_id) {
      /*  Send p(i|z)  */
      MSG_RECV_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBW1_Z, k);
      tag = MSG_TAG (info -> iter, TAG_PROBW1_Z, k);
      result = MPI_Recv (info -> probw1_z_curr + (p * info -> m), (info -> m), MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD, status);

      /*  Send p(j|z)  */
      MSG_RECV_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBW2_Z, k);
      tag = MSG_TAG (info -> iter, TAG_PROBW2_Z, k);
      result = MPI_Recv (info -> probw2_z_curr + (p * info -> n), (info -> n), MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD, status);

      /*  Send p(z)  */
      MSG_RECV_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBZ, k);
      tag = MSG_TAG (info -> iter, TAG_PROBZ, k);
      result = MPI_Recv (&info -> probz_curr[p], 1, MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD, status);
      p++;
    }
  }

  wfree (status);

  return;
}
#else
static void recvProbsFromOthers (INFO *info) {
  return;
}

static void sendProbsToMain (INFO *info) {
  return;
}


static void sendProbsToOthers (INFO *info) {
  return;
}


static void recvProbsFromMain (INFO *info) {
  return;
}
#endif

/*!  MAINPROC sends the initialized (*current*) p(i|z), p(j|z), and p(z) to all other processes  */
void distributeProbs (INFO *info) {
  time_t start;
  time_t end;

  time (&start);
  if (info -> world_size == 1) {
    return;
  }

  if (info -> world_id == MAINPROC) {
    sendProbsToOthers (info);
  }
  else {
    recvProbsFromMain (info);
  }
  time (&end);
  info -> distributeProbs_time += difftime (end, start);

  return;
}


/*!  All other processes send the *current* p(i|z), p(j|z), and p(z) to MAINPROC  */
void gatherProbs (INFO *info) {
  time_t start;
  time_t end;

  time (&start);
  if (info -> world_size == 1) {
    return;
  }

  if (info -> world_id == MAINPROC) {
    recvProbsFromOthers (info);
  }
  else {
    sendProbsToMain (info);
  }
  time (&end);
  info -> gatherProbs_time += difftime (end, start);

  return;
}

