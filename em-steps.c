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
#include <string.h>
#include <limits.h>
#include <math.h>  /*  log10 function  */
#include <stdbool.h>
#include <time.h>

#include "PLSA_MP_Config.h"
#if HAVE_MPI
#include <mpi.h>
#endif

#include "wmalloc.h"
#include "plsa-defn.h"
#include "em-steps.h"


void swapPrevCurr (INFO *info) {
  PROBNODE *probw1_z_temp;
  PROBNODE *probw2_z_temp;
  PROBNODE *probz_temp;
  time_t start;
  time_t end;

  time (&start);
  probw1_z_temp = info -> probw1_z_prev;
  probw2_z_temp = info -> probw2_z_prev;
  probz_temp = info -> probz_prev;

  info -> probw1_z_prev = info -> probw1_z_curr;
  info -> probw2_z_prev = info -> probw2_z_curr;
  info -> probz_prev = info -> probz_curr;

  info -> probw1_z_curr = probw1_z_temp;
  info -> probw2_z_curr = probw2_z_temp;
  info -> probz_curr = probz_temp;

  time (&end);
  info -> swapPrevCurr_time += difftime (end, start);

  return;
}


void initEM (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  unsigned int i;  /*  Index into w1  */
  unsigned int j;  /*  Index into w2  */
  unsigned int k;  /*  Index into clusters  */
  PROBNODE sum;
  time_t start;
  time_t end;

  time (&start);
  PROGRESS_MSG ("Begin initialization...");

  /*  Assign probabilities to probz  */
  sum = 0.0;
  for (k = 0; k < num_clusters; k++) {
    GET_PROBZ_CURR (k) = RANDOM_FLOAT;
    sum += GET_PROBZ_CURR (k);
  }
  for (k = 0; k < num_clusters; k++) {
    GET_PROBZ_CURR (k) = DOLOG (GET_PROBZ_CURR (k) / sum);
  }

  /*  Assign probabilities to probw1_z  */
  for (i = 0; i < (info -> num_clusters * info -> m); i++) {
    info -> probw1_z_curr[i] = RANDOM_FLOAT;
  }
  for (k = 0; k < num_clusters; k++) {
    sum = 0.0;
    for (i = 0; i < info -> m; i++) {
      sum += GET_PROBW1_Z_CURR (k, i);
    }
    for (i = 0; i < info -> m; i++) {
      GET_PROBW1_Z_CURR (k, i) = DOLOG (GET_PROBW1_Z_CURR (k, i) / sum);
    }
  }

  /*  Assign probabilities to probw2_z  */
  for (i = 0; i < (info -> num_clusters * info -> n); i++) {
    info -> probw2_z_curr[i] = RANDOM_FLOAT;
  }
  for (k = 0; k < num_clusters; k++) {
    sum= 0.0;
    for (j = 0; j < info -> n; j++) {
      sum += GET_PROBW2_Z_CURR (k, j);
    }
    for (j = 0; j < info -> n; j++) {
      GET_PROBW2_Z_CURR (k, j) = DOLOG (GET_PROBW2_Z_CURR (k, j) / sum);
    }
  }

  PROGRESS_MSG ("Initialization complete...");

  time (&end);
  info -> initEM_time += difftime (end, start);

  return;
}


void applyEMStep (INFO *info) {
  unsigned int i = 0;  /*  Index into w1  */
  unsigned int j = 0;  /*  Index into w2  */
  signed int k = 0;  /*  Index into clusters, local to this processor  */
  register unsigned int pos_j;  /*  Actual position in the cooccurrence array  */
  register unsigned int cos_count;  /*  Number of cooccurrences in each row  */
  PROBNODE temp;
  PROBNODE cos;
  bool *flag_z = NULL;
  bool **flag_w1_z = NULL;
  bool **flag_w2_z = NULL;

  time_t start;
  time_t end;

  time (&start);

  /*******************************************************/
  /*  Initialize flags that indicate whether the cell is so far untouched   */

  flag_z = wmalloc (info -> block_size * sizeof (bool));
  for (k = 0; k < info -> block_size; k++) {
    flag_z[k] = false;
  }

  flag_w1_z = wmalloc (info -> block_size * sizeof (bool*));
  for (k = 0; k < info -> block_size; k++) {
    flag_w1_z[k] = wmalloc (info -> m * sizeof (bool));
    for (i = 0; i < info -> m; i++) {
      flag_w1_z[k][i] = false;
    }
  }

  flag_w2_z = wmalloc (info -> block_size * sizeof (bool*));
  for (k = 0; k < info -> block_size; k++) {
    flag_w2_z[k] = wmalloc (info -> n * sizeof (bool));
    for (j = 0; j < info -> n; j++) {
      flag_w2_z[k][j] = false;
    }
  }


  /*******************************************************/
  /*  Update probabilities */

#if HAVE_OPENMP
#pragma omp parallel for private(i,cos_count,pos_j,j,cos,temp)
#endif
  for (k = 0; k < info -> block_size; k++) {
    for (i = 0; i < info -> m; i++) {
      cos_count = GET_COS_POSITION (i, 0);
      for (pos_j = 1; pos_j <= cos_count; pos_j++) {
        j = GET_COS_POSITION (i, pos_j);
        cos = GET_COS (i, pos_j);
        temp = (GET_PROBZ_W1W2_PREV (k, i, j)) - GET_PROB_W1W2 (i, j);

        /*  probz  */
        if (flag_z[k]) {
          logSumsInline (GET_PROBZ_CURR (k), cos + temp);
        }
        else {
          GET_PROBZ_CURR (k) = cos + temp;
          flag_z[k] = true;
        }

        /*  probw1_z  */
        if (flag_w1_z[k][i]) {
          logSumsInline (GET_PROBW1_Z_CURR (k, i), cos + temp);
        }
        else {
          GET_PROBW1_Z_CURR (k, i) = cos + temp;
          flag_w1_z[k][i] = true;
        }

        /*  probw2_z  */
        if (flag_w2_z[k][j]) {
          logSumsInline (GET_PROBW2_Z_CURR (k, j), cos + temp);
        }
        else {
          GET_PROBW2_Z_CURR (k, j) = cos + temp;
          flag_w2_z[k][j] = true;
        }
      }
    }
  }

  /*******************************************************/
  /*  Uninitialize flags  */

  wfree (flag_z);
  for (k = 0; k < info -> block_size; k++) {
    wfree (flag_w1_z[k]);
    wfree (flag_w2_z[k]);
  }
  wfree (flag_w1_z);
  wfree (flag_w2_z);

  time (&end);
  info -> applyEMStep_time += difftime (end, start);

  return;
}


PROBNODE calculateML (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  signed int i;  /*  Index into w1  */
  signed int j;  /*  Index into w2  */
  signed int k;  /*  Index into clusters  */
  unsigned int pos_j;  /*  Actual position in the cooccurrence array  */
  unsigned int cos_count;  /*  Number of cooccurrences in each row  */
  PROBNODE total = 0.0;
  PROBNODE temp;
  time_t start;
  time_t end;

  time (&start);

#if HAVE_OPENMP
#pragma omp parallel for private(cos_count,pos_j,j,temp,k) reduction(+:total)
#endif
  for (i = 0; i < info -> m; i++) {
    cos_count = GET_COS_POSITION (i, 0);
    for (pos_j = 1; pos_j <= cos_count; pos_j++) {
      j = GET_COS_POSITION (i, pos_j);

      /*  Initialize with cluster 0  */
      temp = GET_PROBZ_W1W2_CURR (0,i,j);
      /*  Log-likelihood for the co-occurrence of two words  */
      for (k = 1; k < num_clusters; k++) {
        /*  temp stores log values  */
        logSumsInline (temp, (GET_PROBZ_W1W2_CURR (k,i,j)));
      }

      /*  Log-likelihood across all examples  */
      total += (temp * DOEXP (GET_COS (i, pos_j)));
    }
  }

  time (&end);
  info -> calculateML_time += difftime (end, start);

  return (total);
}


void calculateProbW1W2 (INFO *info) {
  signed int i;  /*  Index into w1  */
  unsigned int j;  /*  Index into w2  */
  unsigned int k;  /*  Index into clusters  */
  PROBNODE temp = 0.0;
  PROBNODE *temp_prob_w1w2 = NULL;
  int result = 0;
  unsigned int tag = 0;
  unsigned int owner = 0;
#if HAVE_MPI
  MPI_Status *status = wmalloc (sizeof (MPI_Status));
#endif
  time_t start;
  time_t end;

  time (&start);

#if HAVE_OPENMP
#pragma omp parallel for private(j,temp,k)
#endif
  for (i = 0; i < info -> m; i++) {
    for (j = 0; j < info -> n; j++) {
      temp = GET_PROBZ_W1W2_CURR (0,i,j);
      for (k = 1; k < info -> block_size; k++) {
        /*  temp stores logarithms  */
        logSumsInline (temp, (GET_PROBZ_W1W2_CURR (k,i,j)));
      }
      GET_PROB_W1W2(i,j) = temp;
    }
  }

#if HAVE_MPI
  if (info -> world_id == MAINPROC) {
    temp_prob_w1w2 = wmalloc (info -> m * info -> n * sizeof (PROBNODE));

    /*  Receive each temporary matrix and copy it in  */
    for (owner = 1; owner < info -> world_size; owner++) {
      MSG_RECV_STATUS (info -> world_id, owner, info -> iter, TAG_PROBW1W2, 0);
      tag = MSG_TAG (info -> iter, TAG_PROBW1W2, 0);
      result = MPI_Recv (temp_prob_w1w2, (info -> m * info -> n), MPI_TYPE, owner, tag, MPI_COMM_WORLD, status);

#if HAVE_OPENMP
#pragma omp parallel for private(j)
#endif
      for (i = 0; i < info -> m; i++) {
        for (j = 0; j < info -> n; j++) {
          logSumsInline (GET_PROB_W1W2(i,j), temp_prob_w1w2[i * info -> n + j]);
        }
      }
    }

    wfree (temp_prob_w1w2);
  }
  else {
    MSG_SEND_STATUS (info -> world_id, MAINPROC, info -> iter, TAG_PROBW1W2, 0);
    tag = MSG_TAG (info -> iter, TAG_PROBW1W2, 0);
    result = MPI_Send (info -> prob_w1w2, (info -> m * info -> n), MPI_TYPE, MAINPROC, tag, MPI_COMM_WORLD);
  }
#endif

  time (&end);
  info -> calculateProbW1W2_time += difftime (end, start);

  return;
}


/*!  Normalize probabilities  */
void normalizeProbs (INFO *info) {
  unsigned int i;  /*  Index into w1  */
  unsigned int j;  /*  Index into w2  */
  signed int k;  /*  Index into clusters  */
  PROBNODE sum;
  PROBNODE norm;
  time_t start;
  time_t end;

  time (&start);

#if HAVE_OPENMP
#pragma omp parallel for private(norm,i,j)
#endif
  for (k = 0; k < info -> num_clusters; k++) {
    norm = GET_PROBZ_CURR (k);

    /*  probw1_z  */
    for (i = 0; i < info -> m; i++) {
      GET_PROBW1_Z_CURR (k, i) = GET_PROBW1_Z_CURR (k, i) - norm;
    }

    /*  probw2_z  */
    for (j = 0; j < info -> n; j++) {
      GET_PROBW2_Z_CURR (k, j) = GET_PROBW2_Z_CURR (k, j) - norm;
    }
  }

  /*  probz  */
  sum = GET_PROBZ_CURR (0);
  for (k = 1; k < info -> num_clusters; k++) {
    logSumsInline (sum, GET_PROBZ_CURR (k));
  }
  for (k = 0; k < info -> num_clusters; k++) {
    GET_PROBZ_CURR (k) = GET_PROBZ_CURR (k) - sum;
  }

  time (&end);
  info -> normalizeProbs_time += difftime (end, start);

  return;
}
