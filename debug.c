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
#include <limits.h>                                  /*  UINT_MAX  */
#include <stdbool.h>
#include <math.h>
#include <signal.h>
#include <time.h>

#include "wmalloc.h"
#include "plsa-defn.h"
#include "debug.h"

void handler_sigfpe () {
  signal (SIGFPE, handler_sigfpe); /* reset signal */
  fprintf (stderr, "-->  I have received a SIGFPE!\n");

  return;
}

/*  Used by input.c  */
void debugCheckCo (INFO *info) {
  unsigned int j = 0;
  unsigned int cos_count = 0;
  unsigned int curr_j;

  /*  Check flags in co-occurrence table  */
  for (unsigned int i = 0; i < info -> m; i++) {
    cos_count = GET_COS_POSITION (i, 0);
    curr_j = 0;
    for (unsigned int pos_j = 1; pos_j <= cos_count; pos_j++) {
      j = GET_COS_POSITION (i, pos_j);
      while (curr_j < j) {
        fprintf (stderr, "X");
        curr_j++;
      }
      fprintf (stderr, "O");
    }

    while (curr_j < info -> n) {
      fprintf (stderr, "X");
      curr_j++;
    }

    fprintf (stderr, "\n");
  }

  fprintf (stderr, "-----\n");

  return;
}

/*  Used by output.c  */
void checkCoProb (INFO *info) {
  unsigned int num_clusters = info -> num_clusters;
  unsigned int i = 0;  /*  Index into w1  */
  unsigned int j = 0;  /*  Index into w2  */
  unsigned int k = 0;  /*  Index into clusters  */
  PROBNODE temp;
  PROBNODE tempsum = 0.0;
  unsigned int nonprob = 0;

  for (i = 0; i < info -> m; i++) {
    for (j = 0; j < info -> n; j++) {
      temp = GET_PROBZ_W1W2_CURR (0, i, j);
      for (k = 1; k < num_clusters; k++) {
        /*  temp stores logarithms  */
        logSumsInline (temp, (GET_PROBZ_W1W2_CURR (k, i, j)));
      }
      if (temp > 0) {
        nonprob++;
      }
      tempsum += DOEXP (temp);
    }
  }

  if (info -> verbose) {
    fprintf (stderr, "**\t%u : %f\n", nonprob, tempsum);
  }

  return;
}


void printAllProbsPrev (INFO *info) {
  unsigned int i = 0;
  unsigned int count = 0;

  if ((info -> world_size > 1) && (info -> world_id == MAINPROC)) {
    return;
  }

  count = 0;
  fprintf (stderr, "\n[%u] P[%u] ===== p(w1|z) =====\n", info -> world_id, info -> iter);
  for (i = 0; i < (info -> num_clusters * info -> m); i++) {
    fprintf (stderr, "%f", info -> probw1_z_prev[i]);
    count++;
    if ((count % info -> m) == 0) {
      fprintf (stderr, "\n");
    }
    else {
      fprintf (stderr, "\t");
    }
  }

  count = 0;
  fprintf (stderr, "\n[%u] P[%u] ===== p(w2|z) =====\n", info -> world_id, info -> iter);
  for (i = 0; i < (info -> num_clusters * info -> n); i++) {
    fprintf (stderr, "%f", info -> probw2_z_prev[i]);
    count++;
    if ((count % info -> n) == 0) {
      fprintf (stderr, "\n");
    }
    else {
      fprintf (stderr, "\t");
    }
  }

  fprintf (stderr, "\n[%u] P[%u] ===== p(z) =====\n", info -> world_id, info -> iter);
  for (i = 0; i < info -> num_clusters; i++) {
    fprintf (stderr, "%f\t", info -> probz_prev[i]);
    count++;
  }
  fprintf (stderr, "\n==========\n");

  return;
}


void printAllProbsCurr (INFO *info) {
  unsigned int i = 0;
  unsigned int count = 0;

  if ((info -> world_size > 1) && (info -> world_id == MAINPROC)) {
    return;
  }

  count = 0;
  fprintf (stderr, "\n[%u] C[%u] ===== p(w1|z) =====\n", info -> world_id, info -> iter);
  for (i = 0; i < (info -> num_clusters * info -> m); i++) {
    fprintf (stderr, "%f", info -> probw1_z_curr[i]);
    count++;
    if ((count % info -> m) == 0) {
      fprintf (stderr, "\n");
    }
    else {
      fprintf (stderr, "\t");
    }
  }

  count = 0;
  fprintf (stderr, "\n[%u] C[%u] ===== p(w2|z) =====\n", info -> world_id, info -> iter);
  for (i = 0; i < (info -> num_clusters * info -> n); i++) {
    fprintf (stderr, "%f", info -> probw2_z_curr[i]);
    count++;
    if ((count % info -> n) == 0) {
      fprintf (stderr, "\n");
    }
    else {
      fprintf (stderr, "\t");
    }
  }

  fprintf (stderr, "\n[%u] C[%u] ===== p(z) =====\n", info -> world_id, info -> iter);
  for (i = 0; i < info -> num_clusters; i++) {
    fprintf (stderr, "%f\t", info -> probz_curr[i]);
    count++;
  }
  fprintf (stderr, "\n==========\n");

  return;
}


void printJointProb (INFO *info) {
  unsigned int i = 0;
  unsigned int j = 0;

  for (i = 0; i < info -> m; i++) {
    for (j = 0; j < info -> n; j++) {
      fprintf (stderr, "[%f]\t", info -> prob_w1w2[i * info -> n + j]);
    }
    fprintf (stderr, "\n");
  }

  return;
}


