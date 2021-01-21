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
#include <math.h>  /*  fabs  */
#include <float.h>  /*  DBL_EPSILON  */
#include <time.h>
#include <signal.h>

#include "PLSA_MP_Config.h"
#if HAVE_MPI
#include <mpi.h>
#endif

#include "wmalloc.h"
#include "plsa-defn.h"
#include "em-steps.h"
#include "input.h"
#include "output.h"
#include "parameters.h"
#include "debug.h"
#include "comm.h"
#include "run.h"


INFO *initialize () {
  INFO *info = wmalloc (sizeof (INFO));

  time (&(info -> program_start));
  info -> run_time = 0;
  info -> readCO_time = 0;
  info -> initEM_time = 0;
  info -> calculateProbW1W2_time = 0;
  info -> calculateML_time = 0;
  info -> swapPrevCurr_time = 0;
  info -> applyEMStep_time = 0;
  info -> gatherProbs_time = 0;
  info -> normalizeProbs_time = 0;
  info -> distributeProbs_time = 0;
  info -> printCoProbs_time = 0;

  /*  MPI may or may not be in use; assume it is not and set defaults  */
  info -> world_id = MAINPROC;
  info -> world_size = 1;
  info -> threads = 0;

  /*  Set a handler for floating point exceptions  */
  info -> sigfpe_count = 0;
  signal (SIGFPE, handler_sigfpe);

  return info;
}


void uninitialize (INFO *info) {
  double total_time = 0;

  wfree (info -> cos);
  wfree (info -> probw1_z_curr);
  wfree (info -> probw2_z_curr);
  wfree (info -> probz_curr);
  wfree (info -> probw1_z_prev);
  wfree (info -> probw2_z_prev);
  wfree (info -> probz_prev);
  wfree (info -> base_fn);
  wfree (info -> co_fn);
  wfree (info -> row_ids);
  wfree (info -> column_ids);

  time (&(info -> program_end));

  if (info -> verbose) {
    total_time = difftime (info -> program_end, info -> program_start);
    if (total_time > 60) {
      fprintf (stderr, "==\tProgram execution:                              %.3f mins\n", total_time / 60);
    }
    else {
      fprintf (stderr, "==\tProgram execution:                              %.3f secs\n", total_time);
    }
    if (total_time > 1) {
      fprintf (stderr, "==\t  run() time:                                   %6.2f %%\n", info -> run_time / total_time * 100);

      fprintf (stderr, "==\t    Read data in:                               %6.2f %%\n", info -> readCO_time / total_time * 100);
      fprintf (stderr, "==\t    EM initialization:                          %6.2f %%\n", info -> initEM_time / total_time * 100);
      fprintf (stderr, "==\t    Calculate p(x,y):                           %6.2f %%\n", info -> calculateProbW1W2_time / total_time * 100);
      fprintf (stderr, "==\t    Calculate ML:                               %6.2f %%\n", info -> calculateML_time / total_time * 100);
      fprintf (stderr, "==\t    Swap previous and current:                  %6.2f %%\n", info -> swapPrevCurr_time / total_time * 100);
      fprintf (stderr, "==\t    Apply EM step:                              %6.2f %%\n", info -> applyEMStep_time / total_time * 100);
      fprintf (stderr, "==\t    Gather probabilities:                       %6.2f %%\n", info -> gatherProbs_time / total_time * 100);
      fprintf (stderr, "==\t    Normalize probabilities:                    %6.2f %%\n", info -> normalizeProbs_time / total_time * 100);
      fprintf (stderr, "==\t    Distribute probabilities:                   %6.2f %%\n", info -> distributeProbs_time / total_time * 100);
      fprintf (stderr, "==\t    Print probabilities:                        %6.2f %%\n", info -> printCoProbs_time / total_time * 100);
    }
  }

  wfree (info);

  return;
}


bool run (INFO *info) {
  PROBNODE curr_ML = 0;
  PROBNODE prev_ML = 0;
  PROBNODE diff = 0.0;
  int error_code;

  info -> iter = 0;
  error_code = 0;

  time_t start;
  time_t end;
  time_t loop_start;
  time_t loop_end;
  double timediff = 0.0;

  time (&start);

  /*  All processes read in co-occurrence data  */
  if (!readCO (info)) {
    /*  If there is an error, all processes are terminated  */
    fprintf (stderr, "Error reading co-occurrence data by processor %u.\n", info -> world_id);
#if HAVE_MPI
    MPI_Abort (MPI_COMM_WORLD, 0);
#endif
    return false;
  }

  /*  Only MAINPROC initializes to ensure the random seed only affects it  */
  if (info -> world_id == MAINPROC) {
    /*  Initial probabilties placed in *current*  */
    initEM (info);
    if (info -> verbose) {
      fprintf (stderr, "==\tm = %u; n = %u\n", info -> m, info -> n);
    }
  }

  /*  Send the initial probabilities in *current* for p(w1|z), p(w2|z), and p(z) to all processes  */
  distributeProbs (info);

  time (&loop_start);
  while (true) {
#if HAVE_MPI
    /*  Create a barrier at each iteration start */
    MPI_Barrier (MPI_COMM_WORLD);
#endif

    /*  Calculate p(w1, w2) using *current*  */
    calculateProbW1W2 (info);
    if (info -> world_id == MAINPROC) {
      /*  Calculate the log likelihood  */
      curr_ML = calculateML (info);

      if (info -> iter == 0) {
        if (info -> verbose) {
          fprintf (stderr, "[---]  Initial = %f\n", curr_ML);
        }
        if (info -> snapshot != UINT_MAX) {
          printCoProb (info);
        }
      }
      else {
        diff = (curr_ML - prev_ML) / prev_ML * 100 * -1;
        if (info -> verbose) {
          fprintf (stderr, "[%3u]  %f --> %f\t[%f, %2.4f %%]\n", info -> iter, prev_ML, curr_ML, (curr_ML - prev_ML), diff);
        }
        if ((curr_ML < prev_ML) || (DBL_LESS (fabs (diff), ML_DELTA))) {
          info -> iter = UINT_MAX;  /*  Set an indicator to leave loop  */
        }
      }

      prev_ML = curr_ML;

#if DEBUG
        checkCoProb (info);
#endif

      if (info ->  iter != UINT_MAX) {
        info -> iter++;
      }

      if (info -> iter > (info -> maxiter)) {
        info -> iter = UINT_MAX;  /*  Set an indicator to leave loop  */
      }
    }

#if HAVE_MPI
    /*  Broadcast the iteration number to all processes  */
    error_code = MPI_Bcast (&(info -> iter), 1, MPI_TYPE, MAINPROC, MPI_COMM_WORLD);
    if (error_code != MPI_SUCCESS) {
      fprintf (stderr, "Broadcast iteration from %u result:  %d.\n", info -> world_id, error_code);
    }

    /*  Broadcast p(i,j) from main to all other processes  */
    error_code = MPI_Bcast (info -> prob_w1w2, (info -> m * info -> n), MPI_TYPE, MAINPROC, MPI_COMM_WORLD);
    if (error_code != MPI_SUCCESS) {
      fprintf (stderr, "Second broadcast p(x,y) from %u result:  %d.\n", info -> world_id, error_code);
    }
#endif

    /*  Check if we are suppose to exit this loop  */
    if (info -> iter == UINT_MAX) {
      break;
    }

    /*  Swap the previous with current; *previous* is used to overwrite *current*  */
    swapPrevCurr (info);

    /*  Calculate E- and M-steps together; place results in *current*  */
    applyEMStep (info);

    /*  Transmit *current* to MAINPROC  */
    gatherProbs (info);

    /*  MAINPROC normalizes probabilities in *current*; and decide if a temporary snapshot should be printed  */
    if (info -> world_id == MAINPROC) {
      normalizeProbs (info);
      /*  If snapshots are required, then print it out if this is the first iteration OR
      **  this iteration is a multiple of (info -> snapshot)  */
      if ((info -> snapshot != UINT_MAX) &&
          ((info -> iter % info -> snapshot == 0) || (info -> iter == 1))) {
        if (!info -> no_output) {
          printCoProb (info);
        }
      }
    }

    distributeProbs (info);
  }
  time (&loop_end);
  timediff += difftime (loop_end, loop_start);

  if (info -> maxiter == 1) {
    fprintf (stderr, "==\t  Main loop [one iteration only!]:             %6.2f %% (%f)\n", 0.0, timediff);
  }

  /*  Only MAINPROC prints  */
  if (info -> world_id == MAINPROC) {
    if (!info -> no_output) {
      printCoProb (info);
    }
  }

  time (&end);
  info -> run_time += difftime (end, start);

  return (true);
}

