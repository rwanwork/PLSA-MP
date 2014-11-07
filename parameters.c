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

#define _GNU_SOURCE
#include <getopt.h>

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>                                  /*  UINT_MAX  */
#include <string.h>
#include <stdbool.h>

#include "config.h"
#if HAVE_MPI
#include <mpi.h>
#endif

#if HAVE_OPENMP
#include <omp.h>
#endif

#include "wmalloc.h"
#include "plsa-defn.h"
#include "parameters.h"

/*!  Print out usage information  */
void usage (char *progname) {
  fprintf (stderr, "Probabilistic Latent Semantic Analysis (multiprocessor)\n");
  fprintf (stderr, "=======================================================\n\n");
  fprintf (stderr, "Usage:  %s [options]\n\n", progname);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "--base <file>      :  Base filename for output file.\n");
  fprintf (stderr, "--cooccur <file>   :  Co-occurrence filename.\n");
  fprintf (stderr, "--clusters <int>   :  Number of clusters.\n");
  fprintf (stderr, "--seed <int>       :  Random seed.\n");
  fprintf (stderr, "                   :    (Default:  current time).\n");
  fprintf (stderr, "--maxiter <int>    :  Maximum iterations.\n");
  fprintf (stderr, "--text             :  Text mode (I/O is in text, not binary).\n");
  fprintf (stderr, "--snapshot <int>   :  Output snapshots p(x,y) at regular intervals.\n");
  fprintf (stderr, "                   :    (Default:  Do not output).\n");
  fprintf (stderr, "--openmp <int>     :  Number of OpenMP threads to use.\n");
  fprintf (stderr, "                   :    (Default:  Maximum for PC).\n");
  fprintf (stderr, "--verbose          :  Verbose mode.\n");
  fprintf (stderr, "--debug            :  Debugging output.\n");
  fprintf (stderr, "--rounding         :  Round using %u as the multiplication factor.\n", ROUND_DIGITS);
  fprintf (stderr, "--nooutput         :  Suppress outputting p(x,y) to file.\n");

  fprintf (stderr, "\nCompile-time settings:\n  ");
  fprintf (stderr, "     MPI:                              ");
#if HAVE_MPI
  fprintf (stderr, "Enabled\n");
#else
  fprintf (stderr, "Disabled\n");
#endif
  fprintf (stderr, "    OpenMP:                              ");
#if HAVE_OPENMP
  fprintf (stderr, "Enabled\n");
#else
  fprintf (stderr, "Disabled\n");
#endif

  fprintf (stderr, "\nPLSA version:  %s (%s)\n\n", __DATE__, __TIME__);

  exit (EXIT_SUCCESS);
}


bool checkSettings (INFO *info) {
  if (info -> co_fn == NULL) {
    fprintf (stderr, "==\tError:  Co-occurrence filename required with the --cooccur option.\n");
    return false;
  }

  if (info -> maxiter == 0) {
    fprintf (stderr, "==\tError:  Maximum number of iterations required with the --maxiter option.\n");
    return false;
  }

  if (info -> num_clusters == 0) {
    fprintf (stderr, "==\tError:  Number of clusters required with the --clusters option.\n");
    return false;
  }

  if (info -> num_clusters >= MAX_CLUSTERS) {
    fprintf (stderr, "==\tError:  The maximum number of clusters is %u because of the MPI message tag lengths.  Please change MAX_CLUSTERS in plsa-defn.h and re-compile.\n", MAX_CLUSTERS);
    return false;
  }

  if (info -> base_fn == NULL) {
    fprintf (stderr, "==\tError:  Base filename required with the --base option.\n");
    return false;
  }

  if (info -> world_size > info -> num_clusters) {
    fprintf (stderr, "==\tWarning:  The number of processors is more than the number of clusters.  Increasing the number of clusters.");
    info -> num_clusters = info -> world_size;
  }

  if (info -> verbose) {
    if (info -> world_id == MAINPROC) {
      fprintf (stderr, "Settings\n");
      fprintf (stderr, "--------\n");
      fprintf (stderr, "==\tBase filename:                                  %s\n", info -> base_fn);
      fprintf (stderr, "==\tCo-occurrence filename:                         %s\n", info -> co_fn);
      fprintf (stderr, "==\tProbability data type:                          ");
      if (sizeof (PROBNODE) == 4) {
        fprintf (stderr, "float\n");
      }
      else if (sizeof (PROBNODE) == 8) {
        fprintf (stderr, "double\n");
      }
      else {
        fprintf (stderr, "Unknown!\n");
      }
      fprintf (stderr, "==\tClusters:                                       %u\n", info -> num_clusters);
      if (info -> seed != UINT_MAX) {
        fprintf (stderr, "==\tRandom seed:                                    %u\n", info -> seed);
      }
      else {
        fprintf (stderr, "==\tRandom seed:                                    [from time]\n");
      }
      fprintf (stderr, "==\tExponent difference [utils.h::addLogsFloat]:    %.8f\n", LN_LIMIT);
      fprintf (stderr, "==\tTermination conditions\n");
      fprintf (stderr, "==\t  Maximum EM iterations:                        %u\n", info -> maxiter);
      fprintf (stderr, "==\t  Percentage difference:                        %f\n", ML_DELTA);
      fprintf (stderr, "==\tText mode:                                      %s\n", (info -> textio) ? "yes" : "no");
      fprintf (stderr, "==\tRounding:                                       %s\n", (info -> rounding) ? "yes" : "no");
      if (info -> rounding) {
        fprintf (stderr, "==\tRounding factor:                                %u\n", ROUND_DIGITS);
      }
      fprintf (stderr, "==\tSuppress output to file:                        %s\n", (info -> no_output) ? "yes" : "no");
    }
#if HAVE_MPI
    fprintf (stderr, "==\tMPI:                                            OK\n");
    fprintf (stderr, "==\t  My ID:                                        %d\n", info -> world_id);
    fprintf (stderr, "==\t  Number of processes:                          %d\n", info -> world_size);
    fprintf (stderr, "==\t  Block range:                                  %u - %u\n", info -> block_start, info -> block_end);
    fprintf (stderr, "==\t  Block size:                                   %u\n", info -> block_size);
#else
    fprintf (stderr, "==\tMPI information:                                Not enabled\n");
#endif
#if HAVE_OPENMP
    fprintf (stderr, "==\tOpen MP:                                        OK\n");
    fprintf (stderr, "==\t  Number of threads:                            %u\n", info -> threads);
#else
    fprintf (stderr, "==\tOpen MP:                                        Not enabled\n");
#endif
    fprintf (stderr, "\n\n");
  }

  return true;
}


bool processOptions (int argc, char *argv[], INFO *info) {
  int c = 0;

  char *base_fn = NULL;
  char *co_fn = NULL;
  unsigned int num_clusters = 0;
  unsigned int seed = UINT_MAX;
  unsigned int maxiter = 0;
  unsigned int snapshot = UINT_MAX;
  bool verbose = false;
  bool debug = false;
  bool textio = false;
  bool rounding = false;
  bool no_output = false;

  /*  Usage information if no arguments  */
  if (argc == 1) {
    usage (argv[0]);
  }

  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"base", 1, 0, 0},
      {"cooccur", 1, 0, 0},
      {"clusters", 1, 0, 0},
      {"seed", 1, 0, 0},
      {"maxiter", 1, 0, 0},
      {"snapshot", 1, 0, 0},
      {"openmp", 1, 0, 0},
      {"verbose", 0, 0, 0},
      {"debug", 0, 0, 0},
      {"text", 0, 0, 0},
      {"rounding", 0, 0, 0},
      {"nooutput", 0, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "", long_options, &option_index);
    if (c == -1) {
      break;
    }

    switch (c) {
      case 0:
        if (strcmp (long_options[option_index].name, "cooccur") == 0) {
          co_fn = wmalloc (strlen (optarg) + 1);
          co_fn = strcpy (co_fn, optarg);
        }
        else if (strcmp (long_options[option_index].name, "clusters") == 0) {
          num_clusters = atoi (optarg);
        }
        else if (strcmp (long_options[option_index].name, "seed") == 0) {
          seed = atoi (optarg);
        }
        else if (strcmp (long_options[option_index].name, "base") == 0) {
          base_fn = wmalloc (strlen (optarg) + 1);
          base_fn = strcpy (base_fn, optarg);
        }
        else if (strcmp (long_options[option_index].name, "maxiter") == 0) {
          maxiter = atoi (optarg);
        }
        else if (strcmp (long_options[option_index].name, "snapshot") == 0) {
          snapshot = atoi (optarg);
        }
        else if (strcmp (long_options[option_index].name, "openmp") == 0) {
#if HAVE_OPENMP
          /*  Check the previously set value, which is the maximum for the system  */
          if (atoi (optarg) > info -> threads) {
            fprintf (stderr, "==\tError:  The number of threads requested exceeds the number available in the system (%u).\n", info -> threads);
            exit (-1);
          }
          info -> threads = atoi (optarg);
          omp_set_num_threads (info -> threads);
#else
          fprintf (stderr, "==\tError:  OpenMP is not enabled; --openmp meaningless.\n");
          exit (-1);
#endif
        }
        else if (strcmp (long_options[option_index].name, "verbose") == 0) {
          verbose = true;
        }
        else if (strcmp (long_options[option_index].name, "debug") == 0) {
          debug = true;
        }
        else if (strcmp (long_options[option_index].name, "text") == 0) {
          textio = true;
        }
        else if (strcmp (long_options[option_index].name, "rounding") == 0) {
          rounding = true;
        }
        else if (strcmp (long_options[option_index].name, "nooutput") == 0) {
          no_output = true;
        }
        break;
      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
        exit (EXIT_FAILURE);
    }
  }

  info -> base_fn = base_fn;
  info -> co_fn = co_fn;
  info -> num_clusters = num_clusters;
  info -> seed = seed;
  info -> maxiter = maxiter;
  info -> snapshot = snapshot;
  info -> verbose = verbose;
  info -> debug = debug;
  info -> textio = textio;
  info -> rounding = rounding;
  info -> no_output = no_output;

  /*  Set the range of clusters this process will handle  */
  info -> block_start = BLOCK_LOW (info ->  world_id, info -> world_size, info -> num_clusters);
  info -> block_end = BLOCK_HIGH (info ->  world_id, info -> world_size, info -> num_clusters);
  info -> block_size = BLOCK_SIZE (info ->  world_id, info -> world_size, info -> num_clusters);

  return true;
}

