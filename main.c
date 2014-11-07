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
#include <time.h>

#include "config.h"
#if HAVE_MPI
#include <mpi.h>
#endif

#if HAVE_OPENMP
#include <omp.h>
#endif

#include "plsa-defn.h"
#include "wmalloc.h"
#include "parameters.h"
#include "run.h"


/*!  Main function  */
int main (int argc, char *argv[]) {
  INFO *info;
  bool result = false;

  info = initialize ();

#if HAVE_MPI
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &(info -> world_id));
  MPI_Comm_size (MPI_COMM_WORLD, &(info -> world_size));
#endif

#if HAVE_OPENMP
  info -> threads = omp_get_num_procs ();
  omp_set_num_threads (info -> threads);
#endif

  /*  Process the command line parameters and then check them;
  **  if either fail then print usage information  */
  if ((!processOptions (argc, argv, info)) || (!checkSettings (info))) {
    usage (argv[0]);
  }
  else {
    result = run (info);
  }

#if HAVE_MPI
  MPI_Finalize ();
#endif

  uninitialize (info);

  return (EXIT_SUCCESS);
}

