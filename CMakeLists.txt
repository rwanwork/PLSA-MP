################################################################################
##  Probabilistic latent semantic analysis (PLSA, multiprocessor version)
##  Copyright (C) 2021  by Raymond Wan (rwan.work@gmail.com)
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

##  Set the minimum required CMake version
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

##  Define the project
project (plsa-mp LANGUAGES C)

##  List of C source files
set (PLSA_MP_SRCFILES 
  comm.c
  debug.c
  em-steps.c
  input.c
  main.c
  output.c
  parameters.c
  run.c
  wmalloc.c
)


########################################
##  Detect OpenMP and MPI

find_package (MPI)
if (MPI_FOUND)
  set (HAVE_MPI 1)
endif (MPI_FOUND)

find_package (OpenMP)
if (OPENMP_FOUND)
  set (HAVE_OPENMP 1)
endif (OPENMP_FOUND)


########################################
##  Set various values based on libraries found

IF (MPI_FOUND)
  INCLUDE_DIRECTORIES (${MPI_INCLUDE_PATH})
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
  SET (LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
  SET (CMAKE_C_COMPILER "${MPI_COMPILER}")
ENDIF (MPI_FOUND)

IF (OPENMP_FOUND)
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
ENDIF (OPENMP_FOUND)


########################################
##  Create the configuration file

##  Configure a header file to pass some of the CMake settings
##  to the source code
configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/PLSA_MP_Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/PLSA_MP_Config.h"
  )

##  Include the binary directory so that the created configuration 
##    file can be located
include_directories (${CMAKE_CURRENT_BINARY_DIR})
  

########################################
##  Define the executable

add_executable (plsa ${PLSA_MP_SRCFILES})

##  Link the executable to the math library
target_link_libraries (plsa m)

