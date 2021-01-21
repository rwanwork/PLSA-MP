###########################################################################
##  PLSA-MP
##  Probabilistic latent semantic analysis (PLSA, multiprocessor version)
##  
##  Copyright (C) 2009-2021 by Raymond Wan, All rights reserved.
##  Contact:  rwan.work@gmail.com
##  
##  This file is part of PLSA-MP.
##  
##  PLSA-MP is free software; you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License 
##  as published by the Free Software Foundation; either version 
##  3 of the License, or (at your option) any later version.
##  
##  PLSA-MP is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU Lesser General Public License for more details.
##  
##  You should have received a copy of the GNU General Public 
##  License along with PLSA-MP; if not, see 
##  <http://www.gnu.org/licenses/>.
###########################################################################

##  Set the minimum required CMake version
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)


########################################
##  Define some variables

##  Testing compilation flags, some of which are suggested by the Valgrind 3.3 >
# set (MY_C_FLAGS "-pedantic -Wno-long-long -g -fno-inline -O0 -Wall")
##  Release compilation flags, suggested by the Valgrind 3.3 book
set (MY_C_FLAGS "-O3 -Wall -Wno-unused-result -Wno-unused-but-set-variable")

set (TARGET_NAME_EXEC "plsa")
set (CURR_PROJECT_NAME "PLSA-MP")

##  Define the project
project (${CURR_PROJECT_NAME} LANGUAGES C)


########################################
##  Define the source files
##  Source files for both the test executable and library
set (SRC_FILES
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

if (MPI_FOUND)
  include_directories (${MPI_INCLUDE_PATH})
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
  set (LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
  set (CMAKE_C_COMPILER "${MPI_COMPILER}")
endif (MPI_FOUND)

if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif (OPENMP_FOUND)


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
##  Set initial compilation flags

##  Set compiler flags based on global variable
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")


########################################
##  Define the executable

add_executable (${TARGET_NAME_EXEC} ${SRC_FILES})

##  Link the executable to the math library
target_link_libraries (${TARGET_NAME_EXEC} m)

