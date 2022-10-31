# - Try to find PETSc
# Once done this will define
#
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDES     - the PETSc include directories
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_COMPILER     - Compiler used by PETSc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#  PETSC_MPIEXEC      - Executable for running MPI programs
#  PETSC_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)
#
#  Usage:
#  find_package(PETSc COMPONENTS CXX)  - required if build --with-clanguage=C++ --with-c-support=0
#  find_package(PETSc COMPONENTS C)    - standard behavior of checking build using a C compiler
#  find_package(PETSc)                 - same as above
#
# Setting these changes the behavior of the search
#  PETSC_DIR - directory in which PETSc resides
#  PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#



set(PETSC_DIR /usr/lib/petscdir/3.7.3/x86_64-linux-gnu-real )
set(PETSC_ARCH  )
#set(PETSC_INCLUDES ${PETSC_DIR}/include)
#set(PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc_real.so.3.7.3)
#set(PETSC_MPIEXE ${PETSC_DIR}/bin/petscmpiexec)
#set (PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSc definitions" FORCE)



set (PETSC_EXECUTABLE_RUNS YES CACHE BOOL "Disable checking if this setup works" FORCE)
 
set (PETSC_FOUND YES CACHE BOOL "PETSc was found (manually)" FORCE)
 
set (PETSC_INCLUDES ${PETSC_DIR}/include CACHE STRING
  "Semicolon-delimited list of PETSc include directories" FORCE)
 
set (PETSC_LIBRARIES ${PETSC_DIR}/lib/libpetsc_real.so.3.7.3 CACHE STRING
  "Semicolon-delimited list of PETSc libraries" FORCE)
 
set (PETSC_COMPILER "/usr/lib/openmpi/lib/" CACHE FILEPATH
  "PETSc compiler; helpful to find a compatible MPI" FORCE)
 
set (PETSC_DEFINITIONS "-D__INSDIR__="  CACHE STRING
  "PETSc definitions" FORCE)
 
set (PETSC_MPIEXEC ${PETSC_DIR}"/bin/petscmpiexec" CACHE FILEPATH
  "Executable for running PETSc MPI programs" FORCE)
 
set (PETSC_VERSION "" CACHE STRING
  "PETSc version: MAJOR.MINOR.SUBMINOR" FORCE)
 
mark_as_advanced (PETSC_INCLUDES PETSC_LIBRARIES
  PETSC_COMPILER PETSC_DEFINITIONS
  PETSC_MPIEXEC PETSC_EXECUTABLE_RUNS PETSC_VERSION)


message("dir  = " ${PETSC_DIR})
message("arch  = " ${PETSC_ARCH})
message("includes  = " ${PETSC_INCLUDES})
message("librarie  = " ${PETSC_LIBRARIES})
message("compiler  = " ${PETSC_COMPILER})
message("definition  = " ${PETSC_DEFINITIONS})
message("mpiexe  = " ${PETSC_MPIEXE})
message("version  = " ${PETSC_VERSION})
