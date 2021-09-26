#
# Module to find the library BLACS
#

#
# If BLACS is found, it will set the following variables. Otherwise, 
# BLACS_FOUND will be set to false
#
#  BLACS_FOUND        True if MUMPS is found
#  BLACS_LIBRARIES    METIS_librarie

set (BLACS_DIR "/afs/kth.se/home/f/r/frachon/lib/gcc/lib")

set(BLACS_FOUND "NO")

#if BLACS_DIR is specified
if(BLACS_DIR)
  find_path(BLACS_LIBRARY_DIR 
    NAMES blacs_MPI-LINUX-0.a
    PATHS ${BLACS_DIR} 
    NO_DEFAULT_PATH)
endif()

# otherwise look for standard places
find_path(BLACS_LIBRARY_DIR
  NAMES blacs_MPI-LINUX-0.a
  PATHS 
  /opt/BLACS/lib
  /usr/local/BLACS/lib
  /usr/local/lib
  /usr/lib
  ~/lib/lib)

if(BLACS_LIBRARY_DIR)
  set(BLACS_FOUND YES)

  find_library(BLACS_LIBRARY 
    NAMES blacs_MPI-LINUX-0.a
    PATHS ${BLACS_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  find_library(BLACS_Cinit_LIBRARY 
    NAMES blacsCinit_MPI-LINUX-0.a
    PATHS ${BLACS_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  find_library(BLACS_F77init_LIBRARY 
    NAMES blacsF77init_MPI-LINUX-0.a
    PATHS ${BLACS_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  set(BLACS_LIBRARIES ${BLACS_LIBRARY} ${BLACS_Cinit_LIBRARY} ${BLACS_F77init_LIBRARY})

  mark_as_advanced(BLACS_LIBRARY_DIR BLACS_LIBRARY BLACS_Cinit_LIBRARY BLACS_F77init_LIBRARY)
else()
  if(BLACS_FIND_REQUIRED)
    message( "BLACS_library = ${BLACS_LIBRARY}")
    message(FATAL_ERROR "BLACS not found, please set BLACS_DIR to your BLACS install directory")
  endif()
endif()
