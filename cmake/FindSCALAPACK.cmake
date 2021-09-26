#
# Module to find the library SCALAPACK
#

#
# If SCALAPACK is found, it will set the following variables. Otherwise, 
# SCALAPACK_FOUND will be set to false
#
#  SCALAPACK_FOUND        True if MUMPS is found
#  SCALAPACK_LIBRARIES    METIS_librarie

#set (SCALAPACK_DIR "/NOBACKUP/frachon/lib/gcc/scalapack-2.0.2/build/lib/")
#set (SCALAPACK_DIR "/afs/kth.se/home/f/r/frachon/lib/gcc/lib")


set(SCALAPACK_FOUND "NO")

#if SCALAPACK_DIR is specified
if(SCALAPACK_DIR)
  find_path(SCALAPACK_LIBRARY_DIR 
#    NAMES libscalapack.so
    NAMES libscalapack-openmpi.so
    PATHS ${SCALAPACK_DIR} 
    NO_DEFAULT_PATH)
endif()

# otherwise look for standard places
find_path(SCALAPACK_LIBRARY_DIR
#  NAMES libscalapack.so
NAMES libscalapack-openmpi.so
  PATHS 
  /opt/SCALAPACK/lib
  /usr/local/SCALAPACK/lib
  /usr/local/lib
  /usr/lib
  ~/lib/lib)

if(SCALAPACK_LIBRARY_DIR)
  set(SCALAPACK_FOUND YES)

  find_library(SCALAPACK_LIBRARY 
#    NAMES scalapack
NAMES scalapack-openmpi
    PATHS ${SCALAPACK_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})

  mark_as_advanced(SCALAPACK_LIBRARY_DIR SCALAPACK_LIBRARY)
    message( "-- SCALAPACK_library FOUND")
else()
  if(SCALAPACK_FIND_REQUIRED)
    message( "SCALAPACK_library = ${SCALAPACK_LIBRARY}")
    message(FATAL_ERROR "SCALAPACK not found, please set SCALAPACK_DIR to your SCALAPACK install directory")
  endif()
endif()
