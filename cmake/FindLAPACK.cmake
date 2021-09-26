#
# Module to find the library LAPACK
#

#   LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#     is found
#   LAPACK_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#     and -L).
#   LAPACK_LIBRARIES - uncached list of libraries (using full path name) to
#     link against to use LAPACK

# set(LAPACK_LIBRARIES_DIR "/home/f/r/frachon/lib/lapack-3.7.0/")
set(LAPACK_FOUND "NO")

find_path(LAPACK_INCLUDES
NAMES
lapacke.h
PATHS
/home/f/r/frachon/lib/lapack-3.7.0/LAPACKE/include/
/usr/local/opt/lapack/include/
/usr/include/
/usr/lib/lapack
)

##message(" LAPACK INCLUDE :  ${LAPACK_INCLUDES}")

# otherwise look for standard places
find_path(LAPACK_LIBRARY_DIR
NAMES liblapacke.dylib
PATHS
/home/f/r/frachon/lib/lapack-3.7.0/
/usr/lib/lapack
/usr/local/opt/lapack/lib/
/opt
/usr/local/lib
/usr/lib
~/lib/lib)

##message(" LAPACK LIB DIR :  ${LAPACK_LIBRARY_DIR}")

if(LAPACK_LIBRARY_DIR)
set(LAPACK_FOUND YES)

find_library(LAPACK_LIBRARY
NAMES lapacke
PATHS ${LAPACK_LIBRARY_DIR}
NO_DEFAULT_PATH)

set(LAPACK_LIBRARIES ${LAPACK_LIBRARY})
set(LAPACK_MODULES_DIR ${LAPACK_LIBRARY_DIR})
##message( "LAPACK_library = ${LAPACK_LIBRARIES}")

mark_as_advanced(LAPACK_LIBRARY_DIR LAPACK_LIBRARY LAPACK_MODULES_DIR)
else()
if(LAPACK_FIND_REQUIRED)
##message( "LAPACK_library = ${LAPACK_LIBRARY}")
message(FATAL_ERROR "LAPACK not found, please set LAPACK_DIR to your LAPACK install directory")
endif()
endif()
