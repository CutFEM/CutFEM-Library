# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

#set(INCLUDE_INSTALL_DIR "/NOBACKUP/frachon/lib/SuiteSparse/include")
#set(INCLUDE_INSTALL_DIR "/usr/local/Cellar/suite-sparse/5.6.0/include")
set(INCLUDE_INSTALL_DIR "/usr/include/suitesparse")

if (UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif (UMFPACK_INCLUDES AND UMFPACK_LIBRARIES)


find_path(UMFPACK_INCLUDES
  NAMES
  umfpack.h
  PATHS
  /usr/local/Cellar/suite-sparse/5.6.0/include
  /usr/include/suitesparse/
  /opt/homebrew/include
  $ENV{UMFPACKDIR}
  ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
  SuiteSparse
  ufsparse
  suitesparse
)


find_library(UMFPACK_LIBRARIES umfpack PATHS $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
if(UMFPACK_LIBRARIES)
  if (NOT UMFPACK_LIBDIR)
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARIES} PATH)
  endif(NOT UMFPACK_LIBDIR)
  find_library(COLAMD_LIBRARY colamd PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if (COLAMD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${COLAMD_LIBRARY})
  endif (COLAMD_LIBRARY)

  find_library(AMD_LIBRARY amd PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if (AMD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${AMD_LIBRARY})
  endif (AMD_LIBRARY)
  find_library(SUITESPARSE_LIBRARY SuiteSparse PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if (SUITESPARSE_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${SUITESPARSE_LIBRARY})
  endif (SUITESPARSE_LIBRARY)
endif(UMFPACK_LIBRARIES)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK DEFAULT_MSG
                                  UMFPACK_INCLUDES UMFPACK_LIBRARIES)

message( "umfpack_include = ${UMFPACK_INCLUDES}")
message( "umfpack_libraries = ${UMFPACK_LIBRARIES}")
