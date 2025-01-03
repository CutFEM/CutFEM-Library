#
# Module to find the library MUMPS
#

#
# If MUMPS is found, it will set the following variables. Otherwise,
# MUMPS_FOUND will be set to false
#
#  MUMPS_FOUND        True if MUMPS is found
#  MUMPS_LIBRARIES    MUMPS_librarie
#  MUMPS_INCLUDE_DIR  where to find mumps_compat.h

set(MUMPS_FOUND "NO")

#if MUMPS_DIR is specified
if(MUMPS_DIR)
  set(MUMPS_INCLUDE_DIR ${MUMPS_DIR}/include)
  set(MUMPS_LIBRARY_DIR ${MUMPS_DIR}/lib)
else()
find_path (MUMPS_INCLUDE_DIR
  NAMES dmumps_struc.h
  PATHS
  /usr/local/Cellar/brewsci-mumps/5.2.1/include
  /opt/homebrew/Cellar/brewsci-mumps/5.3.5/include
  # /Users/thomasfrachon/lib/MUMPS_5.6.2/include
  /opt/homebrew/Cellar/brewsci-mumps/5.6.2/include
  /usr/local/MUMPS/include
  /usr/local/include
  /opt/MUMPS/include
  /usr/include/MUMPS
  /usr/include
  ~/lib/include)

if(USE_MPI)
message(" Search for MPI MUMPS")
find_path(MUMPS_LIBRARY_DIR
  NAMES libmumps_common.dylib libdmumps.dylib libpord.dylib
  libmumps_common.a libdmumps.a libpord.a
#  NAMES libmumps_common.a libdmumps.a libpord.a

  PATHS
  /usr/local/Cellar/brewsci-mumps/5.2.1/lib
  /opt/homebrew/Cellar/brewsci-mumps/5.3.5/lib
  # /Users/thomasfrachon/lib/MUMPS_5.6.2/lib
  /opt/homebrew/Cellar/brewsci-mumps/5.6.2/lib
  /opt/MUMPS/lib
  /usr/local/MUMPS/lib
  /usr/local/lib
  /usr/lib/x86_64-linux-gnu
  /usr/lib
  ~/lib/lib)
else()
message(" Search for SEQ MUMPS")
find_path(MUMPS_LIBRARY_DIR
  NAMES 
  libmumps_common_seq.a libdmumps_seq.a libpord_seq.a
# libmumps_common.a libdmumps.a libpord.a
  PATHS
  /usr/local/Cellar/brewsci-mumps/5.2.1/lib
  /opt/homebrew/Cellar/brewsci-mumps/5.3.5/lib
  # /Users/thomasfrachon/lib/MUMPS_5.6.2/libs
  /opt/homebrew/Cellar/brewsci-mumps/5.6.2/lib
  /opt/MUMPS/lib
  /usr/local/MUMPS/lib
  /usr/local/lib
  /usr/lib/x86_64-linux-gnu/
  /usr/lib
  ~/lib/lib)

  message(" library directory = ${MUMPS_LIBRARY_DIR}")

endif()

endif()

if(MUMPS_INCLUDE_DIR AND MUMPS_LIBRARY_DIR)
  set(MUMPS_FOUND YES)

  if(USE_MPI)
  find_library(MUMPS_COMMON_LIBRARY
    NAMES mumps_common
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)

  find_library(MUMPS_D_LIBRARY
    NAMES dmumps
    PATHS ${MUMPS_LIBRARY_DIR}
    NO_DEFAULT_PATH)
  find_library(MUMPS_PORD_LIBRARY
    NAMES pord
    PATHS ${MUMPS_LIBRARY_DIR}
    /usr/lib/x86_64-linux-gnu
    NO_DEFAULT_PATH)    
  else()
  find_library(MUMPS_COMMON_LIBRARY
    NAMES mumps_common_seq
    PATHS ${MUMPS_LIBRARY_DIR}
    /usr/lib/x86_64-linux-gnu
    NO_DEFAULT_PATH)

  find_library(MUMPS_D_LIBRARY
    NAMES dmumps_seq
    PATHS ${MUMPS_LIBRARY_DIR}
    /usr/lib/x86_64-linux-gnu
    NO_DEFAULT_PATH)
  find_library(MUMPS_PORD_LIBRARY
    NAMES pord_seq
    PATHS ${MUMPS_LIBRARY_DIR}
    /usr/lib/x86_64-linux-gnu
    NO_DEFAULT_PATH)
  endif()


  find_library(MUMPS_PARMETIS_LIBRARY
    NAMES parmetis
    PATHS 
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    /usr/local/Cellar/brewsci-parmetis/4.0.3_1/lib
		/opt/homebrew/Cellar/brewsci-parmetis/4.0.3_1/lib
    NO_DEFAULT_PATH)

#  set(SCOTCH_LIBRARY_DIR /usr/lib )
#  find_library(SCOTCH_esmumps_LIBRARY
#    NAMES esmumps
#    PATHS ${SCOTCH_LIBRARY_DIR}
#    NO_DEFAULT_PATH)

  find_library(SCOTCH_scotch_LIBRARY
    NAMES scotch scotch-6
    PATHS /usr/lib
    /opt/homebrew/Cellar/scotch/7.0.2/lib
		/opt/homebrew/Cellar/scotch/7.0.4/lib
    /opt/homebrew/Cellar/scotch/7.0.5/lib
    /usr/lib/x86_64-linux-gnu
    /usr/local/Cellar/brewsci-scotch/6.0.4/lib
    NO_DEFAULT_PATH)

  find_library(SCOTCH_scotcherr_LIBRARY
    NAMES scotcherr scotcherr-6
    PATHS /usr/lib
    /usr/lib/x86_64-linux-gnu
		/opt/homebrew/Cellar/scotch/7.0.2/lib
		/opt/homebrew/Cellar/scotch/7.0.4/lib
    /opt/homebrew/Cellar/scotch/7.0.5/lib
		/usr/local/Cellar/brewsci-scotch/6.0.4/lib
    NO_DEFAULT_PATH)

  set(SCOTCH_LIBRARIES ${SCOTCH_scotcherr_LIBRARY} ${SCOTCH_scotch_LIBRARY})
  # ${SCOTCH_esmumps_LIBRARY})

  set(MUMPS_LIBRARIES ${MUMPS_COMMON_LIBRARY} ${MUMPS_D_LIBRARY}   ${SCOTCH_LIBRARIES} ${MUMPS_PARMETIS_LIBRARY} ${MUMPS_PORD_LIBRARY})

  set(MUMPS_INCLUDES ${MUMPS_INCLUDE_DIR})

  message( "-- MUMPS_library FOUND")
  message( " MUMPS_includes = ${MUMPS_INCLUDES}")
  message( " MUMPS_library = ${MUMPS_LIBRARIES}")

else()
  if(MUMPS_FIND_REQUIRED)
    message( "MUMPS_include_dirs = ${MUMPS_INCLUDE_DIR}")
    message( "MUMPS_library = ${MUMPS_LIBRARIES}")
    message(FATAL_ERROR "MUMPS not found, please set MUMPS_DIR to your MUMPS install directory")
  endif()
endif()