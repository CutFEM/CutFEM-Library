#
# Module to find the library SCOTCH
#

#
# If SCOTTH is found, it will set the following variables. Otherwise, 
# SCOTCH_FOUND will be set to false
#
#  SCOTCH_FOUND        True if SCOTCH is found
#  SCOTCH_LIBRARIES    SCOTCH_librarie

#set (SCOTCH_DIR "/afs/kth.se/home/f/r/frachon/lib/scotch_6.0.4/lib")
set(SCOTCH_DIR "/opt/homebrew/Cellar/scotch/7.0.1")
set(SCOTCH_FOUND "NO")

#if SCOTCH_DIR is specified
if(SCOTCH_DIR)
  find_path(SCOTCH_LIBRARY_DIR 
    NAMES libesmumps.a libscotch.a libscotcherr.a
    PATHS ${SCOTCH_DIR}/lib 
    NO_DEFAULT_PATH)
endif()

# otherwise look for standard places
find_path(SCOTCH_LIBRARY_DIR
  NAMES libesmumps.a libscotch.a libscotcherr.a
  PATHS
  /opt/SCOTCH/lib
  /opt/homebrew/Cellar/scotch/7.0.1
  /usr/local/lib
  /usr/lib
  ~/lib/lib)

if(SCOTCH_LIBRARY_DIR)
  set(SCOTCH_FOUND YES)

  find_library(SCOTCH_esmumps_LIBRARY 
    NAMES esmumps
    PATHS ${SCOTCH_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  find_library(SCOTCH_scotch_LIBRARY 
    NAMES scotch
    PATHS ${SCOTCH_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  find_library(SCOTCH_scotcherr_LIBRARY 
    NAMES scotcherr
    PATHS ${SCOTCH_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  set(SCOTCH_LIBRARIES ${SCOTCH_scotcherr_LIBRARY} ${SCOTCH_scotch_LIBRARY} ${SCOTCH_esmumps_LIBRARY})

  mark_as_advanced(SCOTCH_LIBRARY_DIR SCOTCH_scotcherr_LIBRARY SCOTCH_scotch_LIBRARY SCOTCH_esmumps_LIBRARY)
    message( "-- SCOTCH_library FOUND")
else()
  if(SCOTCH_FIND_REQUIRED)
    message( "SCOTCH_library = ${SCOTCH_LIBRARY}")
    message(FATAL_ERROR "SCOTCH not found, please set SCOTCH_DIR to your SCOTCH install directory")
  endif()
endif()
