#
# Module to find the library PORD
#

#
# If PORD is found, it will set the following variables. Otherwise, 
# PORD_FOUND will be set to false
#
#  PORD_FOUND        True if PORD is found
#  PORD_LIBRARIES    PORD_librarie


set(PORD_FOUND "NO")

#if PORD_DIR is specified
if(PORD_DIR)
  find_path(PORD_LIBRARY_DIR 
    NAMES libpord.so
    PATHS ${PORD_DIR}/lib 
    NO_DEFAULT_PATH)
endif()

# otherwise look for standard places
find_path(PORD_LIBRARY_DIR
  NAMES libpord.so
  PATHS 
  /opt/PORD/lib
  /usr/local/PORD/lib
  /usr/local/lib
  /usr/lib
  ~/lib/lib)

if(PORD_LIBRARY_DIR)
  set(PORD_FOUND YES)

  find_library(PORD_LIBRARY 
    NAMES pord
    PATHS ${PORD_LIBRARY_DIR} 
    NO_DEFAULT_PATH)

  set(PORD_LIBRARIES ${PORD_LIBRARY})

  mark_as_advanced(PORD_LIBRARY_DIR PORD_LIBRARY)
else()
  if(PORD_FIND_REQUIRED)
    message( "PORD_library = ${PORD_LIBRARY}")
    message(FATAL_ERROR "PORD not found, please set PORD_DIR to your PORD install directory")
  endif()
endif()
