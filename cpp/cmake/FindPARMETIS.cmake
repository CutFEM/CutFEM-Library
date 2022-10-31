#
# Module to find the library METIS
#

#
# If METIS is found, it will set the following variables. Otherwise, 
# METIS_FOUND will be set to false
#
#  METIS_FOUND        True if MUMPS is found
#  METIS_LIBRARIES    METIS_librarie

set(PARMETIS_INCLUDES "/home/f/r/frachon/lib/parmetis/4.0.3/include")
set(PARMETIS_LIBRARIES "/home/f/r/frachon/lib/parmetis/4.0.3/build/Linux-x86_64/libparmetis/libparmetis.a")
set(PARMETIS_FOUND YES)

message( "PARMETIS_include   = ${PARMETIS_INCLUDES}")
message( "PARMETIS_libraries = ${PARMETIS_LIBRARIES}")
