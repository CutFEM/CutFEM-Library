#
# Module to find the library METIS
#

#
# If METIS is found, it will set the following variables. Otherwise, 
# METIS_FOUND will be set to false
#
#  METIS_FOUND        True if MUMPS is found
#  METIS_LIBRARIES    METIS_librarie

set(METIS_INCLUDES "/home/f/r/frachon/lib/metis/5.1.0/include")
set(METIS_LIBRARIES "/home/f/r/frachon/lib/metis/5.1.0/build/libmetis/libmetis.a")
set(METIS_FOUND YES)

message( "METIS_include   = ${METIS_INCLUDES}")
message( "METIS_libraries = ${METIS_LIBRARIES}")
