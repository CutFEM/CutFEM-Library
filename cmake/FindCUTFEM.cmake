set(CUTFEM_FOUND "NO")

set(CUTFEM_COMPONENENT libparallel.so 
                        libcommon.so  
                        libFESpace.so 
                        libsolver.so 
                        libproblem.so)

find_path(CUTFEM_INCLUDE_DIR
    NAMES
        cutfem.hpp
    PATHS
        /Users/thomasfrachon/Documents/code/cutfem/libcutfem/cpp/        
        /afs/kth.se/home/f/r/frachon/Documents/code/CutFEM-Library/cpp
)

find_path(CUTFEM_LIBRARY_DIR
    NAMES 
    libcommon.dylib
    libcommon.so
    PATHS
        /Users/thomasfrachon/Documents/code/cutfem/libcutfem/build/lib
        /afs/kth.se/home/f/r/frachon/Documents/code/CutFEM-Library/build/lib
)

if(CUTFEM_LIBRARY_DIR)

    set(CUTFEM_FOUND "YES")
    foreach( lib_name ${CUTFEM_COMPONENENT})
        find_library(COMP_LIBRARY
            NAMES 
                ${lib_name}
            PATHS 
                ${CUTFEM_LIBRARY_DIR}
                NO_DEFAULT_PATH
        )
        list(APPEND CUTFEM_LIBRARIES ${COMP_LIBRARY})
        SET(COMP_LIBRARY "COMP_LIBRARY-NOTFOUND")
    endforeach(lib_name)

    message(STATUS "CutFEM Found ")
    message(STATUS "libs are ${CUTFEM_LIBRARIES}")
else()
    message(FATAL_ERROR " CutFEM Not Found ")
endif()
message( STATUS "inlude directory: ${CUTFEM_INCLUDE_DIR}" )
message( STATUS "lib directory: ${CUTFEM_LIBRARY_DIR} ")