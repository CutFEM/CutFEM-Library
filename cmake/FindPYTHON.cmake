

# set(PYTHON_INCLUDES /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10/include/python3.10/
# /usr/local/lib/python3.10/site-packages/numpy/core/include)
# set(PYTHON_LIBRARIES /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10/lib/libpython3.10.dylib) 
# set(PYTHON_DIR /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10)


set(PYTHON_FOUND "NO")


find_path (PYTHON_INCLUDE_DIR
  NAMES Python.h
  PATHS
    /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10/include/python3.10/
    /opt/homebrew/opt/python@3.10/Frameworks/Python.framework/Versions/3.10/include/python3.10
    /usr/include/python3.8/
  )

  find_path(MATPLOTLIBCPP_DIR 
    NAMES
        matplotlibcpp.h
    PATHS
        /Users/thomasfrachon/lib/matplotlib-cpp
)

find_path (NUMPY_INCLUDE_DIR
  NAMES numpy/arrayobject.h
  PATHS
    /usr/local/lib/python3.10/site-packages/numpy/core/include
    /usr/lib/python3/dist-packages/numpy/core/include/
    /opt/homebrew/lib/python3.10/site-packages/numpy/core/include
)

find_path(PYTHON_LIBRARY_DIR
  NAMES 
        libpython3.10.dylib 
        libpython3.8.so
  PATHS
    /usr/lib/x86_64-linux-gnu/
    /usr/local/Cellar/python@3.10/3.10.9/Frameworks/Python.framework/Versions/3.10/lib/
    /opt/homebrew/opt/python@3.10/Frameworks/Python.framework/Versions/3.10/lib
)


if(PYTHON_INCLUDE_DIR AND NUMPY_INCLUDE_DIR AND PYTHON_LIBRARY_DIR AND MATPLOTLIBCPP_DIR)
  set(PYTHON_FOUND YES)

  find_library(PYTHON_LIBRARY
    NAMES 
        libpython3.10.dylib
        libpython3.8.so
    PATHS ${PYTHON_LIBRARY_DIR}
    NO_DEFAULT_PATH
  )


  set(PYTHON_LIBRARIES ${PYTHON_LIBRARY} )

  set(PYTHON_INCLUDES ${PYTHON_INCLUDE_DIR} ${NUMPY_INCLUDE_DIR} ${MATPLOTLIBCPP_DIR})

  message( "-- PYTHON_library FOUND")
  message( STATUS "PYTHON & NUMPY includes = ${PYTHON_INCLUDES}")
  message( STATUS "PYTHON library = ${PYTHON_LIBRARIES}")

  else()
    message( STATUS "PYTHON include dirs = ${PYHTON_INCLUDE_DIR}")
    message( STATUS "PYTHON library = ${PYTHON_LIBRARIES}")
endif()

