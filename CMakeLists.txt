cmake_minimum_required(VERSION 2.8)

# Change policy to allow a version to be set
cmake_policy(SET CMP0048 NEW)
project(gridppOI VERSION 0.1.0)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wno-unused-parameter -fPIC")
if (CMAKE_BUILD_TYPE MATCHES "DEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -pg -fprofile-arcs -ftest-coverage --coverage")
elseif(CMAKE_BUILD_TYPE MATCHES "RELEASE")
    # find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_BUILD_TYPE Debug)

# Add libgridppOI
include_directories(./include)
file(GLOB SOURCES "src/*.c*")
add_library(gridppOI SHARED ${SOURCES})
set_target_properties(gridppOI PROPERTIES PUBLIC_HEADER "include/gridppOI.h")
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Armadillo 6.5 REQUIRED)
target_link_libraries(gridppOI "${GSL_LIBRARIES}")
target_link_libraries(gridppOI "${ARMADILLO_LIBRARIES}")
target_link_libraries(gridppOI "${Boost_LIBRARIES}")
target_link_libraries(gridppOI "${BLAS_LIBRARIES}")
target_include_directories(gridppOI PUBLIC "${Boost_INCLUDE_DIRS}")
target_include_directories(gridppOI PUBLIC "${GSL_INCLUDE_DIRS}")
target_include_directories(gridppOI PUBLIC "${ALGLIB_INCLUDE_DIRS}")
target_include_directories(gridppOI PUBLIC "${ARMADILLO_INCLUDE_DIRS}")
install(TARGETS gridppOI
      EXPORT CMakeSwigTargets
      PUBLIC_HEADER DESTINATION include
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib
      )

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Add swig subdirectory
add_subdirectory(SWIG)

option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/config.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with
    # the application
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API
        documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be
    installed to generate the
    doxygen documentation")
endif (DOXYGEN_FOUND)

