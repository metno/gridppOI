#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -Wall -Wextra -Wno-unused-parameter -fPIC")

find_package(Boost  REQUIRED)
find_package(gsl  REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(gridppOI ${SOURCES})

# Set up public headers for liba
set_target_properties(gridppOI PROPERTIES PUBLIC_HEADER "gridppOI.h")

# Create main and link to liba
# add_library(${PROJECT_NAME} gridppOI.cpp)
#target_link_libraries(${PROJECT_NAME} gridppOI)

# Set up public headers for swig-example
# set_target_properties(gridppOI PROPERTIES PUBLIC_HEADER "gridppOI.h")

# Add install target
Install(TARGETS ${PROJECT_NAME} gridppOI
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static
    PUBLIC_HEADER DESTINATION include)
add_subdirectory(swig)
