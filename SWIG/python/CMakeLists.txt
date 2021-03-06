find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(../${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
include_directories(../../include)


# Add swig module
#set(CMAKE_SWIG_FLAGS -interface _pygridppOI)
swig_add_library(gridppOI TYPE SHARED LANGUAGE python SOURCES ../${PROJECT_NAME}.i)
swig_link_libraries(gridppOI ${PROJECT_NAME} ${PYTHON_LIBRARIES})
#set_property(TARGET _pygridppOI PROPERTY OUTPUT_NAME gridppOI)

# Files to install with Python
set(PYTHON_INSTALL_FILES
        ${CMAKE_CURRENT_BINARY_DIR}/gridppOI.py
        ${CMAKE_CURRENT_BINARY_DIR}/_gridppOI.so)

# Configure setup.py and copy to output directory
set(SETUP_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/setup.in.py)
set(SETUP_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${SETUP_PY_IN} ${SETUP_PY_OUT})

#set(INIT_PY_IN ${CMAKE_CURRENT_SOURCE_DIR}/version.in.py)
#set(INIT_PY_OUT ${CMAKE_CURRENT_BINARY_DIR}/gridppOI/__init__.py)
#configure_file(${INIT_PY_IN} ${INIT_PY_OUT})

# Install target to call setup.py
add_custom_target(install-python
        DEPENDS _gridppOI
        COMMAND python3 ${SETUP_PY_OUT} install)

add_custom_target(build-python
        DEPENDS _gridppOI
        COMMAND python3 ${SETUP_PY_OUT} build)

add_custom_target(develop-python
        DEPENDS _gridppOI
        COMMAND python3 ${SETUP_PY_OUT} develop)

