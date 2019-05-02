cmake_minimum_required (VERSION 3.5)

#####################################
#  Find python include directories  #
#####################################
find_package(PythonLibs REQUIRED)
set(python_deps ${PYTHON_LIBRARIES})
include_directories($PYTHON_INCLUDE_DIR)

#####################################
#  Adding binding code directories  #
#####################################

include_directories("pybind11")
add_subdirectory(pybind11)
add_subdirectory(src)
