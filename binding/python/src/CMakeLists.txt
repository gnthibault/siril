cmake_minimum_required (VERSION 3.5)

##################
#  Generic Part  #
##################
project(pyPDWT)
set(OUTPUT_NAME pysiril)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#Compiler must support some option for pybind11 code optimization
if (UNIX)
  add_definitions(${CMAKE_CXX_FLAGS} "-flto -fvisibility=hidden -fopenmp -std=c++14")
endif ()

#########################
#  Python binding code  #
#########################
set(sources
  pysiril.cpp
)

set(headers
  BufferWrapper.h
)

####################
#  Lib Dependancy  #
####################

include_directories("${PRJ_ROOT_DIR}/src")

###############################
#  Equivalent of add_library  #
###############################
##pybind11_add_module(${OUTPUT_NAME} ${sources})
#add_library(${OUTPUT_NAME} MODULE ${sources})
#target_link_libraries(${OUTPUT_NAME} PRIVATE
#  ${openmp_deps} pybind11::module)
#set_target_properties(${OUTPUT_NAME} PROPERTIES 
#  PREFIX "${PYTHON_MODULE_PREFIX}")
#  #SUFFIX "${PYTHON_MODULE_EXTENSION}")
#install(TARGETS ${OUTPUT_NAME}
#  LIBRARY DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})
