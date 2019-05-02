cmake_minimum_required(VERSION 3.5)

##################
#  Generic Part  #
##################
project(siril)
set(OUTPUT_NAME siril)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

###############
#  Src files  #
###############
set(sources
  siril.c
)

####################
#  Lib Dependancy  #
####################
include_directories("${PRJ_ROOT_DIR}/src")

################
#  Executable  #
################
add_executable(${OUTPUT_NAME} ${sources})
target_link_libraries(${OUTPUT_NAME}
  libsiril
  ${openmp_deps}
)
set_property(TARGET ${OUTPUT_NAME}
  PROPERTY OUTPUT_NAME ${OUTPUT_NAME})
install(TARGETS ${OUTPUT_NAME}
  RUNTIME DESTINATION ${RUNTIME_OUTPUT_DIRECTORY})

