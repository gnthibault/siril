cmake_minimum_required (VERSION 3.5)

#####################################
#  Adding binding code directories  #
#####################################

if (WITH_PYTHON)
  add_subdirectory(python)
endif ()
