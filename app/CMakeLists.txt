cmake_minimum_required (VERSION 3.5)

#####################################
#  Adding various code directories  #
#####################################

add_subdirectory(siril)

if (CUDA_FOUND AND USE_CUDA)
#  add_subdirectory(sirilcuda)
endif ()
