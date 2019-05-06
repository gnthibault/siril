cmake_minimum_required(VERSION 3.5)

project(siril)
set(PROJECT_NAME siril)
add_definitions(-DPACKAGE="${PROJECT_NAME}")
add_definitions(-DVERSION="0.0")

#Default build destination
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(DOC_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc")
#Explicit make install destination
if (INSTALLDIR)
  set(RUNTIME_OUTPUT_DIRECTORY ${INSTALL_DIR}/bin)
  set(LIBRARY_OUTPUT_DIRECTORY ${INSTALL_DIR}/lib)
else ()
  set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
  set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
endif ()

add_definitions(-DPACKAGE_DATA_DIR="${LIBRARY_OUTPUT_DIRECTORY}")
add_definitions(-DPACKAGE_DOC_DIR="${DOC_INSTALL_DIR}")
add_definitions(-DLOCALEDIR="${PRJ_ROOT_DIR}/po")

#cmake ref
set(PRJ_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PRJ_CMAKE_DIR})

#local ref
set(PRJ_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#################################
#  Generic Compilation options  #
#################################

#Compiler must support c++14
if (CMAKE_COMPILER_IS_GNUCXX)
  set(CXX_FLAGS -std=c++14)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CXX_FLAGS ${CXX_FLAGS} -Wall -O0 -ggdb)
  else ()
    set(CXX_FLAGS ${CXX_FLAGS} -O3 -mtune=native -g -march=native )
    #valgrind --tool=callgrind ./program
    #kcachegrind callgrind.out.<pid>
    FIND_PACKAGE(OpenMP)
    if (OPENMP_FOUND)
      set(CXX_FLAGS ${CXX_FLAGS} ${OpenMP_CXX_FLAGS})
      set(openmp_deps gomp)
    endif ()
  endif ()
  add_compile_options(${CXX_FLAGS})
endif ()

##############################
#  Cuda Compilation options  #
##############################

find_package(CUDA 10.0)

if (CUDA_FOUND AND USE_CUDA)
  add_definitions(-DUSE_CUDA)
  
  #Eventually add some code instrumentation
  if (USE_NVTX)
    add_definitions(-DUSE_NVTX)
  endif ()

  #generic compilation flag, defines target architecture, type of build, ...
  list(APPEND CUDA_NVCC_FLAGS
    -gencode arch=compute_30,code=sm_30
    -gencode arch=compute_35,code=sm_35
    -gencode arch=compute_50,code=sm_50
    -gencode arch=compute_52,code=sm_52
    -gencode arch=compute_60,code=sm_60
    -gencode arch=compute_61,code=sm_61
    -std=c++11
    --expt-extended-lambda
    -rdc=true
    -Xcompiler -fPIC)

  set(CUDA_SEPARABLE_COMPILATION ON)

  if (CMAKE_BUILD_TYPE STREQUAL "Debug") 
    list(APPEND CUDA_NVCC_FLAGS 
      --device-debug
      --generate-line-info
      -g
      -G
      -DTHRUST_DEBUG )
  else ()
    list(APPEND CUDA_NVCC_FLAGS
      -use_fast_math
      -O3)
  endif ()
  #runtime libs for cuda
  list(APPEND cuda_deps cudadevrt nvToolsExt)
  link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib64)
endif ()

###############
# Boost part  #
###############
find_package(Boost COMPONENTS  program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    #runtime libs for boost
    set(boost_deps ${Boost_LIBRARIES})
endif ()
################
# /Boost part  #
################
###############
# Boost part  #
###############
find_package(Boost COMPONENTS  program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    #runtime libs for boost
    set(boost_deps ${Boost_LIBRARIES})
endif ()
################
# /Boost part  #
################
#############
# RAW part  #
#############
find_path(LIBRAW_INCLUDE_DIR libraw/libraw.h)
find_library(LIBRAW_LIBRARY raw)
include_directories(${LIBRAW_INCLUDE_DIR})
set(libraw_deps ${LIBRAW_LIBRARY})
##############
# /RAW part  #
##############
#################
# CFITSIO part  #
#################
find_path(CFITSIO_INCLUDE_DIR fitsio.h)
find_library(CFITSIO_LIBRARY cfitsio)
include_directories(${CFITSIO_INCLUDE_DIR})
set(cfitsio_deps ${CFITSIO_LIBRARY})
##################
# /CFITSIO part  #
##################
##############
# FFTW part  #
##############
find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARY fftw3)
include_directories(${FFTW_INCLUDE_DIR})
set(fftw_deps ${FFTW_LIBRARY})
###############
# /FFTW part  #
###############
#############
# GSL part  #
#############
find_package(GSL)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    #link_directories(${GSL_LIBRARY_DIRS})
    #runtime libs for gsl
    set(gsl_deps ${GSL_LIBRARIES})
endif ()
##############
# /GSL part  #
##############
###############
# Boost part  #
###############
find_package(Boost COMPONENTS  program_options)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    #runtime libs for boost
    set(boost_deps ${Boost_LIBRARIES})
endif ()
################
# /Boost part  #
################
#################
# Libtiff part  #
#################
find_package(TIFF)
if (TIFF_FOUND)
    include_directories(${TIFF_INCLUDE_DIRS})
    #link_directories(${TIFF_LIBRARY_DIRS})
    #runtime libs for tiff
    set(tiff_deps ${TIFF_LIBRARIES})
endif ()
##################
# /Libtiff part  #
##################
##############
# Jpeg part  #
##############
find_package(JPEG)
if (JPEG_FOUND)
    include_directories(${JPEG_INCLUDE_DIR})
    #link_directories(${JPEG_LIBRARY_DIRS})
    #runtime libs for jpeg
    set(jpeg_deps ${JPEG_LIBRARIES})
endif ()
###############
# /Jpeg part  #
###############
#############
# PNG part  #
#############
find_package(PNG)
if (PNG_FOUND)
    include_directories(${PNG_INCLUDE_DIRS})
    #link_directories(${PNG_LIBRARY_DIRS})
    #runtime libs for boost
    set(png_deps ${PNG_LIBRARIES})
endif ()
##############
# /PNG part  #
##############
################
# FFMPEG part  #
################
find_path(AVCODEC_INCLUDE_DIR libavcodec/avcodec.h)
find_library(AVCODEC_LIBRARY avcodec)
find_path(SWSCALE_INCLUDE_DIR libswscale/swscale.h)
find_library(SWSCALE_LIBRARY swscale)
find_path(SWRESAMPLE_INCLUDE_DIR libswresample/swresample.h)
find_library(SWRESAMPLE_LIBRARY swresample)
include_directories(${AVCODEC_INCLUDE_DIR} ${SWSCALE_INCLUDE_DIR} ${SWRESAMPLE_INCLUDE_DIR})
set(ffmpeg_deps ${AVCODEC_LIBRARY} ${SWSCALE_LIBRARY} ${SWRESAMPLE_LIBRARY})
#################
# /FFMPEG part  #
#################
###############
# FFMS2 part  #
###############
find_path(FFMS_INCLUDE_DIR ffms.h)
find_library(FFMS_LIBRARY ffms2)
include_directories(${FFMS_INCLUDE_DIR})
set(ffms_deps ${FFMS_LIBRARY})
################
# /FFMS2 part  #
################
##############
# CURL part  #
##############
find_package(CURL)
if (CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIRS})
    #link_directories(${CURL_LIBRARY_DIRS})
    set(curl_deps ${CURL_LIBRARIES})
endif ()
################
# /CURL part  #
################
################
# OpenCV part  #
################
find_package(OpenCV REQUIRED COMPONENTS core imgproc)
include_directories(${OpenCV_INCLUDE_DIRS})
set(opencv_deps ${OpenCV_LIBS})
#################
# /OpenCV part  #
#################
###################
# libconfig part  #
###################
find_path(LIBCONFIG_INCLUDE_DIR libconfig.h)
find_library(LIBCONFIG_LIBRARY config)
include_directories(${LIBCONFIG_INCLUDE_DIR})
set(libconfig_deps ${LIBCONFIG_LIBRARY})
####################
# /libconfig part  #
####################
#############
# GTK part  #
#############
# Use the package PkgConfig to detect GTK+ headers/library files
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTK3 REQUIRED gtk+-3.0)
include_directories(${GTK3_INCLUDE_DIRS})
#link_directories(${GTK3_LIBRARY_DIRS})
add_definitions(${GTK3_CFLAGS_OTHER})
#message("GTK DEPENDENCIES ARE ${GTK3_LIBRARIES}")
set(gtk_deps "pangoft2-1.0;cairo;${GTK3_LIBRARIES}")
##############
# /GTK part  #
##############
###############
# KPLOT part  #
###############
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/deps/kplot DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/deps)
add_custom_command(
    OUTPUT libkplot.a
    COMMAND make
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/deps/kplot
    COMMENT "Building libkplot"
    VERBATIM
)
add_custom_target(
  kplot ALL
  DEPENDS libkplot.a
)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/deps/kplot/libkplot.a
  DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
  PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ
)
include_directories("deps/kplot")
link_directories(${LIBRARY_OUTPUT_DIRECTORY})
#find_library(LIB_KPLOT libkplot.a
#  PATHS  ${CMAKE_CURRENT_BINARY_DIR}/deps/kplot/
#)
#set(kplot_deps ${LIB_KPLOT})
set(kplot_deps ${CMAKE_CURRENT_BINARY_DIR}/deps/kplot/libkplot.a m)

################
# /KPLOT part  #
################

#############
#  Testing  #
#############
if (TESTING_ENABLED)
  enable_testing()
endif ()

########################################
####        Doxygen Part            ####
########################################

# Add a target to generate API documentation with Doxygen
find_package(Doxygen)
if (DOXYGEN_FOUND AND USE_DOC)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/siril.doxyfile
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(
    doc
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DOC_INSTALL_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif ()

##############
#  Cppcheck  #
##############

# Add a target to generate Cppcheck report
include(${PRJ_CMAKE_DIR}/gen_cppcheck.cmake)
GENERATE_CPPCHECK(SOURCEDIRS src 
  ENABLE_IDS warning style performance portability information
  INCLUDES ${CURRENT_SOURCE_DIR}
  PROJECT_NAME ${PROJECT_NAME}
)

#####################################
#  Adding various code directories  #
#####################################

add_subdirectory(app)
add_subdirectory(binding)
add_subdirectory(deps)
add_subdirectory(src)

if (TESTING_ENABLED)
  add_subdirectory(test)
endif ()
