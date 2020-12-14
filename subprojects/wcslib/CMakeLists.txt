project(wcslib C CXX)
cmake_minimum_required(VERSION 3.10)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

find_package(CFITSIO REQUIRED)
include_directories(${CFITSIO_INCLUDE_DIR})

# add sources
file(GLOB SRC_FILES *.c)
#build a shared library
ADD_LIBRARY(wcslib STATIC ${SRC_FILES})
set_target_properties(wcslib PROPERTIES VERSION 7.3 SOVERSION 7)

# Glob headers
file(GLOB SRC_HEADERS *.h)
#install headers
INSTALL(FILES ${SRC_HEADERS} DESTINATION include/wcslib)

#install library
INSTALL(TARGETS wcslib RUNTIME DESTINATION bin
                       LIBRARY DESTINATION lib
                       ARCHIVE DESTINATION lib)
