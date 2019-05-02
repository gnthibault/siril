cmake_minimum_required(VERSION 3.5)

project(test)

include_directories(../src/)

set(TEST_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TEST_BIN_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test)

#Add handmade macro to cmake path
list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(siriltest)

#Add various app tests
add_subdirectory(src)
