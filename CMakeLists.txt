cmake_minimum_required(VERSION 3.14)
project(
  SPHEXPANSION
  VERSION 0.1
  HOMEPAGE_URL https://github.com/michael-petersen/sphexpansion
  LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# compiler flags
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using intel
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qno-offload")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# required compiler features
find_package(Boost)
find_package(Eigen3)
find_package(yaml-cpp)

if(Eigen3_FOUND)
  set(HAVEEIGEN TRUE)
  message(STATUS "Do we have eigen? ${HAVEEIGEN}")
  add_compile_definitions(HAVEEIGEN=${HAVEEIGEN})
endif()

if(yaml-cpp_FOUND)
  set(HAVEYAML TRUE)
  message(STATUS "Do we have yaml-cpp? ${HAVEYAML}")
endif()

#configure_file(${CMAKE_SOURCE_DIR}/config_cmake.h.in ${CMAKE_SOURCE_DIR}/config.h)

add_subdirectory(examples)
