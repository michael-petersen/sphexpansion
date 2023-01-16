cmake_minimum_required(VERSION 3.14)
project(
  SPHEXPANSION
  VERSION 0.2.2
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

# this may or may not be needed, depending on if you are using anaconda (needed)
# or if you are on an M1/2 Mac using native python installs (not needed)
#set(CMAKE_OSX_ARCHITECTURES "x86_64")

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# check which python this actually picks up!
find_package(Python COMPONENTS Interpreter Development)


# required compiler features
find_package(Eigen3 REQUIRED)
add_subdirectory(extern/pybind11)
add_subdirectory(extern/yaml-cpp)


if(pybind11_FOUND)
  set(HAVEPYBIND11 TRUE)
  message(STATUS "Do we have pybind11? ${HAVEPYBIND11}")
  add_compile_definitions(HAVEPYBIND11=${HAVEPYBIND11})
endif()

if(Eigen3_FOUND)
  set(HAVEEIGEN TRUE)
  message(STATUS "Do we have eigen? ${HAVEEIGEN}")
  add_compile_definitions(HAVEEIGEN=${HAVEEIGEN})
endif()

if(yaml-cpp_FOUND)
  set(HAVEYAML TRUE)
  message(STATUS "Do we have yaml-cpp? ${HAVEYAML}")
endif()

#configure_file(${CMAKE_SOURCE_DIR}/examples/config_flags.h.in ${CMAKE_SOURCE_DIR}/flags.h)

add_subdirectory(examples)
