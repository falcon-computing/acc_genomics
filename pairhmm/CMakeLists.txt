include(ExternalProject)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

cmake_minimum_required(VERSION 2.8.12)

project(pairhmm_tb)

find_package(Boost 1.53.0 COMPONENTS
	  system thread iostreams filesystem regex REQUIRED)
find_package(Protobuf 2.5.0 REQUIRED)
find_package(Threads)

set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Google REQUIRED)
find_package(Blaze REQUIRED)
find_package(KSight REQUIRED)
find_package(IntelGKL REQUIRED)
find_package(XilinxOCL)
find_package(IntelOCL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DBOOST_NO_CXX11_SCOPED_ENUMS -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKERNEL_NAME=\\\"pmm_core_top0\\\"")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKERNEL_BANKID=0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DKERNEL_NUM_PE=64")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIMER_VERBOSE=0")
set(CMAKE_CXX_FLAGS_DEBUG "-O2 -g")

if (NO_PROFILE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_PROFILE=1")
endif()

include_directories(
    ${Google_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/interface
    ${CMAKE_SOURCE_DIR}/client
    ${CMAKE_SOURCE_DIR}/../common/include
    ${Blaze_INCLUDE_DIRS}
    ${Google_INCLUDE_DIRS}
    ${KSight_INCLUDE_DIRS}
    ${PMMClient_INCLUDE_DIRS}
    ${IntelGKL_INCLUDE_DIRS})

link_directories(
    ${Google_LIBRARY_DIRS}
    ${Blaze_LIBRARY_DIRS}
    ${KSight_LIBRARY_DIRS}
    ${IntelGKL_LIBRARY_DIRS})

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

add_subdirectory(interface)
add_subdirectory(client)
add_subdirectory(task)
add_subdirectory(host)
