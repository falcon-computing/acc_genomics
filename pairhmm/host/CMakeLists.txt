set(SRC_LIST main.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(HOST_MODE MATCHES "Local")
  if(XilinxOCL_FOUND)
    set(FPGA_INCLUDE_DIRS "${XilinxOCL_INCLUDE_DIRS}" )
    set(FPGA_LIBRARY_DIRS "${XilinxOCL_LIBRARY_DIRS}")
    set(FPGA_LIBRARIES "${XilinxOCL_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXILINX")

  elseif(IntelOCL_FOUND)
    set(FPGA_INCLUDE_DIRS "${IntelOCL_INCLUDE_DIRS}" )
    set(FPGA_LIBRARY_DIRS "${IntelOCL_LIBRARY_DIRS}")
    set(FPGA_LIBRARIES "${IntelOCL_LIBRARIES}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINTEL")

  else()
    message(SEND_ERROR "No OpenCL module loaded. Build without using FPGA." )
  endif()
 
  include_directories(${FPGA_INCLUDE_DIRS})
  link_directories(${FPGA_LIBRARY_DIRS})
  
  if(USE_BLAZE)
    message(STATUS "Local mode with Blaze client")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCLIENT_MODE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLOCAL_BLAZE")
  else()
    list(APPEND SRC_LIST PairHMMFpga.cpp)
  endif()
else()
  message(STATUS "Remote mode with Blaze client")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCLIENT_MODE")
endif()

add_executable(host_tb ${SRC_LIST})
add_dependencies(host_tb Google)
add_dependencies(host_tb KSight)
add_dependencies(host_tb Blaze)
add_dependencies(host_tb IntelGKL)
add_dependencies(host_tb pmm_interface)
add_dependencies(host_tb pmm_client)

target_link_libraries(host_tb
    pmm_client 
    pmm_interface
    ${IntelGKL_LIBRARIES}
    ${Blaze_LIBRARIES}
    ${KSight_LIBRARIES}
    ${Google_LIBRARIES} 
    ${Boost_LIBRARIES} 
    ${CMAKE_DL_LIBS} 
    ${FPGA_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS host_tb 
    RUNTIME DESTINATION ./)
