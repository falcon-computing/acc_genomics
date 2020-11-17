file(GLOB SRC_LIST *.cpp)
file(GLOB HDR_LIST *.h)

if(DEFINED DEVICE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFPGA_DEVICE=${DEVICE}")
  message(STATUS "Use FPGA Device: ${DEVICE}")
endif()

add_library(pmm_interface STATIC ${SRC_LIST})
add_dependencies(pmm_interface Google)
add_dependencies(pmm_interface KSight)
add_dependencies(pmm_interface Blaze)
set_target_properties(pmm_interface PROPERTIES PUBLIC_HEADER "${HDR_LIST}")

install(TARGETS pmm_interface
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/pairhmm)
