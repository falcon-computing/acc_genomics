include_directories(
    ${XilinxOCL_INCLUDE_DIRS}
    ${BLAZE_DIR}/platforms/xlnx_opencl/include)

link_directories(${XilinxOCL_LIBRARY_DIRS})

set(proj_name PairHMMTask_xlnx)
add_library(${proj_name} SHARED PairHMMTask.cpp)
add_dependencies(${proj_name} Blaze)
add_dependencies(${proj_name} Google)
add_dependencies(${proj_name} KSight)

target_link_libraries(${proj_name}
  pmm_interface
  blaze
  blaze_message
  ${XilinxOCL_LIBRARIES})

set_target_properties(${proj_name}
  PROPERTIES OUTPUT_NAME PairHMMTask)

install(TARGETS PairHMMTask_xlnx
    LIBRARY DESTINATION pairhmm-task/xlnx)
