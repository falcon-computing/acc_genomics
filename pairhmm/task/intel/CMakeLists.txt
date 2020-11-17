include_directories(
    ${IntelOCL_INCLUDE_DIRS}
    ${BLAZE_DIR}/platforms/altr_opencl/include)

link_directories(${IntelOCL_LIBRARY_DIRS})

add_library(PairHMMTask_intel SHARED PairHMMTask.cpp)
add_dependencies(PairHMMTask_intel Blaze)
add_dependencies(PairHMMTask_intel Google)
add_dependencies(PairHMMTask_intel KSight)

target_link_libraries(PairHMMTask_intel
  pmm_interface
  blaze
  blaze_message
  ${IntelOCL_LIBRARIES}
  ${KSight_LIBRARIES})

set_target_properties(PairHMMTask_intel
  PROPERTIES OUTPUT_NAME PairHMMTask)

install(TARGETS PairHMMTask_intel
  LIBRARY DESTINATION pairhmm-task/intel)
