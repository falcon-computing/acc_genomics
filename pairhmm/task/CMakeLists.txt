if(XilinxOCL_FOUND)
  add_subdirectory(xlnx)
endif()

if(IntelOCL_FOUND)
  add_subdirectory(intel)
endif()
