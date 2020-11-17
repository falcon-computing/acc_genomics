add_custom_target(IntelOCL)

if(DEFINED ENV{ALTERAOCLSDKROOT} OR DEFINED ENV{INTELFPGAOCLSDKROOT})
  set(IntelOCL_DIR "$ENV{ALTERAOCLSDKROOT}")
  execute_process(COMMAND aocl compile-config OUTPUT_VARIABLE IntelOCL_INCLUDE_DIRS)
  string(REPLACE "-I" "" IntelOCL_INCLUDE_DIRS ${IntelOCL_INCLUDE_DIRS})
  
  execute_process(COMMAND aocl ldflags OUTPUT_VARIABLE IntelOCL_Link_Dir)
  string(REPLACE "\n" " " IntelOCL_Link_Dir "${IntelOCL_Link_Dir}")
  string(REPLACE " " ";" IntelOCL_Link_Dir "${IntelOCL_Link_Dir}")
  foreach (subflag ${IntelOCL_Link_Dir})
    string(SUBSTRING ${subflag} 2 -1 iflag)
    list(APPEND IntelOCL_LIBRARY_DIRS ${iflag})
  endforeach()
  
  execute_process(COMMAND aocl ldlibs OUTPUT_VARIABLE IntelOCL_Link_Libs)
  string(REPLACE "\n" " " IntelOCL_Link_Libs ${IntelOCL_Link_Libs})
  string(REPLACE " " ";" IntelOCL_Link_Libs ${IntelOCL_Link_Libs})
  list(FILTER IntelOCL_Link_Libs INCLUDE REGEX "-l.[a-zA-Z0-9-_]+" )
  foreach (sublib ${IntelOCL_Link_Libs})
    string(SUBSTRING ${sublib} 2 -1 ilib)
    list(APPEND IntelOCL_LIBRARIES ${ilib})
  endforeach()

#  execute_process(COMMAND aocl compile-config OUTPUT_VARIABLE IntelOCL_INCLUDE_CONFIG)
#  execute_process(COMMAND python -c "import sys; print(';'.join([ token[2:] for token in sys.argv[1].split() if token[:2]==\"-I\" ])); " "${IntelOCL_INCLUDE_CONFIG}" OUTPUT_VARIABLE IntelOCL_INCLUDE_DIRS)
#  string(REPLACE "\n" "" IntelOCL_INCLUDE_DIRS "${IntelOCL_INCLUDE_DIRS}")
#  
#  execute_process(COMMAND aocl link-config OUTPUT_VARIABLE IntelOCL_LINK_CONFIG)
#  execute_process(COMMAND python -c "import sys; print(';'.join([ token[2:] for token in sys.argv[1].split() if token[:2]==\"-L\" ])); " "${IntelOCL_LINK_CONFIG}" OUTPUT_VARIABLE IntelOCL_LIBRARY_DIRS)
#  string(REPLACE "\n" "" IntelOCL_LIBRARY_DIRS "${IntelOCL_LIBRARY_DIRS}")
#  
#  execute_process(COMMAND python -c "import sys; print(';'.join([ token[2:] for token in sys.argv[1].split() if token[:2]==\"-l\" ])); " "${IntelOCL_LINK_CONFIG}" OUTPUT_VARIABLE IntelOCL_LIBRARIES)
#  string(REPLACE "\n" "" IntelOCL_LIBRARIES "${IntelOCL_LIBRARIES}")
#  
#  execute_process(COMMAND python -c "import sys; print(' '.join([ token for token in sys.argv[1].split() if token[:2]!=\"-L\" and token[:2]!=\"-l\" ])); " "${IntelOCL_LINK_CONFIG}" OUTPUT_VARIABLE IntelOCL_LINKER_FLAGS)
#  string(REPLACE "\n" "" IntelOCL_LINKER_FLAGS "${IntelOCL_LINKER_FLAGS}")
#
#  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${IntelOCL_LINKER_FLAGS}")
  set(IntelOCL_FOUND TRUE)

  message(STATUS "Found Intel OpenCL Environment")
endif()
