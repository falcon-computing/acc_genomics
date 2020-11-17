ExternalProject_Add(glog-download
    PREFIX "glog"
    URL https://s3.amazonaws.com/fcs-build-public/glog-falcon.tar.gz
    URL_MD5 2b1bb4285ef4c8963d5e0e338f1952b8
    SOURCE_DIR "${CMAKE_BINARY_DIR}/glog/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(gflags-download
    PREFIX "gflags"
    URL https://s3.amazonaws.com/fcs-build-public/gflags.tar.gz
    URL_MD5 1de8187489fbced5cc86c2ba241440e4
    SOURCE_DIR "${CMAKE_BINARY_DIR}/gflags/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(googletest-download
    PREFIX "googletest"
    URL https://s3.amazonaws.com/fcs-build-public/googletest.tar.gz
    URL_MD5 18fda945045354e264e3cca5428525d6
    SOURCE_DIR "${CMAKE_BINARY_DIR}/googletest/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(protobuf-download
    PREFIX "protobuf"
    URL https://s3.amazonaws.com/fcs-build-public/protobuf-2.5.0.tar.gz
    URL_MD5 a68c7ee81a65a84e57287af2a0738a75
    SOURCE_DIR "${CMAKE_BINARY_DIR}/protobuf/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(Google)
add_dependencies(Google
    glog-download
    gflags-download
    googletest-download
    protobuf-download)

set(Google_INCLUDE_DIRS
    ${CMAKE_BINARY_DIR}/glog/install/include
    ${CMAKE_BINARY_DIR}/gflags/install/include
    ${CMAKE_BINARY_DIR}/googletest/install/include
    ${CMAKE_BINARY_DIR}/gtest/install/include
    ${CMAKE_BINARY_DIR}/protobuf/install/include)

set(Google_LIBRARY_DIRS
    ${CMAKE_BINARY_DIR}/glog/install/lib
    ${CMAKE_BINARY_DIR}/gflags/install/lib
    ${CMAKE_BINARY_DIR}/googletest/install/lib
    ${CMAKE_BINARY_DIR}/protobuf/install/lib)

set(Google_LIBRARIES
    glog
    gflags
    protobuf
    gtest)
