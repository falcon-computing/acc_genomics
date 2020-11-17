ExternalProject_Add(intelgkl-download
    PREFIX "intel-gkl"
    URL https://s3.amazonaws.com/fcs-build-public/gkl-v0.8.5.tgz
    URL_MD5 f7b85946b18d6efe5149ffa24566c5bf
    SOURCE_DIR "${CMAKE_BINARY_DIR}/intel-gkl/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(IntelGKL)
add_dependencies(IntelGKL intelgkl-download)

set(IntelGKL_DIR "${CMAKE_BINARY_DIR}/intel-gkl/install")
set(IntelGKL_INCLUDE_DIRS "${IntelGKL_DIR}/include")
set(IntelGKL_LIBRARY_DIRS "${IntelGKL_DIR}/lib")
set(IntelGKL_LIBRARIES "gkl-pairhmm")
