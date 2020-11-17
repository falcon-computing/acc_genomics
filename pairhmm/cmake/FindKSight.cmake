ExternalProject_Add(ksight-download
    PREFIX "falconlm"
    URL https://s3.amazonaws.com/fcs-build-public/ksight.tgz
    CONFIGURE_COMMAND ""
    SOURCE_DIR "${CMAKE_BINARY_DIR}/ksight/install"
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(KSight)
add_dependencies(KSight ksight-download)

set(KSight_DIR "${CMAKE_BINARY_DIR}/ksight/install")
set(KSight_INCLUDE_DIRS "${KSight_DIR}/include")
set(KSight_LIBRARY_DIRS "${KSight_DIR}/lib")
set(KSight_LIBRARIES "ksight")
