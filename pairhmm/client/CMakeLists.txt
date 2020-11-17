file(GLOB SRC_LIST *.cpp)
file(GLOB HDR_LIST *.h)

add_library(pmm_client STATIC ${SRC_LIST})
add_dependencies(pmm_client pmm_interface)
add_dependencies(pmm_client Google)
add_dependencies(pmm_client Blaze)
add_dependencies(pmm_client KSight)
set_target_properties(pmm_client PROPERTIES PUBLIC_HEADER "${HDR_LIST}")

set(PMM_CLIENT_LIB pmm_client CACHE INTERNAL "client lib")

install(TARGETS pmm_client
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER DESTINATION include/pairhmm)
