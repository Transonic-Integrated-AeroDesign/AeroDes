#[[
if(USE_PRANDTLINE)
    add_subdirectory(prandtline)
    list(APPEND EXTRA_LIBS prandtline)
endif()
]]

add_library(aerolib SHARED
        aerodes.cpp
        prandtline/prandtline.cpp
        canareq/canareq.cpp
        propulsion/cfmactu.cpp
        wake/wake.cpp)

# assign target info after adding library
set_target_properties(aerolib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(aerolib PROPERTIES SOVERSION 1)

set_target_properties(aerolib PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1
        PUBLIC_HEADER "aerodes.hpp;prandtline/prandtline.hpp;canareq/canareq.hpp;wake/wake.hpp")

# configure a header file to pass the version number only
target_include_directories(aerolib PRIVATE .)

include_directories("${CMAKE_SOURCE_DIR}/src"
        "${CMAKE_SOURCE_DIR}/src/prandtline"
        "${CMAKE_SOURCE_DIR}/src/canareq"
        "${CMAKE_SOURCE_DIR}/src/wake")

#install(TARGETS aerolib DESTINATION bin
install(TARGETS aerolib
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
