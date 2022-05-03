#[[
if(USE_PRANDTLINE)
    add_subdirectory(ADprandtline)
    list(APPEND EXTRA_LIBS ADprandtline)
endif()
]]

add_library(aerolib SHARED
        AD.cpp
        AD.hpp
        ADvariables.hpp
        ADprandtline/ADprandtline.cpp
        ADcanareq/ADcanareq.cpp
        propulsion/cfmactu.cpp
        ADwake/ADwake.cpp)

# assign target info after adding library
set_target_properties(aerolib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(aerolib PROPERTIES SOVERSION 1)

set_target_properties(aerolib PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1
        PUBLIC_HEADER "AD.hpp;ADvariables.hpp;ADprandtline/ADprandtline.hpp;ADcanareq/ADcanareq.hpp;ADwake/ADwake.hpp;ADwake/config.hpp")

# configure a header file to pass the version number only
target_include_directories(aerolib PRIVATE .)

include_directories("${CMAKE_SOURCE_DIR}/src"
        "${CMAKE_SOURCE_DIR}/src/ADprandtline"
        "${CMAKE_SOURCE_DIR}/src/ADcanareq"
        "${CMAKE_SOURCE_DIR}/src/ADwake")

#install(TARGETS aerolib DESTINATION bin
install(TARGETS aerolib
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
