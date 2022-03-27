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
        propulsion/cfmactu.cpp)

# assign target info after adding library
set_target_properties(aerolib PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(aerolib PROPERTIES SOVERSION 1)

set_target_properties(aerolib PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1
        PUBLIC_HEADER "aerodes.hpp;prandtline/prandtline.hpp;canareq/canareq.hpp")

#[[target_link_libraries(aerolib PUBLIC ${EXTRA_LIBS})]]

# configure a header file to pass the version number only
#configure_file(aerolib.pc.in aerolib.pc @ONLY)

target_include_directories(aerolib PRIVATE .)
#target_include_directories(aerolib PUBLIC_HEADER ${MyHeaders})
#target_include_directories(aerolib PRIVATE src)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
#target_include_directories(aerolib PUBLIC
#        "${PROJECT_BINARY_DIR}")

include_directories("${CMAKE_SOURCE_DIR}/src"
        "${CMAKE_SOURCE_DIR}/src/prandtline"
        "${CMAKE_SOURCE_DIR}/src/canareq")

#install(TARGETS aerolib DESTINATION bin
install(TARGETS aerolib
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
