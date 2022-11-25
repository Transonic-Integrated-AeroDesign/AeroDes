# ADD LIBRARY: vlib
add_library(vlib ADvariables.hpp)
set_target_properties(vlib PROPERTIES LINKER_LANGUAGE CXX)

# ADD LIBRARY: inlib
add_library(inlib ADinput.hpp)
set_target_properties(inlib PROPERTIES LINKER_LANGUAGE CXX)

# ADD LIBRARY: memlib
add_library(memlib ADmemory.hpp)
set_target_properties(memlib PROPERTIES LINKER_LANGUAGE CXX)

# ADD LIBRARY: errlib
add_library(errlib ADerror.hpp)
set_target_properties(errlib PROPERTIES LINKER_LANGUAGE CXX)

# ADD LIBRARY: ad_pran, ad_wake, ad_canary
add_subdirectory(ADprandtline)
add_subdirectory(ADwake)
add_subdirectory(ADcanareq)

# CREATE AERODES STATIC LIBRARY
add_library(adlib STATIC AD.cpp;ADoutput.cpp)
target_sources(adlib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/AD.hpp") # HEADER FILES
target_include_directories(adlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# CREATE AERODES OBJECT
add_library(ad_obj OBJECT AD.cpp;ADoutput.cpp;ADerror.hpp)
target_link_libraries(ad_obj vlib adlib errlib)

install(TARGETS adlib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ad_pran
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ad_wake
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ad_canary
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# CREATE AERODES SHARED LIBRARY
add_library(aerodes SHARED $<TARGET_OBJECTS:ad_obj> $<TARGET_OBJECTS:ad_pran_obj> $<TARGET_OBJECTS:ad_wake_obj> $<TARGET_OBJECTS:ad_canary_obj>)
install(TARGETS aerodes
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
