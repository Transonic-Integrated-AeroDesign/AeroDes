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

# TARGET:
# * ACQUIRE SOURCE FILE NAMES
# * CREATE AERODES OBJECT
file(GLOB SOURCE "AD*.cpp")
add_library(ad_obj OBJECT ${SOURCE})
set_target_properties(ad_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)

# CREATE AERODES SHARED LIBRARY
if(SHARED)
    add_library(aerolib SHARED $<TARGET_OBJECTS:ad_obj>)
    install(TARGETS aerolib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endif()
