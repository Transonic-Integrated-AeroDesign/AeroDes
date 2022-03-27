cmake_minimum_required(VERSION 3.17)

# set the project name and version
project(aerolib VERSION 1.0.0 DESCRIPTION "open source aero design")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# should we use our own math functions
option(USE_PRANDTLINE "Use prandtline" ON)

# add module libraries
add_subdirectory(src)

# add the whole library
#[[add_subdirectory(src)
add_subdirectory(src/prandtline)
add_subdirectory(src/canareq)
add_subdirectory(src/propulsion)]]

# does not function without pkgconfig
#configure_file(aerolib.pc.in aerolib.pc @ONLY)
#[[install(FILES ${CMAKE_BINARY_DIR}/aerolib.h
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)]]
