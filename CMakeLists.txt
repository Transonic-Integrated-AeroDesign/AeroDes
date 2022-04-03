cmake_minimum_required(VERSION 3.17)

# set the project name and version
project(aerolib VERSION 1.0.0 DESCRIPTION "open source aero design")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# should we use debug mode
option(DBG_ENABLE "debug mode" OFF)

if(DBG_ENABLE)
    set(DBG_VALUE 1)
endif()

configure_file(config.hpp.in "${CMAKE_SOURCE_DIR}/src/wake/config.hpp" @ONLY)

# add module libraries
add_subdirectory(src)
#add_subdirectory(src/wake)

include_directories(src/prandtline)
include_directories(src/wake)

# executables
add_executable(prandtl
        src/prandtline/prandtline.cpp
        src/prandtline/prandtline.hpp
        src/prandtline/main.cpp)

add_executable(wake
        src/wake/wake.cpp
        src/wake/wake.hpp
        src/wake/main.cpp)

add_executable(canary
        src/prandtline/prandtline.cpp
        src/prandtline/prandtline.hpp
        src/wake/wake.cpp
        src/wake/wake.hpp
        src/canareq/canareq.cpp
        src/canareq/canareq.hpp
        src/canareq/main.cpp)