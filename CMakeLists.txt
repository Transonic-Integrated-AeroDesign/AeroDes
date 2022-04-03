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

include_directories(src/prandtline)
include_directories(src/wake)

# executables
add_executable(prandtl
        src/prandtline/prandtline.cpp
        src/prandtline/prandtline.hpp
        src/prandtline/main.cpp)

add_executable(canary
        src/canareq/canareq.cpp
        src/canareq/canareq.hpp
        src/prandtline/prandtline.cpp
        src/prandtline/prandtline.hpp
        src/canareq/main.cpp)

add_executable(wake
        src/wake/wake.cpp
        src/wake/wake.hpp
        src/wake/main.cpp)

