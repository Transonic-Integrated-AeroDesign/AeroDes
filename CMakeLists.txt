cmake_minimum_required(VERSION 3.05)

# set the project name and version
project(aerolib VERSION 1.0.1 DESCRIPTION "open source aero design")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# should we use debug mode
option(DBG "debug mode" OFF)

if(DBG)
    set(DBG_VALUE 1)
endif()

configure_file(config.hpp.in "${CMAKE_SOURCE_DIR}/src/config.hpp" @ONLY)
configure_file(config.hpp.in "${CMAKE_SOURCE_DIR}/src/ADwake/config.hpp" @ONLY)
configure_file(config.hpp.in "${CMAKE_SOURCE_DIR}/src/ADcanareq/config.hpp" @ONLY)
configure_file(config.hpp.in "${CMAKE_SOURCE_DIR}/src/ADprandtline/config.hpp" @ONLY)

# add module libraries
add_subdirectory(src)

#
#include_directories(src)
#include_directories(src/ADprandtline)
#include_directories(src/ADwake)
#include_directories(src/ADcanareq)

# executables
#add_executable(ADcfmactu
#        src/ADoutput.cpp
#        src/ADvariables.hpp
#        src/ADmemory.hpp
#        src/propulsion/cfmactu.cpp
#        src/propulsion/cfmactu.hpp
#        src/propulsion/main.cpp
#        src/AD.cpp
#        src/AD.hpp
#        src/ADcanareq/ADcanareq.cpp
#        src/ADcanareq/ADcanareq.hpp
#        src/ADwake/ADwake.cpp
#        src/ADwake/ADwake.hpp
#        src/ADprandtline/ADprandtline.cpp
#        src/ADprandtline/ADprandtline.hpp src/ADerror.hpp)
#
#add_executable(ADprandtline
#        src/ADoutput.cpp
#        src/ADprandtline/main.cpp
#        src/ADprandtline/ADprandtline.cpp
#        src/ADprandtline/ADprandtline.hpp
#        src/ADvariables.hpp
#        src/ADmemory.hpp
#        src/AD.cpp
#        src/AD.hpp
#        src/ADwake/ADwake.cpp
#        src/ADwake/ADwake.hpp
#        src/ADcanareq/ADcanareq.cpp
#        src/ADcanareq/ADcanareq.hpp src/ADerror.hpp)
#
#add_executable(ADwake
#        src/ADoutput.cpp
#        src/ADvariables.hpp
#        src/ADmemory.hpp
#        src/ADwake/main.cpp
#        src/ADwake/ADwake.cpp
#        src/ADwake/ADwake.hpp
#        src/AD.cpp
#        src/AD.hpp
#        src/ADprandtline/ADprandtline.cpp
#        src/ADprandtline/ADprandtline.hpp
#        src/ADcanareq/ADcanareq.cpp
#        src/ADcanareq/ADcanareq.hpp src/ADerror.hpp)
#
#add_executable(ADcanary
#        src/AD.cpp
#        src/AD.hpp
#        src/ADoutput.cpp
#        src/ADprandtline/ADprandtline.cpp
#        src/ADvariables.hpp
#        src/ADmemory.hpp
#        src/ADwake/ADwake.cpp
#        src/ADwake/ADwake.hpp
#        src/ADcanareq/ADcanareq.cpp
#        src/ADcanareq/ADcanareq.hpp
#        src/ADcanareq/main.cpp src/ADerror.hpp)
#
#add_executable(ADtsd
#        src/AD.cpp
#        src/AD.hpp
#        src/ADoutput.cpp
#        src/ADvariables.hpp
#        src/ADmemory.hpp
#        src/ADprandtline/ADprandtline.cpp
#        src/ADwake/ADwake.cpp
#        src/ADcanareq/ADcanareq.cpp
#        src/tsd/tsd.cpp
#        src/tsd/tsd.hpp
#        src/tsd/main.cpp src/ADerror.hpp)
#
#install(TARGETS ADcfmactu)
#install(TARGETS ADprandtline)
#install(TARGETS ADwake)
#install(TARGETS ADcanary)
#install(TARGETS ADtsd)

#add_executable(smooth
#        src/AD.cpp
#        src/AD.hpp
#        src/ADvariables.hpp
#        src/smoothpolar/smoothpolar.cpp
#        src/smoothpolar/smoothpolar.hpp
#        src/smoothpolar/main.cpp)