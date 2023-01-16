cmake_minimum_required(VERSION 3.05)

# set the project name and version
project(aerolib VERSION 1.0.1 DESCRIPTION "open source aero design")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-fPIC")

# should we use debug mode
option(DBG "debug mode" OFF)
option(SHARED "compile shared libs" ON)

if(DBG)
    set(DBG_VALUE 1)
endif()

configure_file(config.hpp.in "${CMAKE_SOURCE_DIR}/src/config.hpp" @ONLY)

# main source code
add_subdirectory(src)
list(APPEND EXTRA_LIBS ad_obj)

# main aerodes executable
add_executable(aerodes src/run.cpp)
target_link_libraries(aerodes PUBLIC ${EXTRA_LIBS})