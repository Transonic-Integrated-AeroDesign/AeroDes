cmake_minimum_required(VERSION 3.17)
project(AeroDes)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/canareq)
include_directories(src/propulsion)

add_executable(canareq
        src/canareq/canareq.cpp
        src/canareq/canareq.hpp
        src/canareq/main.cpp)

add_executable(cfmactu
        src/propulsion/cfmactu.cpp
        src/propulsion/cfmactu.hpp
        src/propulsion/main.cpp)

add_executable(prandtline
        src/prandtline/prandtline.cpp
        src/prandtline/prandtline.hpp
        src/prandtline/main.cpp)
