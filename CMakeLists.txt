cmake_minimum_required(VERSION 3.17)
project(AeroDes)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(AeroDes
        src/canareq.cpp
        src/canareq.hpp
        src/main.cpp)
