add_library(ad_wake STATIC ADwake.cpp)
target_sources(ad_wake PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ADwake.hpp)
target_include_directories(ad_wake PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ad_wake vlib adlib)