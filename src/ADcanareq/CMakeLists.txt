add_library(ad_canary STATIC ADcanareq.cpp)
target_sources(ad_canary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ADcanareq.hpp)
target_include_directories(ad_canary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ad_canary vlib adlib)