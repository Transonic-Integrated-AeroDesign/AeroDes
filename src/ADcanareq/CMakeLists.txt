add_library(ad_canary STATIC ../ADcanareq.cpp)
target_sources(ad_canary PUBLIC ../ADcanareq.hpp)
target_include_directories(ad_canary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ad_canary vlib adlib)

add_library(ad_canary_obj OBJECT ../ADcanareq.cpp)
target_sources(ad_canary_obj PUBLIC ../ADcanareq.hpp)
target_include_directories(ad_canary_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ad_canary_obj vlib adlib)