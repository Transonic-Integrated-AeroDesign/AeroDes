add_library(ad_wake STATIC ../ADwake.cpp)
target_sources(ad_wake PUBLIC ../ADwake.hpp)
target_include_directories(ad_wake PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ad_wake vlib adlib)

add_library(ad_wake_obj OBJECT ../ADwake.cpp)
target_sources(ad_wake_obj PUBLIC ../ADwake.hpp)
target_include_directories(ad_wake_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ad_wake_obj vlib adlib)