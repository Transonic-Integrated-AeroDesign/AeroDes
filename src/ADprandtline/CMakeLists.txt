add_library(ad_pran STATIC ../ADprandtline.cpp)
target_sources(ad_pran PUBLIC ../ADprandtline.hpp)
target_include_directories(ad_pran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ad_pran vlib adlib)

add_library(ad_pran_obj OBJECT ../ADprandtline.cpp)
target_sources(ad_pran_obj PUBLIC ../ADprandtline.hpp)
target_include_directories(ad_pran_obj PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ad_pran_obj vlib adlib)
