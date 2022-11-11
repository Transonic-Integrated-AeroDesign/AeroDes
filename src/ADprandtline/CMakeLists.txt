add_library(ad_pran STATIC ADprandtline.cpp)
target_sources(ad_pran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ADprandtline.hpp)
target_include_directories(ad_pran PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(ad_pran vlib adlib)