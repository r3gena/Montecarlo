cmake_minimum_required(VERSION 3.13)

project(Montecarlo)

include_directories(${PROJECT_SOURCE_DIR})
add_library(applib progress_bar.cpp potenciales.cpp Sistema_Montecarlo.cpp)
add_executable(app run.cpp)
target_link_libraries(app applib)
