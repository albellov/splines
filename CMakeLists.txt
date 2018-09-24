cmake_minimum_required(VERSION 3.10)
project(splines)

set(CMAKE_CXX_STANDARD 11)

add_executable(data_gen tools/data_gen.cpp)

add_executable(splines_2d spline_2d/src/main.cpp spline_2d/src/spline_2d.cpp tools/tools.cpp)
add_executable(splines_3d spline_3d/src/main.cpp spline_3d/src/spline_3d.cpp tools/tools.cpp)