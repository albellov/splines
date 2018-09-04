cmake_minimum_required(VERSION 3.10)
project(spline_2d)

set(CMAKE_CXX_STANDARD 11)

add_executable(data_gen spline_2d/data_gen.cpp)
add_executable(spline spline_2d/main.cpp spline_2d/spline.cpp spline_2d/spline.h)