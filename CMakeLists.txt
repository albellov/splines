cmake_minimum_required(VERSION 3.10)
project(splines)

set(CMAKE_CXX_STANDARD 11)

add_executable(data_gen tools/data_gen.cpp)

add_executable(univariate_spline univariate_spline/src/main.cpp univariate_spline/src/univariate_spline.cpp univariate_spline/src/univariate_spline.h tools/tools_2d.cpp tools/tools_2d.h)
add_executable(bivariate_spline bivariate_spline/src/main.cpp bivariate_spline/src/bivariate_spline.cpp bivariate_spline/src/bivariate_spline.h tools/tools_3d.cpp tools/tools_3d.h)