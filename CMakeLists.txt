cmake_minimum_required(VERSION 3.10)
project(splines)

set(CMAKE_CXX_STANDARD 11)

add_executable(data_gen tools/data_gen.cpp)

add_executable(splines_2d univariate_spline/src/main.cpp univariate_spline/src/univariate_spline.cpp tools/tools.cpp)
add_executable(splines_3d bivariate_spline/src/main.cpp bivariate_spline/src/bivariate_spline.cpp tools/tools.cpp)