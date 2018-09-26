CC=g++

data_gen : ..\tools\data_gen.cpp
	$(CC) -o data_gen ..\tools\data_gen.cpp -std=c++11

univariate_spline : ..\tools\tools_2d.cpp src\univariate_spline.cpp src\main.cpp
	$(CC) -o univariate_spline src\main.cpp src\univariate_spline.cpp ..\tools\tools_2d.cpp -std=c++11

all : data_gen univariate_spline
	